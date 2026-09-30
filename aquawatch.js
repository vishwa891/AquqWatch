/*****************************************************************************
 * PROJECT:     AquaWatch - Satellite-Based Surface Water Change Detection
 * PROBLEM:     Surface Water Change Detection - Problem 2.1
 * DOMAIN:      Water Resource Management
 * TECHNOLOGY:  Google Earth Engine | Sentinel-2 | NDWI | Geospatial Analysis |
 *              Interactive Dashboard
 *
 * PROBLEM STATEMENT
 *   Monitor the expansion and contraction of lakes, ponds and reservoirs by:
 *   (1) extracting water with NDWI/MNDWI, (2) comparing multiple years,
 *   (3) calculating water-area change, (4) identifying loss/gain hotspots.
 *
 * WORKFLOW
 *   Satellite data -> image processing -> geospatial analysis -> indicator
 *   -> interactive map/dashboard -> decision/recommendation
 *
 * DATASET
 *   COPERNICUS/S2_SR_HARMONIZED (Sentinel-2 Level-2A surface reflectance).
 *   Filtered by study area, year + season window, and scene cloud % (<60).
 *   Cloudy pixels (Scene Classification Layer classes 8, 9, 10 and saturated
 *   class 1) are masked. A MEDIAN composite is built per year so that
 *   remaining clouds/noise are suppressed.
 *
 * NDWI (McFeeters)   NDWI  = (Green - NIR)   / (Green + NIR)    = (B3 - B8)  / (B3 + B8)
 * MNDWI (Xu)         MNDWI = (Green - SWIR1) / (Green + SWIR1)  = (B3 - B11) / (B3 + B11)
 *   Higher positive values generally indicate open water. Vegetation and soil
 *   are usually near or below zero. MNDWI is often better in built-up areas.
 *
 * WATER THRESHOLD
 *   Pixel = water if index > threshold. Default NDWI 0.05 / MNDWI 0.00
 *   (adjustable in the UI). The SAME threshold is used for both years, so the
 *   comparison is consistent. Theoretical boundary is 0; a small positive
 *   value reduces false positives from wet soil and shadows.
 *
 * CHANGE DETECTION (post-classification comparison)
 *   change = waterHistorical + 2 * waterRecent
 *     0 = no significant water
 *     1 = historical water only  -> WATER LOSS
 *     2 = recent water only      -> WATER GAIN
 *     3 = water in both periods  -> STABLE WATER
 *   Only pixels with valid (cloud-free) data in BOTH years are compared, so
 *   clouds are never mistaken for water loss.
 *
 * AREA CALCULATION
 *   ee.Image.pixelArea() (m2) / 10000 = hectares per pixel. Each class mask
 *   is multiplied by pixel area and summed with reduceRegion().
 *   Net change = recent - historical (= gain - loss).
 *   Percentage change = net change / historical area * 100.
 *
 * HOTSPOTS
 *   The study area is divided into a grid (100-1000 m cells, chosen by area).
 *   A cell is a LOSS (or GAIN) hotspot when its lost (or gained) water area
 *   is at least 15% of the cell area.
 *
 * INTERPRETATION
 *   The recommendations are rule-based on the % change. They describe what
 *   the satellite data shows ("detected water-surface reduction"); they do NOT
 *   prove a cause (drought, encroachment, construction, dam operation...).
 *
 * LIMITATIONS
 *   - Same-season composites are used, but rainfall differs between years.
 *   - Spectral indices can confuse shadows, wet soil, algae or turbid water.
 *   - 10 m pixels cannot resolve very small ponds (<~0.02 ha).
 *   - Areas above ~300 km2 are analysed at 20-30 m to stay within EE limits.
 *   - No field verification / ground truth is used.
 *****************************************************************************/


// ===========================================================================
// 1. CONFIGURATION
// ===========================================================================
var CFG = {
  defaultYearA: '2019',
  defaultYearB: '2026',
  years: ['2018', '2019', '2020', '2021', '2022', '2023', '2024', '2025', '2026'],
  cloudMax: 60,                         // max scene cloud % (CLOUDY_PIXEL_PERCENTAGE)
  thrDefault: {NDWI: 0.05, MNDWI: 0.0},
  minAreaKm2: 0.05,                     // smaller than this -> rejected
  maxAreaKm2: 2000,                     // larger than this -> rejected
  minCoverage: 50,                      // % valid pixels needed for a good result
  hotspotFrac: 0.15,                    // cell is hotspot if change >= 15% of cell
  pctHighLoss: -20, pctLoss: -5,        // recommendation thresholds (% change)
  pctGain: 5, pctHighGain: 20
};

// Season windows: same months in both years => fair comparison.
var SEASONS = {
  'Jan-Apr (dry season, recommended)': {m0: 1, n: 4},
  'Jan-Dec (full year)': {m0: 1, n: 12},
  'Oct-Dec (NE monsoon)': {m0: 10, n: 3},
  'Jun-Sep (SW monsoon)': {m0: 6, n: 4}
};
var DEFAULT_SEASON = 'Jan-Apr (dry season, recommended)';

// Default demo area: Puzhal (Red Hills) & Sholavaram lakes, north Chennai.
var DEFAULT_AOI = ee.Geometry.Rectangle([80.12, 13.10, 80.27, 13.25]);

var COL = {
  loss: '#E53935', gain: '#43A047', stable: '#1E88E5',
  hist: '#0D47A1', rec: '#00ACC1', lossHot: 'FF9100', gainHot: 'C6FF00'
};


// ===========================================================================
// 2. IMAGE PROCESSING FUNCTIONS
// ===========================================================================
function maskClouds(img) {
  var scl = img.select('SCL');
  var bad = scl.eq(1).or(scl.eq(8)).or(scl.eq(9)).or(scl.eq(10));
  return img.select(['B2', 'B3', 'B4', 'B8', 'B11'])
    .divide(10000)
    .updateMask(bad.not())
    .copyProperties(img, ['system:time_start']);
}

function s2Collection(aoi, year, season) {
  var start = ee.Date.fromYMD(year, season.m0, 1);
  var end = start.advance(season.n, 'month');
  return ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(aoi)
    .filterDate(start, end)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CFG.cloudMax))
    .map(maskClouds);
}

function waterIndex(img, idx) {
  var bands = (idx === 'MNDWI') ? ['B3', 'B11'] : ['B3', 'B8'];
  return img.normalizedDifference(bands).rename('idx');
}

// Grid of square cells (size in metres) covering the AOI.
function makeGrid(aoi, cellM, lat) {
  var dy = cellM / 111320;
  var dx = cellM / (111320 * Math.cos(lat * Math.PI / 180));
  var ring = ee.List(aoi.bounds().coordinates().get(0));
  var ll = ee.List(ring.get(0));
  var ur = ee.List(ring.get(2));
  var xs = ee.List.sequence(ee.Number(ll.get(0)), ee.Number(ur.get(0)), dx);
  var ys = ee.List.sequence(ee.Number(ll.get(1)), ee.Number(ur.get(1)), dy);
  var cells = xs.map(function (x) {
    return ys.map(function (y) {
      var x0 = ee.Number(x);
      var y0 = ee.Number(y);
      return ee.Feature(ee.Geometry.Rectangle([x0, y0, x0.add(dx), y0.add(dy)]));
    });
  }).flatten();
  return ee.FeatureCollection(cells).filterBounds(aoi);
}

// Builds every Earth Engine object needed for the analysis (nothing runs yet).
function buildAnalysis(ctx) {
  var aoi = ctx.aoi;
  var compA = ctx.colA.median().clip(aoi);
  var compB = ctx.colB.median().clip(aoi);
  var iA = waterIndex(compA, ctx.idx);
  var iB = waterIndex(compB, ctx.idx);

  // Valid = cloud-free data available in BOTH years.
  var valid = iA.mask().and(iB.mask()).unmask(0);
  var wA = iA.gt(ctx.thr).unmask(0).and(valid);   // historical water 0/1
  var wB = iB.gt(ctx.thr).unmask(0).and(valid);   // recent water 0/1
  var lossM = wA.and(wB.not());
  var gainM = wB.and(wA.not());
  var stableM = wA.and(wB);
  var change = wA.add(wB.multiply(2)).rename('change');   // classes 0-3

  var px = ee.Image.pixelArea().divide(10000);            // hectares per pixel

  var stack = ee.Image.cat([
    wA.multiply(px).rename('hist'),
    wB.multiply(px).rename('rec'),
    lossM.multiply(px).rename('loss'),
    gainM.multiply(px).rename('gain'),
    stableM.multiply(px).rename('stable'),
    valid.multiply(px).rename('valid'),
    px.rename('total')
  ]);

  var stats = stack.reduceRegion({
    reducer: ee.Reducer.sum(), geometry: aoi, scale: ctx.scale,
    maxPixels: 1e10, bestEffort: true, tileScale: 4
  });

  // Hotspots on a grid
  var grid = makeGrid(aoi, ctx.cellM, ctx.lat);
  var cellStats = ee.Image.cat([
    lossM.multiply(px).rename('loss'),
    gainM.multiply(px).rename('gain')
  ]).reduceRegions({
    collection: grid, reducer: ee.Reducer.sum(), scale: ctx.scale, tileScale: 4
  });
  var lossHot = cellStats.filter(ee.Filter.gte('loss', ctx.minHa));
  var gainHot = cellStats.filter(ee.Filter.gte('gain', ctx.minHa));

  function topList(fc, prop) {
    return fc.sort(prop, false).limit(5).map(function (f) {
      var c = f.geometry().centroid(1).coordinates();
      return ee.Feature(null, {lon: c.get(0), lat: c.get(1), ha: f.get(prop)});
    }).toList(5);
  }
  var hot = ee.Dictionary({
    nLoss: lossHot.size(), nGain: gainHot.size(),
    topLoss: topList(lossHot, 'loss'), topGain: topList(gainHot, 'gain')
  });

  return {
    compA: compA, compB: compB, wA: wA, wB: wB, lossM: lossM, gainM: gainM,
    stableM: stableM, change: change, valid: valid,
    stats: stats, lossHot: lossHot, gainHot: gainHot, hot: hot
  };
}


// ===========================================================================
// 3. DECISION SUPPORT (rule-based, data-driven, no claims about causes)
// ===========================================================================
function recommend(s, ctx) {
  var level, text;
  if (s.hist <= 0 && s.rec <= 0) {
    return {level: 'warn', text: 'No surface water was detected in either year with the current index and threshold. ' +
      'Try lowering the threshold, switching to MNDWI, or selecting an area that contains water bodies.'};
  }
  if (s.hist <= 0) {
    return {level: 'gain', text: 'Observation: no water was detected in ' + ctx.yA + ', but ' + s.rec.toFixed(1) +
      ' ha was detected in ' + ctx.yB + '. Percentage change is undefined. Inspect the water-gain layer and verify with other data.'};
  }
  var gross = (s.gain + s.loss) / s.hist * 100;
  if (s.pct <= CFG.pctHighLoss) {
    level = 'high';
    text = 'Priority: detected water-surface reduction of ' + Math.abs(s.pct).toFixed(1) + '% (' + s.loss.toFixed(1) +
      ' ha lost). Investigate the water-loss hotspot areas and consider conservation or restoration measures. ' +
      'Confirm with rainfall records, field visits or higher-resolution imagery before deciding on causes.';
  } else if (s.pct <= CFG.pctLoss) {
    level = 'medium';
    text = 'Watch: a moderate water-surface reduction of ' + Math.abs(s.pct).toFixed(1) + '% was detected. ' +
      'Continue monitoring and review the loss hotspots. Check whether rainfall differed between the two years.';
  } else if (s.pct >= CFG.pctHighGain) {
    level = 'gain';
    text = 'Observation: surface-water extent has increased substantially (+' + s.pct.toFixed(1) + '%, ' + s.gain.toFixed(1) +
      ' ha gained). Check the gain hotspots; increases can reflect wetter conditions, new storage or seasonal timing.';
  } else if (s.pct >= CFG.pctGain) {
    level = 'gain';
    text = 'Observation: surface-water extent has increased (+' + s.pct.toFixed(1) + '%). Review the water-gain layer.';
  } else {
    level = 'stable';
    text = 'Observation: surface-water extent is relatively stable (' + s.pct.toFixed(1) + '% net change).';
    if (gross > 20) {
      text += ' However, gross change (gain + loss) is ' + gross.toFixed(0) + '% of the historical area, ' +
        'so water may have shifted location. Inspect both hotspot layers.';
    }
  }
  return {level: level, text: text};
}


// ===========================================================================
// 4. UI: HELPERS + STATE
// ===========================================================================
var useDrawn = false;
var overlays = [];
var lastCtx = null;

function num(x) { return (x === null || x === undefined || isNaN(x)) ? 0 : Number(x); }

var statusLabel = ui.Label('', {fontSize: '12px', margin: '6px 0'});
function setStatus(msg, kind) {
  var color = kind === 'error' ? '#C62828' : (kind === 'warn' ? '#E65100' :
             (kind === 'ok' ? '#2E7D32' : '#455A64'));
  statusLabel.setValue(msg);
  statusLabel.style().set('color', color);
}

function heading(text) {
  return ui.Label(text, {fontWeight: 'bold', fontSize: '15px', color: '#0D47A1',
    margin: '16px 0 4px 0', border: '0px 0px 1px 0px solid #B0BEC5'});
}
function note(text) {
  return ui.Label(text, {fontSize: '12px', color: '#546E7A', margin: '2px 0 4px 0', whiteSpace: 'pre-wrap'});
}
function hrow(widgets) { return ui.Panel(widgets, ui.Panel.Layout.flow('horizontal')); }

function makeCard(title, color) {
  var value = ui.Label('-', {fontSize: '19px', fontWeight: 'bold', color: color, margin: '2px 0 0 0'});
  var p = ui.Panel([ui.Label(title, {fontSize: '11px', color: '#607D8B', margin: '0'}), value],
    ui.Panel.Layout.flow('vertical'),
    {border: '1px solid #CFD8DC', padding: '6px 8px', margin: '3px', width: '188px', backgroundColor: '#FAFAFA'});
  return {panel: p, value: value};
}


// ===========================================================================
// 5. UI: MAP + DRAWING TOOLS
// ===========================================================================
var map = ui.Map();
map.setOptions('HYBRID');

var drawingTools = map.drawingTools();
drawingTools.setShown(false);
while (drawingTools.layers().length() > 0) {
  drawingTools.layers().remove(drawingTools.layers().get(0));
}
drawingTools.layers().add(ui.Map.GeometryLayer({
  geometries: null, name: 'Study area (drawn)', color: 'FFD600'
}));

function clearDrawn() {
  drawingTools.stop();
  var geoms = drawingTools.layers().get(0).geometries();
  while (geoms.length() > 0) { geoms.remove(geoms.get(0)); }
}

function startDraw(shape) {
  clearDrawn();
  useDrawn = false;
  drawingTools.setShape(shape);
  drawingTools.draw();
  setStatus('Click on the map to draw your study area. For a polygon, double-click to finish.', 'info');
}

drawingTools.onDraw(function () {
  useDrawn = true;
  aoiLabel.setValue('Study area: your drawn area');
  setStatus('Area captured. Choose years and press "Run Analysis".', 'ok');
});

function getAOI() {
  if (useDrawn) {
    var layer = drawingTools.layers().get(0);
    if (layer.geometries().length() > 0) { return layer.getEeObject(); }
  }
  return DEFAULT_AOI;
}

function outlineOf(aoi) {
  return ee.FeatureCollection([ee.Feature(aoi)])
    .style({color: 'FFEB3B', fillColor: '00000000', width: 2});
}

function showIdleMap() {
  map.layers().reset();
  map.addLayer(outlineOf(DEFAULT_AOI), {}, 'Default study area');
  map.centerObject(DEFAULT_AOI, 12);
}


// ===========================================================================
// 6. UI: LEFT PANEL
// ===========================================================================
var panel = ui.Panel({style: {width: '430px', padding: '10px'}});

panel.add(ui.Label('AquaWatch', {fontSize: '30px', fontWeight: 'bold', color: '#0D47A1', margin: '0'}));
panel.add(ui.Label('Satellite-Based Surface Water Change Detection',
  {fontSize: '15px', fontWeight: 'bold', color: '#37474F', margin: '2px 0'}));
panel.add(ui.Label('Monitoring water-body expansion and contraction using Sentinel-2 satellite imagery',
  {fontSize: '12px', color: '#607D8B', margin: '0 0 6px 0'}));

// --- Instructions
panel.add(heading('How to use'));
panel.add(note('1. Use the default area or draw your own.\n' +
  '2. Pick a Historical year and a Recent year.\n' +
  '3. Press "Run Analysis" and wait (10-60 s).\n' +
  '4. Read the map, statistics, hotspots and recommendation.'));

// --- Study area
panel.add(heading('1. Study area'));
var aoiLabel = ui.Label('Study area: default demo area (Puzhal & Sholavaram lakes, Chennai region, Tamil Nadu)',
  {fontSize: '12px', margin: '2px 0 4px 0'});
panel.add(aoiLabel);
panel.add(hrow([
  ui.Button({label: 'Use default area', onClick: function () {
    clearDrawn(); useDrawn = false;
    aoiLabel.setValue('Study area: default demo area (Puzhal & Sholavaram lakes, Chennai region, Tamil Nadu)');
    showIdleMap(); setStatus('Default area selected.', 'info');
  }}),
  ui.Button({label: 'Draw rectangle', onClick: function () { startDraw('rectangle'); }}),
  ui.Button({label: 'Draw polygon', onClick: function () { startDraw('polygon'); }})
]));

// --- Settings
panel.add(heading('2. Analysis settings'));
var yearASel = ui.Select({items: CFG.years, value: CFG.defaultYearA, style: {stretch: 'horizontal'}});
var yearBSel = ui.Select({items: CFG.years, value: CFG.defaultYearB, style: {stretch: 'horizontal'}});
panel.add(hrow([
  ui.Panel([ui.Label('Historical year', {fontSize: '12px'}), yearASel], null, {width: '195px'}),
  ui.Panel([ui.Label('Recent year', {fontSize: '12px'}), yearBSel], null, {width: '195px'})
]));

var seasonSel = ui.Select({items: Object.keys(SEASONS), value: DEFAULT_SEASON, style: {stretch: 'horizontal'}});
panel.add(ui.Label('Season window (same months used for both years)', {fontSize: '12px'}));
panel.add(seasonSel);

var thrSlider = ui.Slider({min: -0.3, max: 0.6, value: CFG.thrDefault.NDWI, step: 0.01,
  style: {stretch: 'horizontal'}});
var indexSel = ui.Select({items: ['NDWI', 'MNDWI'], value: 'NDWI', onChange: function (v) {
  thrSlider.setValue(CFG.thrDefault[v]);
}});
panel.add(ui.Label('Water index', {fontSize: '12px'}));
panel.add(indexSel);
panel.add(ui.Label('Water threshold (index > threshold = water)', {fontSize: '12px'}));
panel.add(thrSlider);
panel.add(note('NDWI = (Green - NIR) / (Green + NIR) = (B3 - B8) / (B3 + B8). ' +
  'Higher positive NDWI values generally indicate water. MNDWI = (B3 - B11) / (B3 + B11) ' +
  'often works better near buildings.'));

var opacitySlider = ui.Slider({min: 0.1, max: 1, value: 0.8, step: 0.05, style: {stretch: 'horizontal'},
  onChange: function (v) { overlays.forEach(function (l) { l.setOpacity(v); }); }});
panel.add(ui.Label('Layer opacity', {fontSize: '12px'}));
panel.add(opacitySlider);

panel.add(hrow([
  ui.Button({label: 'Run Analysis', onClick: function () { runAnalysis(); }}),
  ui.Button({label: 'Reset', onClick: function () { resetAll(); }})
]));
panel.add(statusLabel);

// --- Statistics
panel.add(heading('3. Water-area statistics'));
var cHist = makeCard('Historical water area', COL.hist);
var cRec = makeCard('Recent water area', COL.rec);
var cGain = makeCard('Water gain', COL.gain);
var cLoss = makeCard('Water loss', COL.loss);
var cNet = makeCard('Net change', '#37474F');
var cPct = makeCard('Percentage change', '#37474F');
panel.add(hrow([cHist.panel, cRec.panel]));
panel.add(hrow([cGain.panel, cLoss.panel]));
panel.add(hrow([cNet.panel, cPct.panel]));
var metaLabel = note('');
panel.add(metaLabel);

// --- Chart
panel.add(heading('4. Chart'));
var chartPanel = ui.Panel();
panel.add(chartPanel);

// --- Recommendation
panel.add(heading('5. Decision support'));
var recPanel = ui.Panel([], null, {padding: '8px', margin: '2px 0'});
panel.add(recPanel);

// --- Hotspots
panel.add(heading('6. Hotspots (grid cells)'));
var hotPanel = ui.Panel();
panel.add(hotPanel);

// --- Time series
panel.add(heading('7. Optional: yearly water trend'));
panel.add(note('Water area for every year 2018-2026 in the selected season (20-30 m resolution, approximate). ' +
  'Years with under 50% valid cloud-free coverage are left out.'));
panel.add(ui.Button({label: 'Compute yearly trend', onClick: function () { runTimeSeries(); }}));
var tsPanel = ui.Panel();
panel.add(tsPanel);

// --- Limitations
panel.add(heading('Limitations'));
panel.add(note(
  '- Results show detected water-surface change only; they do not prove its cause.\n' +
  '- Rainfall, season and reservoir operation differ between years.\n' +
  '- Shadows, wet soil, algae or turbid water can be mis-classified by NDWI/MNDWI.\n' +
  '- 10 m pixels miss very small ponds; large areas are analysed at 20-30 m.\n' +
  '- Verify important findings with field data and other imagery.'));
panel.add(note('Data: Copernicus Sentinel-2 (COPERNICUS/S2_SR_HARMONIZED) via Google Earth Engine.'));

// --- Legend (on the map)
function legendRow(color, text, outline) {
  var box = ui.Label('', {backgroundColor: color, padding: '8px', margin: '0 8px 4px 0',
    border: '2px solid ' + (outline || '#9E9E9E')});
  return ui.Panel([box, ui.Label(text, {fontSize: '12px', margin: '0 0 4px 0'})],
    ui.Panel.Layout.flow('horizontal'));
}
var legend = ui.Panel({style: {position: 'bottom-left', padding: '8px 12px', backgroundColor: 'white'}});
legend.add(ui.Label('Water change classes', {fontWeight: 'bold', margin: '0 0 4px 0'}));
legend.add(legendRow(COL.loss, '1 - Water loss (water in historical year only)'));
legend.add(legendRow(COL.gain, '2 - Water gain (water in recent year only)'));
legend.add(legendRow(COL.stable, '3 - Stable water (both years)'));
legend.add(legendRow('#FFFFFF', '0 - No significant water (transparent)'));
legend.add(legendRow('#FF910040', 'Loss hotspot cell', '#FF9100'));
legend.add(legendRow('#C6FF0040', 'Gain hotspot cell', '#C6FF00'));
map.add(legend);


// ===========================================================================
// 7. ANALYSIS FLOW
// ===========================================================================
var runBtnBusy = false;
function setBusy(b) { runBtnBusy = b; }

function fail(msg) {
  setBusy(false);
  setStatus(msg, 'error');
}

function clearResults() {
  [cHist, cRec, cGain, cLoss, cNet, cPct].forEach(function (c) { c.value.setValue('-'); });
  metaLabel.setValue('');
  chartPanel.clear();
  recPanel.clear();
  hotPanel.clear();
}

function runAnalysis() {
  if (runBtnBusy) { setStatus('An analysis is already running. Please wait.', 'warn'); return; }
  var yA = parseInt(yearASel.getValue(), 10);
  var yB = parseInt(yearBSel.getValue(), 10);
  var idx = indexSel.getValue();
  var thr = Number(thrSlider.getValue());
  var seasonName = seasonSel.getValue();
  var season = SEASONS[seasonName];

  if (isNaN(yA) || isNaN(yB)) { setStatus('Please select both years.', 'error'); return; }
  if (yA >= yB) { setStatus('Invalid years: the Historical year must be earlier than the Recent year.', 'error'); return; }
  if (isNaN(thr)) { setStatus('Invalid threshold value.', 'error'); return; }

  var aoi = getAOI();
  setBusy(true);
  clearResults();
  setStatus('Step 1/3: checking study area and Sentinel-2 imagery...', 'info');

  var colA = s2Collection(aoi, yA, season);
  var colB = s2Collection(aoi, yB, season);

  ee.Dictionary({
    nA: colA.size(),
    nB: colB.size(),
    areaKm2: aoi.area(1).divide(1e6),
    lat: aoi.bounds().centroid(1).coordinates().get(1)
  }).evaluate(function (v, err) {
    if (err || !v) { fail('Could not read the study area or imagery. ' + (err || '')); return; }
    var area = num(v.areaKm2);
    if (area < CFG.minAreaKm2) {
      fail('Study area is too small (' + area.toFixed(3) + ' km2). Draw an area of at least ' +
        CFG.minAreaKm2 + ' km2 (a polygon, not a point or line).');
      return;
    }
    if (area > CFG.maxAreaKm2) {
      fail('Study area is too large (' + area.toFixed(0) + ' km2). Please draw an area below ' +
        CFG.maxAreaKm2 + ' km2.');
      return;
    }
    if (num(v.nA) === 0) {
      fail('No usable Sentinel-2 images for ' + yA + ' (' + seasonName + ') in this area with cloud < ' +
        CFG.cloudMax + '%. Try another year, a longer season window, or a different area.');
      return;
    }
    if (num(v.nB) === 0) {
      fail('No usable Sentinel-2 images for ' + yB + ' (' + seasonName + ') in this area with cloud < ' +
        CFG.cloudMax + '%. Try another year (the current year may not have this season yet) or a longer season.');
      return;
    }

    var scale = area > 800 ? 30 : (area > 300 ? 20 : 10);
    var cellM = area < 5 ? 100 : (area < 300 ? 250 : (area < 1000 ? 500 : 1000));
    var ctx = {
      aoi: aoi, yA: yA, yB: yB, idx: idx, thr: thr, season: season, seasonName: seasonName,
      colA: colA, colB: colB, nA: num(v.nA), nB: num(v.nB), areaKm2: area, lat: num(v.lat),
      scale: scale, cellM: cellM, minHa: (cellM * cellM / 10000) * CFG.hotspotFrac
    };
    lastCtx = ctx;

    var a = buildAnalysis(ctx);
    showLayers(ctx, a);
    map.centerObject(aoi);

    var warn = '';
    if (area < 1) { warn += ' Note: very small area - few pixels, interpret with care.'; }
    if (ctx.nA < 3 || ctx.nB < 3) { warn += ' Note: few images used in a composite (' + ctx.nA + ' / ' + ctx.nB + ').'; }
    setStatus('Step 2/3: computing water-area statistics at ' + scale + ' m (10-60 s)...' + warn, 'info');

    a.stats.evaluate(function (s0, err2) {
      if (err2 || !s0) { fail('Statistics could not be computed. ' + (err2 || 'Try a smaller area.')); return; }
      showStats(ctx, s0);
      setStatus('Step 3/3: finding hotspots...', 'info');
      a.hot.evaluate(function (h, err3) {
        setBusy(false);
        if (err3 || !h) {
          setStatus('Statistics done, but hotspots could not be computed. ' + (err3 || ''), 'warn');
          return;
        }
        showHotspots(ctx, h);
        setStatus('Done. Analysed ' + ctx.yA + ' vs ' + ctx.yB + ' using ' + ctx.idx +
          ' (threshold ' + ctx.thr.toFixed(2) + ').' + warn, 'ok');
      });
    });
  });
}


// ===========================================================================
// 8. SHOW RESULTS
// ===========================================================================
function showLayers(ctx, a) {
  map.layers().reset();
  overlays = [];
  var op = opacitySlider.getValue();
  var tc = {bands: ['B4', 'B3', 'B2'], min: 0.02, max: 0.3};

  function addL(img, vis, name, shown, overlay) {
    var l = ui.Map.Layer(img, vis, name, shown, overlay ? op : 1);
    map.layers().add(l);
    if (overlay) { overlays.push(l); }
  }
  addL(a.compA, tc, 'Satellite (true colour) ' + ctx.yA, false, false);
  addL(a.compB, tc, 'Satellite (true colour) ' + ctx.yB, false, false);
  addL(a.wA.selfMask(), {min: 0, max: 1, palette: [COL.hist]}, 'Historical water ' + ctx.yA, false, true);
  addL(a.wB.selfMask(), {min: 0, max: 1, palette: [COL.rec]}, 'Recent water ' + ctx.yB, false, true);
  addL(a.stableM.selfMask(), {min: 0, max: 1, palette: [COL.stable]}, 'Stable water', false, true);
  addL(a.lossM.selfMask(), {min: 0, max: 1, palette: [COL.loss]}, 'Water loss', false, true);
  addL(a.gainM.selfMask(), {min: 0, max: 1, palette: [COL.gain]}, 'Water gain', false, true);
  addL(a.change.updateMask(a.valid).updateMask(a.change.gt(0)),
    {min: 1, max: 3, palette: [COL.loss, COL.gain, COL.stable]}, 'Water change map (classes 1-3)', true, true);
  addL(a.lossHot.style({color: COL.lossHot, fillColor: COL.lossHot + '40', width: 2}), {},
    'Water-loss hotspots', true, false);
  addL(a.gainHot.style({color: COL.gainHot, fillColor: COL.gainHot + '40', width: 2}), {},
    'Water-gain hotspots', true, false);
  addL(outlineOf(ctx.aoi), {}, 'Study area', true, false);
}

function showStats(ctx, s0) {
  var s = {
    hist: num(s0.hist), rec: num(s0.rec), loss: num(s0.loss), gain: num(s0.gain),
    stable: num(s0.stable), valid: num(s0.valid), total: num(s0.total)
  };
  s.net = s.rec - s.hist;
  s.pct = s.hist > 0 ? (s.net / s.hist) * 100 : null;
  s.cov = s.total > 0 ? (s.valid / s.total) * 100 : 0;

  cHist.panel.widgets().get(0).setValue('Historical water area (' + ctx.yA + ')');
  cRec.panel.widgets().get(0).setValue('Recent water area (' + ctx.yB + ')');
  cHist.value.setValue(s.hist.toFixed(2) + ' ha');
  cRec.value.setValue(s.rec.toFixed(2) + ' ha');
  cGain.value.setValue('+' + s.gain.toFixed(2) + ' ha');
  cLoss.value.setValue('-' + s.loss.toFixed(2) + ' ha');
  cNet.value.setValue((s.net >= 0 ? '+' : '') + s.net.toFixed(2) + ' ha');
  cPct.value.setValue(s.pct === null ? 'N/A' : ((s.pct >= 0 ? '+' : '') + s.pct.toFixed(1) + '%'));
  metaLabel.setValue('Study area: ' + ctx.areaKm2.toFixed(1) + ' km2 | Valid cloud-free coverage (both years): ' +
    s.cov.toFixed(0) + '% | Images used: ' + ctx.nA + ' (' + ctx.yA + '), ' + ctx.nB + ' (' + ctx.yB +
    ') | Resolution: ' + ctx.scale + ' m | Season: ' + ctx.seasonName);

  // Chart: water area summary
  var rows = [
    ['Category', 'Area (ha)', {role: 'style'}],
    ['Historical ' + ctx.yA, s.hist, 'color:' + COL.hist],
    ['Recent ' + ctx.yB, s.rec, 'color:' + COL.rec],
    ['Stable', s.stable, 'color:' + COL.stable],
    ['Gain', s.gain, 'color:' + COL.gain],
    ['Loss', s.loss, 'color:' + COL.loss]
  ];
  var chart = ui.Chart(rows, 'ColumnChart', {
    title: 'Water area summary (hectares)', legend: {position: 'none'},
    vAxis: {title: 'Hectares', minValue: 0}, hAxis: {slantedText: true}
  });
  chart.style().set({height: '260px', stretch: 'horizontal'});
  chartPanel.clear();
  chartPanel.add(chart);

  // Recommendation
  var rec = recommend(s, ctx);
  var text = rec.text;
  if (s.cov < CFG.minCoverage) {
    text += '\n\nData-quality warning: only ' + s.cov.toFixed(0) + '% of the area has cloud-free data in both years. ' +
      'Results may be unreliable. Try a longer season window or different years.';
  }
  var style = {
    high: ['#FFEBEE', '#E53935'], medium: ['#FFF3E0', '#FB8C00'], gain: ['#E3F2FD', '#1E88E5'],
    stable: ['#ECEFF1', '#78909C'], warn: ['#FFFDE7', '#FBC02D']
  }[rec.level];
  recPanel.clear();
  recPanel.style().set({backgroundColor: style[0], border: '2px solid ' + style[1]});
  recPanel.add(ui.Label(text, {fontSize: '13px', whiteSpace: 'pre-wrap', backgroundColor: style[0], margin: '2px'}));
  recPanel.add(ui.Label('Automatic, rule-based interpretation of detected change - not proof of cause.',
    {fontSize: '10px', color: '#78909C', backgroundColor: style[0], margin: '2px'}));
}

function showHotspots(ctx, h) {
  hotPanel.clear();
  hotPanel.add(note('Grid cell size: ' + ctx.cellM + ' m. A cell is a hotspot if the changed water area is at least ' +
    ctx.minHa.toFixed(2) + ' ha (' + (CFG.hotspotFrac * 100) + '% of the cell).'));
  hotPanel.add(ui.Label('Loss hotspot cells: ' + num(h.nLoss) + '   |   Gain hotspot cells: ' + num(h.nGain),
    {fontSize: '13px', fontWeight: 'bold'}));

  function rows(title, list, color) {
    hotPanel.add(ui.Label(title, {fontSize: '12px', fontWeight: 'bold', color: color, margin: '8px 0 2px 0'}));
    if (!list || list.length === 0) { hotPanel.add(note('None found.')); return; }
    list.forEach(function (f, i) {
      var p = f.properties;
      var lat = num(p.lat), lon = num(p.lon);
      hotPanel.add(hrow([
        ui.Label((i + 1) + '. ' + lat.toFixed(4) + ', ' + lon.toFixed(4) + '  -  ' + num(p.ha).toFixed(2) + ' ha',
          {fontSize: '12px', margin: '6px 8px 0 0'}),
        ui.Button({label: 'Zoom', onClick: function () { map.setCenter(lon, lat, 16); }})
      ]));
    });
  }
  rows('Top detected water-loss hotspots (lat, lon)', h.topLoss, COL.loss);
  rows('Top detected water-gain hotspots (lat, lon)', h.topGain, COL.gain);
}


// ===========================================================================
// 9. OPTIONAL: YEARLY TIME SERIES
// ===========================================================================
function runTimeSeries() {
  if (!lastCtx) { setStatus('Run the analysis first, then compute the yearly trend.', 'warn'); return; }
  var ctx = lastCtx;
  var tsScale = ctx.scale < 20 ? 20 : ctx.scale;
  var px = ee.Image.pixelArea().divide(10000);
  tsPanel.clear();
  tsPanel.add(ui.Label('Computing yearly trend... (can take up to a minute)', {fontSize: '12px', color: '#546E7A'}));

  var years = ee.List(CFG.years.map(function (y) { return parseInt(y, 10); }));
  var fc = ee.FeatureCollection(years.map(function (y) {
    y = ee.Number(y);
    var col = s2Collection(ctx.aoi, y, ctx.season);
    var label = y.format('%d');
    var empty = ee.Feature(null, {label: label, ha: 0, cov: 0});
    var good = ee.Feature(null, {}).set('label', label);
    var idxImg = waterIndex(col.median(), ctx.idx);
    var v = idxImg.mask().unmask(0);
    var w = idxImg.gt(ctx.thr).unmask(0).and(v);
    var st = ee.Image.cat([
      w.multiply(px).rename('ha'), v.multiply(px).rename('valid'), px.rename('total')
    ]).reduceRegion({
      reducer: ee.Reducer.sum(), geometry: ctx.aoi, scale: tsScale,
      maxPixels: 1e10, bestEffort: true, tileScale: 4
    });
    var full = good.set({
      ha: st.get('ha'),
      cov: ee.Number(st.get('valid')).divide(ee.Number(st.get('total'))).multiply(100)
    });
    return ee.Feature(ee.Algorithms.If(col.size().gt(0), full, empty));
  })).filter(ee.Filter.gte('cov', CFG.minCoverage));

  var chart = ui.Chart.feature.byFeature(fc, 'label', ['ha'])
    .setChartType('ColumnChart')
    .setOptions({
      title: 'Detected water area per year (' + ctx.idx + ' > ' + ctx.thr.toFixed(2) + ', ' + ctx.seasonName + ')',
      legend: {position: 'none'}, colors: [COL.stable],
      hAxis: {title: 'Year'}, vAxis: {title: 'Water area (ha)', minValue: 0}
    });
  chart.style().set({height: '260px', stretch: 'horizontal'});
  tsPanel.clear();
  tsPanel.add(chart);
}


// ===========================================================================
// 10. RESET + START
// ===========================================================================
function resetAll() {
  clearDrawn();
  useDrawn = false;
  aoiLabel.setValue('Study area: default demo area (Puzhal & Sholavaram lakes, Chennai region, Tamil Nadu)');
  yearASel.setValue(CFG.defaultYearA);
  yearBSel.setValue(CFG.defaultYearB);
  seasonSel.setValue(DEFAULT_SEASON);
  indexSel.setValue('NDWI');
  thrSlider.setValue(CFG.thrDefault.NDWI);
  opacitySlider.setValue(0.8);
  clearResults();
  tsPanel.clear();
  lastCtx = null;
  overlays = [];
  setBusy(false);
  setStatus('Reset complete. Press "Run Analysis" to start.', 'info');
  showIdleMap();
}

ui.root.clear();
ui.root.add(ui.SplitPanel(panel, map));
showIdleMap();
setStatus('Ready. Press "Run Analysis" to compare the two years.', 'info');
