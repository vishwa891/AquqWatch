# AquaWatch – Satellite-Based Surface Water Change Detection

Monitoring water-body expansion and contraction using Sentinel-2 satellite imagery.

**Hackathon problem:** Domain 2 – Water Resource Management, Problem 2.1 – Surface Water Change Detection
**Platform:** Google Earth Engine (JavaScript Code Editor + Earth Engine App)

## What it does
AquaWatch lets a user pick an Area of Interest (AOI) and compare surface-water extent between two years. It extracts water with NDWI/MNDWI, classifies change, calculates areas in hectares, highlights loss/gain hotspots and gives a rule-based recommendation.

## Features
- Default demo area: Puzhal & Sholavaram lakes, Chennai region (Tamil Nadu), or draw your own rectangle/polygon
- Historical and recent year selectors, season window, NDWI/MNDWI switch, threshold slider, opacity slider
- Change classes: 0 no significant water, 1 water loss, 2 water gain, 3 stable water
- Statistics cards: historical area, recent area, gain, loss, net change, % change
- Bar chart (area summary) and optional yearly water-area trend
- Loss/gain hotspot grid cells with a "top 5" list and zoom buttons
- Automatic decision-support text (data-driven, no claims about causes)
- Error handling: empty collections, invalid years, tiny/huge/invalid AOI, no water, low cloud-free coverage

## Workflow
Sentinel-2 SR → cloud masking + median composite → NDWI/MNDWI → water mask → change classification → area statistics + hotspot grid → interactive dashboard → recommendation

## Data
`COPERNICUS/S2_SR_HARMONIZED` (Sentinel-2 Level-2A). Bands: B2, B3, B4, B8, B11, SCL. Scene cloud filter < 60 %; SCL cloud/cirrus/saturated pixels masked.

## Method
- NDWI = (Green − NIR) / (Green + NIR) = (B3 − B8) / (B3 + B8)
- MNDWI = (Green − SWIR1) / (Green + SWIR1) = (B3 − B11) / (B3 + B11)
- Water if index > threshold (default NDWI 0.05, MNDWI 0.00, adjustable)
- `change = waterHistorical + 2 × waterRecent` (only where both years have valid data)
- Area = `ee.Image.pixelArea()` ÷ 10 000 (hectares), summed with `reduceRegion()`
- Net change = recent − historical; % change = net ÷ historical × 100
- Hotspot = grid cell where lost (or gained) water ≥ 15 % of the cell area

## How to run
1. Open https://code.earthengine.google.com/
2. Create a new script, paste `aquawatch.js`, click **Run**
3. Press **Run Analysis** in the left panel
4. To publish: **Apps → NEW APP** (see the guide)

Live app: https://natural-nebula-510204-f6.projects.earthengine.app/view/aquawatch

## Limitations
- Detects water-surface change only; it does **not** prove causes (drought, encroachment, construction, dam operation)
- Rainfall and season differ between years; same-season windows reduce but do not remove this
- NDWI/MNDWI can confuse shadows, wet soil, algae and turbid water
- 10 m pixels miss very small ponds; areas above ~300 km² are analysed at 20–30 m
- No ground-truth validation

## Future enhancements
JRC Global Surface Water validation, Sentinel-1 SAR for cloudy monsoon scenes, Otsu automatic thresholds, per-lake statistics from a waterbody layer, rainfall (CHIRPS) overlay, PDF/CSV export, monthly time series.

## Team
_Add your name / team name here_

## License
MIT

Live Project App Link:https://natural-nebula-510204-f6.projects.earthengine.app/view/aquawatch


