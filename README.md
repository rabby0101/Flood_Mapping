# Sentinel-1 SAR Flood Mapping using Google Earth Engine (GEE)

This repository contains a Google Earth Engine (GEE) script that performs flood mapping using Sentinel-1 SAR imagery. The script filters Sentinel-1 data, applies preprocessing steps including speckle filtering, and generates flood and water masks based on a predefined threshold.

## Features
- **Sentinel-1 SAR Image Filtering**: Filters Sentinel-1 SAR GRD images based on area of interest, polarization, and orbit pass.
- **Preprocessing**: Converts SAR images to natural units, applies a refined Lee speckle filter, and converts the images back to dB scale.
- **Flood and Water Mask Generation**: Uses a predefined threshold to classify flooded and water areas.
- **Visualization**: Displays processed images and masks on the GEE map.
- **Data Export**: Exports classified water and flood masks as GeoTIFF files, divided into four tiles for better resolution.

## Data Processing Workflow
### 1. Sentinel-1 SAR Image Collection Filtering
```javascript
var filtered_collection = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(study_area)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .filter(ee.Filter.or(
        ee.Filter.eq('orbitProperties_pass', 'DESCENDING'),
        ee.Filter.eq('orbitProperties_pass', 'ASCENDING')
    ));
```

### 2. Image Selection by Date
```javascript
function getFilteredImage(startDate, endDate) {
    return filtered_collection.filterDate(startDate, endDate)
                              .select('VH')
                              .qualityMosaic('VH')
                              .clip(study_area);
}
```

- `normal_image`: Represents normal conditions.
- `flood_image`: Represents the flooded conditions.

### 3. Speckle Filtering and Conversion
```javascript
var normal_image_ft = toDB(RefinedLee(toNatural(normal_image)));
var flood_image_ft = toDB(RefinedLee(toNatural(flood_image)));
```

Functions used:
- `toNatural(img)`: Converts dB to natural scale.
- `toDB(img)`: Converts natural scale to dB.
- `RefinedLee(img)`: Applies the refined Lee speckle filter.

### 4. Flood and Water Mask Generation
```javascript
var threshold = -15;
var normal_threshold = normal_image_ft.gt(threshold);
var flood_threshold = flood_image_ft.lt(threshold);

var floodMask = normal_threshold.and(flood_threshold).updateMask(normal_threshold.and(flood_threshold));
var waterMask = normal_image_ft.lt(threshold).and(flood_image_ft.lt(threshold)).updateMask(normal_image_ft.lt(threshold));
```

### 5. Map Visualization
```javascript
Map.centerObject(study_area);
Map.addLayer(study_area, {}, "Study Area");
Map.addLayer(normal_image_ft, {min: -25, max: 0}, "Normal Image");
Map.addLayer(flood_image_ft, {min: -25, max: 0}, "Flood Image");
Map.addLayer(waterMask, {palette: ['Blue']}, "Water Body");
Map.addLayer(floodMask, {palette: ['Red']}, "Flood Water");
```

### 6. Exporting Results
The flood and water masks are exported as GeoTIFF files, divided into four tiles.

```javascript
Export.image.toDrive({
  image: waterMask.clip(tile1),
  description: "Water_Mask_Tile1",
  scale: 10,
  region: tile1,
  fileFormat: 'GeoTIFF'
});
```

This process is repeated for `tile2`, `tile3`, and `tile4` to export all tiles separately.

## Utility Functions
The repository includes utility functions for SAR image processing:
- **`toNatural(img)`**: Converts SAR image from dB to natural scale.
- **`toDB(img)`**: Converts SAR image from natural scale to dB.
- **`RefinedLee(img)`**: Implements the refined Lee speckle filter.

## Requirements
To run this script, you need:
- A Google Earth Engine (GEE) account.
- A defined `study_area` in your GEE script.

## Usage
1. Open [Google Earth Engine](https://code.earthengine.google.com/).
2. Copy and paste the script into the code editor.
3. Define your `study_area`.
4. Run the script to visualize the results.
5. Download the exported GeoTIFFs from Google Drive.

## References
- [Google Earth Engine Documentation](https://developers.google.com/earth-engine/)
- [Sentinel-1 Data](https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S1_GRD)
- [Speckle Filtering in SAR](https://github.com/senbox-org/s1tbx)

## License
This project is licensed under the MIT License.

