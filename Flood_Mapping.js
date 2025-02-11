// Sentinel-1 SAR Image Collection Filtering
var filtered_collection = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(study_area)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .filter(ee.Filter.or(
        ee.Filter.eq('orbitProperties_pass', 'DESCENDING'),
        ee.Filter.eq('orbitProperties_pass', 'ASCENDING')
    ));
    
// Function to filter by date
function getFilteredImage(startDate, endDate) {
    return filtered_collection.filterDate(startDate, endDate)
                              .select('VH')
                              .qualityMosaic('VH') // Picks the best image
                              .clip(study_area);
}


// Get Normal and Flood Images
var normal_image = getFilteredImage('2024-07-08', '2024-08-15');
var flood_image = getFilteredImage('2024-08-21', '2024-08-30');


// Convert to natural and apply speckle filter
var normal_image_ft = toDB(RefinedLee(toNatural(normal_image)));
var flood_image_ft = toDB(RefinedLee(toNatural(flood_image)));                        
                            
                            
// Compute Water and Flood Masks efficiently
var threshold = -15;
var normal_threshold = normal_image_ft.gt(threshold);
var flood_threshold = flood_image_ft.lt(threshold);

var floodMask = normal_threshold.and(flood_threshold).updateMask(normal_threshold.and(flood_threshold));
var waterMask = normal_image_ft.lt(threshold).and(flood_image_ft.lt(threshold)).updateMask(normal_image_ft.lt(threshold));

// Map Visualization
Map.centerObject(study_area);
Map.addLayer(study_area, {}, "Study Area");
Map.addLayer(normal_image_ft, {min: -25, max: 0}, "Normal Image");
Map.addLayer(flood_image_ft, {min: -25, max: 0}, "Flood Image");
Map.addLayer(waterMask, {palette: ['Blue']}, "Water Body");
Map.addLayer(floodMask, {palette: ['Red']}, "Flood Water");

// Export the layers
var bounds = study_area.geometry().bounds();
var coords = ee.List(bounds.coordinates().get(0));

// Extract min/max values
var xVals = coords.map(function(coord) { return ee.Number(ee.List(coord).get(0)); }); // Get all X values (longitudes)
var yVals = coords.map(function(coord) { return ee.Number(ee.List(coord).get(1)); }); // Get all Y values (latitudes)

var xmin = ee.Number(xVals.reduce(ee.Reducer.min()));
var xmax = ee.Number(xVals.reduce(ee.Reducer.max()));
var ymin = ee.Number(yVals.reduce(ee.Reducer.min()));
var ymax = ee.Number(yVals.reduce(ee.Reducer.max()));

print("Bounding Box:", xmin, ymin, xmax, ymax);

// Compute midpoint
var xmid = xmin.add(xmax).divide(2);
var ymid = ymin.add(ymax).divide(2);

// Define four tiles
var tile1 = ee.Geometry.Rectangle([xmin, ymid, xmid, ymax]);
var tile2 = ee.Geometry.Rectangle([xmid, ymid, xmax, ymax]);
var tile3 = ee.Geometry.Rectangle([xmin, ymin, xmid, ymid]);
var tile4 = ee.Geometry.Rectangle([xmid, ymin, xmax, ymid]);

// Export each tile separately
var scale = 10; // Adjust scale if needed

Export.image.toDrive({
  image: waterMask.clip(tile1),
  description: "Water_Mask_Tile1",
  scale: scale,
  region: tile1,
  fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image: waterMask.clip(tile2),
  description: "Water_Mask_Tile2",
  scale: scale,
  region: tile2,
  fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image: waterMask.clip(tile3),
  description: "Water_Mask_Tile3",
  scale: scale,
  region: tile3,
  fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image: waterMask.clip(tile4),
  description: "water_Mask_Tile4",
  scale: scale,
  region: tile4,
  fileFormat: 'GeoTIFF'
});

//Utility Function

//declared function to convert the images types and perform the functions 
// Function to convert from d
function toNatural(img) {
  return ee.Image(10.0).pow(img.select(0).divide(10.0));
}

//Function to convert to dB
function toDB(img) {
  return ee.Image(img).log10().multiply(10.0);
}

//Apllying a Refined Lee Speckle filter as coded in the SNAP 3.0 S1TBX:

//https://github.com/senbox-org/s1tbx/blob/master/s1tbx-op-sar-processing/src/main/java/org/esa/s1tbx/sar/gpf/filtering/SpeckleFilters/RefinedLee.java
//Adapted by Guido Lemoine

// by Guido Lemoine
function RefinedLee(img) {
  // img must be in natural units, i.e. not in dB!
  // Set up 3x3 kernels 
  var weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
  var kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, false);

  var mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3);
  var variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3);

  // Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
  var sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]]);

  var sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, false);

  // Calculate mean and variance for the sampled windows and store as 9 bands
  var sample_mean = mean3.neighborhoodToBands(sample_kernel); 
  var sample_var = variance3.neighborhoodToBands(sample_kernel);

  // Determine the 4 gradients for the sampled windows
  var gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs();
  gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs());
  gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs());
  gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs());

  // And find the maximum gradient amongst gradient bands
  var max_gradient = gradients.reduce(ee.Reducer.max());

  // Create a mask for band pixels that are the maximum gradient
  var gradmask = gradients.eq(max_gradient);

  // duplicate gradmask bands: each gradient represents 2 directions
  gradmask = gradmask.addBands(gradmask);

  // Determine the 8 directions
  var directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1);
  directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2));
  directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3));
  directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4));
  // The next 4 are the not() of the previous 4
  directions = directions.addBands(directions.select(0).not().multiply(5));
  directions = directions.addBands(directions.select(1).not().multiply(6));
  directions = directions.addBands(directions.select(2).not().multiply(7));
  directions = directions.addBands(directions.select(3).not().multiply(8));

  // Mask all values that are not 1-8
  directions = directions.updateMask(gradmask);

  // "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
  directions = directions.reduce(ee.Reducer.sum());  

  //var pal = ['ffffff','ff0000','ffff00', '00ff00', '00ffff', '0000ff', 'ff00ff', '000000'];
  //Map.addLayer(directions.reduce(ee.Reducer.sum()), {min:1, max:8, palette: pal}, 'Directions', false);

  var sample_stats = sample_var.divide(sample_mean.multiply(sample_mean));

  // Calculate localNoiseVariance
  var sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]);

  // Set up the 7*7 kernels for directional statistics
  var rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4));

  var diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0], 
    [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]]);

  var rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, false);
  var diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, false);

  // Create stacks for mean and variance using the original kernels. Mask with relevant direction.
  var dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1));
  var dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1));

  dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)));
  dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)));

  // and add the bands for rotated kernels
  for (var i=1; i<4; i++) {
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
  }

  // "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
  dir_mean = dir_mean.reduce(ee.Reducer.sum());
  dir_var = dir_var.reduce(ee.Reducer.sum());

  // A finally generate the filtered value
  var varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0));

  var b = varX.divide(dir_var);

  var result = dir_mean.add(b.multiply(img.subtract(dir_mean)));
  return(result.arrayFlatten([['sum']]));
}