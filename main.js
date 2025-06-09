// ----- Code to obtain images, perform linear regression, calculate residuals, downscale LST -----

// Define the center point and creating ROI
var lon = -123.08;
var lat = 49.24;
var point = ee.Geometry.Point([lon, lat]);
var roi = point.buffer(15000).bounds();

// Defining start and end dates
var startDate = '2020-07-28';
var endDate = '2020-07-31';

// Obtain correct images for L8 and S2
var L8_collection = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
                  .filterBounds(roi)
                  .filterDate(startDate, endDate)
                  .filterMetadata("CLOUD_COVER", "less_than", 10);
var l8Projection = L8_collection.first().projection();
var L8Image = L8_collection.first();

var tile_ids = ['10UDV', '10UEV'];
var S2_collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                  .filterDate('2020-07-29', '2020-07-30')
                  .filter(ee.Filter.inList('MGRS_TILE', tile_ids));
var s2Projection = S2_collection.first().select('B4').projection();
var S2Image = S2_collection.mosaic();

// Landsat image scaling
function landsatScaling(image) {
  var opticalBands = image.select("SR_B.").multiply(0.0000275).add(-0.2).multiply(10000);
  var thermalBands = image.select("ST_B.*").multiply(0.00341802).add(149.0).subtract(273.15);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
}
var L8ImageScaled = landsatScaling(L8Image);
var L8ImageScaledClipped = L8ImageScaled.clip(roi);

// Sentinel image scaling
var S2ImageScaled = S2Image.clip(roi).divide(10000);

// Landsat indices and LST
var L8_NDVI = L8ImageScaledClipped.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');
var L8_NDWI = L8ImageScaledClipped.normalizedDifference(['SR_B3', 'SR_B5']).rename('NDWI');
var L8_NDBI = L8ImageScaledClipped.normalizedDifference(['SR_B6', 'SR_B5']).rename('NDBI');
var L8_LST = L8ImageScaledClipped.select('ST_B10').rename('L8_LST');

// Sentinel indices
var S2_NDVI = S2ImageScaled.normalizedDifference(['B8', 'B4']).rename('NDVI');
var S2_NDWI = S2ImageScaled.normalizedDifference(['B3', 'B11']).rename('NDWI');
var S2_NDBI = S2ImageScaled.normalizedDifference(['B11', 'B8']).rename('NDBI');

// Projections
var L8_projection = L8_NDVI.projection();
var S2_projection = S2_NDVI.projection();

// Creating image chart and random points, with seed
var L8ImageChart = L8_LST.addBands([L8_NDVI, L8_NDWI, L8_NDBI]);
var randomPoints = ee.FeatureCollection.randomPoints(roi, 5000, 40);

var pointValues = L8ImageChart.reduceRegions({
      collection: randomPoints,
      reducer: ee.Reducer.mean(),
      scale: 30,
  });
  
var getValues = ee.FeatureCollection(pointValues.filter(ee.Filter.notNull(['L8_LST', 'NDVI', 'NDWI', 'NDBI'])));

var L8_LST_values = getValues.aggregate_array('L8_LST');
var L8_NDVI_values = getValues.aggregate_array('NDVI');
var L8_NDWI_values = getValues.aggregate_array('NDWI');
var L8_NDBI_values = getValues.aggregate_array('NDBI');

// Linear Regression
var bands = ee.Image(1).addBands(L8_NDVI).addBands(L8_NDWI).addBands(L8_NDBI).addBands(L8_LST).rename(['constant', 'NDVI', 'NDWI', 'NDBI', 'L8']);

var imageRegression = bands.reduceRegion({
  reducer: ee.Reducer.linearRegression({numX:4, numY:1}),
  geometry: roi,
  scale: 30,
});

var coefficientsList = ee.Array(imageRegression.get('coefficients')).toList();
var intercept = ee.Image(ee.Number(ee.List(coefficientsList.get(0)).get(0)));
var interceptList = ee.List(coefficientsList.get(0)).get(0);
var slopeNDVI = ee.Image(ee.Number(ee.List(coefficientsList.get(1)).get(0)));
var slopeNDVIList = ee.List(coefficientsList.get(1)).get(0);
var slopeNDWI = ee.Image(ee.Number(ee.List(coefficientsList.get(2)).get(0)));
var slopeNDWIList = ee.List(coefficientsList.get(2)).get(0);
var slopeNDBI = ee.Image(ee.Number(ee.List(coefficientsList.get(3)).get(0)));
var slopeNDBIList = ee.List(coefficientsList.get(3)).get(0);

print(interceptList, "intercept", "",
      slopeNDVIList, "slope NDVI", "", 
      slopeNDWIList, "slope NDWI", "",
      slopeNDBIList, "slope NDWBI");

var modelLST_30m_unprojected = intercept.add(slopeNDVI.multiply(L8_NDVI))
                                     .add(slopeNDWI.multiply(L8_NDWI))
                                     .add(slopeNDBI.multiply(L8_NDBI))
                                     .clip(roi);
var modelLST_30m = modelLST_30m_unprojected.reproject({crs: l8Projection});

var downscaledLST_10m = intercept.add(slopeNDVI.multiply(S2_NDVI))
                                     .add(slopeNDWI.multiply(S2_NDWI))
                                     .add(slopeNDBI.multiply(S2_NDBI))
                                     .clip(roi)
                                     .reproject({crs: s2Projection});

// Residuals
var modelResiduals = L8_LST.subtract(modelLST_30m);
var modelResiduals_30m = modelResiduals.reproject({crs: l8Projection});

var gaussianKernel = ee.Kernel.gaussian({radius: 1.5, units: 'pixels'});
var downscaledResiduals = modelResiduals_30m.resample('bicubic').reproject({crs: s2Projection}).convolve(gaussianKernel);

// Final downscaled LST
var downscaledLST_w_residuals = downscaledLST_10m.add(downscaledResiduals);
