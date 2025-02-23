var UScounty_shp = ee.FeatureCollection("projects/operating-edge-294605/assets/UScounty");
var AOD_datasets = ee.ImageCollection("MODIS/006/MCD19A2_GRANULES");
var geometry = ee.Geometry.Polygon(
    [[[-126.78896181946607, 49.707049503176094],
      [-126.78896181946607, 23.458674275459153],
      [-66.05653994446607, 23.458674275459153],
      [-66.05653994446607, 49.707049503176094]]], null, false);

// Define the date range
for (var year=2001; year<=2022; year++) {
    AOD(year);
  }
  
  function AOD(year) {
    
    var startDate = ee.Date(year+'-01-01');
    var endDate = ee.Date(year+1+'-01-01');
    var processedAerosol = AOD_datasets.filterDate(startDate, endDate).select('Optical_Depth_055');
    var numberOfDays = endDate.difference(startDate, 'days');
    
    var CalculationDaily = function (dayOffset) {
      var start = startDate.advance(dayOffset, 'days');
      var end = start.advance(1, 'days');
      return processedAerosol
        .filterDate(start, end)
        .mean().rename(start.format('YYYY-MM-dd'))};
    
    var daily = ee.ImageCollection(
      ee.List.sequence(0, numberOfDays.subtract(1))
      .map(CalculationDaily)
    );
    
    
    // Define a function to extract the GEOID from the feature properties
    var extractGEOID = function(feature) {
      var geoid = feature.get('GEOID');
      return feature.set('GEOID', geoid);
    };
    
    // Map the function over the county feature collection
    var countyWithGEOID = UScounty_shp.map(extractGEOID);
    
    var pm25Mean = daily.toBands().reduceRegions({
      collection: countyWithGEOID,
      reducer: ee.Reducer.mean(),
      scale: 10000
    });
    
    // Export the results as a CSV file
    Export.table.toDrive({
      collection: pm25Mean,
      description: 'UScounty_AOD_' + year,
      folder: "US heatwave",
      fileFormat: 'CSV'
    });
  }
  
// add visualisation functions
var band_viz = {
    min: 0,
    max: 500,
    palette: ['black', 'blue', 'purple', 'cyan', 'green', 'yellow', 'red']
};
var AOD_image_2023 = AOD_datasets.select('Optical_Depth_055').filterDate('2023-01-01', '2023-12-31').mean().clip(geometry);

Map.centerObject(geometry, 4);
Map.addLayer(AOD_image_2023, band_viz, 'Aerosol optical depth in 2023');
Map.addLayer(UScounty_shp, {color: 'black', fillColor: '00000000', width:0.3}, 'AOI');
