//Locate bounding cordinates of interested test site
var MB = /* color: #d63000 */ee.Geometry.Polygon(
        [[[-98.14276611583256, 49.81478289632037],
          [-98.14276611583256, 49.34106317137146],
          [-97.76373779552006, 49.34106317137146],
          [-97.76373779552006, 49.81478289632037]]]);
          

// Load Sentinel-1 C-band SAR Ground Range collection 
// var S1collection = ee.ImageCollection('COPERNICUS/S1_GRD').filterBounds(MB)
//     .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))  
//     .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
//     // Filter to get images collected in interferometric wide swath mode
//     .filter(ee.Filter.eq('instrumentMode', 'IW'))
//     // Filter to get images collected in ASCENDING or DECENDING orbit
//     .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
//     .filterDate('2016-06-12', '2016-07-20').sort('system:time');

var s1 = ee.ImageCollection('COPERNICUS/S1_GRD').filterMetadata('instrumentMode', 'equals', 'IW').filterMetadata('orbitProperties_pass', 'equals', 'ASCENDING').
  filter(ee.Filter.eq('transmitterReceiverPolarisation', ['VV', 'VH'])).
  filterBounds(MB).filterDate('2016-06-12', '2016-07-20').
  sort('system:time');
  
//print avaliable list
print('All collection VV+VH IW mode all dates', s1);


//--------------------------------------------------------------------------------------
// Functions to convert from/to dB
function toNatural(img) {
  return ee.Image(10.0).pow(img.select('..').divide(10.0)).copyProperties(img, ['system:time_start'])
}

function toDB(img) {
  return ee.Image(img).log10().multiply(10.0);
}
//--------------------------------------------------------------------------------------
//Apply function on dB image
var S1collection = s1.map(toNatural)
//--------------------------------------------------------------------------------------

// Function to create boxcar 3 x 3 pixel filter
var boxcar = ee.Kernel.square({
 radius: 1.5, units: 'pixels', normalize: true
});
// Function to apply boxcar filter
var fltr = function(image) {
 return image.convolve(boxcar);
};
//--------------------------------------------------------------------------------------
// Apply 3 x 3 pixel mean filter to Pol image
var S1collection = S1collection.map(fltr);


//Create sample points --------------------------------------------
// var p1 = ee.Geometry.Point([-97.9565128,49.6753044]);
// var p2 = ee.Geometry.Point([-97.9640659,49.6704164]);
// var pts = ee.FeatureCollection(ee.List([ee.Feature(p1),ee.Feature(p2)]));


//-----------------------------------------------------------------------------------
//Sampling site Locations--Soybean Fields (SITE_ID)
var Q1 = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point(-97.98110,49.69720), {'label': '41-3'}),
  ee.Feature(ee.Geometry.Point(-97.97830,49.69990), {'label': '41-10'}),
  ee.Feature(ee.Geometry.Point(-97.97830,49.69790), {'label': '41-13'}),
  ee.Feature(ee.Geometry.Point(-97.95620,49.62180), {'label': '65-3'}),
  ee.Feature(ee.Geometry.Point(-97.95200,49.62350), {'label': '65-10'}),
  ee.Feature(ee.Geometry.Point(-97.95520,49.62350), {'label': '65-13'}),
  ee.Feature(ee.Geometry.Point(-97.97250,49.58910), {'label': '71-3'}),
  ee.Feature(ee.Geometry.Point(-97.97530,49.58640), {'label': '71-10'}),
  ee.Feature(ee.Geometry.Point(-97.97530,49.58840), {'label': '71-13'}),
  ee.Feature(ee.Geometry.Point(-97.98090,49.58190), {'label': '72-3'}),
  ee.Feature(ee.Geometry.Point(-97.97810,49.57910), {'label': '72-10'}),
  ee.Feature(ee.Geometry.Point(-97.97810,49.58120), {'label': '72-13'}),
  ee.Feature(ee.Geometry.Point(-97.99540,49.55970), {'label': '82-3'}),
  ee.Feature(ee.Geometry.Point(-97.99820,49.55700), {'label': '82-10'}),
  ee.Feature(ee.Geometry.Point(-97.99830,49.55880), {'label': '82-13'}),
  ee.Feature(ee.Geometry.Point(-97.89950,49.67890), {'label': '101-3'}),
  ee.Feature(ee.Geometry.Point(-97.89850,49.67710), {'label': '101-10'}),
  ee.Feature(ee.Geometry.Point(-97.89540,49.67710), {'label': '101-13'})
  ]);
Map.addLayer(Q1);

var pts = Q1;

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
// Empty Collection to fill
var ft = ee.FeatureCollection(ee.List([]));

//With removal of null values ------------------------------------------
//Function to extract values from image collection based on point file and export as a table 
var fill = function(img, ini) {
// type cast
var inift = ee.FeatureCollection(ini);

// gets the values for the points in the current img
var ft2 = img.reduceRegions(pts, ee.Reducer.first(),30);

// gets the date of the img
var date = img.date().format();

// writes the date in each feature
var ft3 = ft2.map(function(f){return f.set("date", date)});

// merges the FeatureCollections

var ft3a = ft3.filter(ee.Filter.neq('VV', null));//filter first to remove null values
return inift.merge(ft3a);
};
//--------------------------------------------------------------------------------------
// Iterates over the ImageCollection
var newft_remove_null = ee.FeatureCollection(S1collection.iterate(fill, ft));
//print(newft_remove_null);

//--------------------------------------------------------------------------------------
//plot scene and sample points ------------------------------------
Map.setCenter(-97.9565128,49.6753044, 12);
//Map.addLayer(pts);

//--------------------------------------------------------------------------------------
// Export table as .csv (default) in drive
Export.table.toDrive(newft_remove_null,
"BackscatterIntensitiesS1",
"EarthEngine",
"BackscatterIntensities");



// End