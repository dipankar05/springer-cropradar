//GEE4Bio--Carman test site Canada
//Crop biophysical parameter (Plant area index--PAI or wetbiomass--BIOM) estimation from Sentinel-1 data

//--------------------------------------------------------------------------------------
//Reference: Mandal, D., Kumar, V., McNairn, H., Bhattacharya, A. and Rao, Y.S., 2019.
//Joint estimation of Plant Area Index (PAI) and wet biomass in wheat and soybean from C-band polarimetric SAR data. 
//International Journal of Applied Earth Observation and Geoinformation, 79, pp.24-34.
//https://doi.org/10.1016/j.jag.2019.02.007
//--------------------------------------------------------------------------------------

//Test site: Carman (Red river watershed), Manitoba, Canada
//Joint Experiment for Crop Assessment and Monitoring--SAR Intercomparison experiment (Managed by AAFC, CANADA)

//Satellite dataset: Sentinel-1 dual-pol (VV-VH) time series June2016-Aug2016

//Vegetation model: Water Cloud Model
//Inversion approach: Random Forest Regression
//--------------------------------------------------------------------------------------

//Start code
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
var wheatMB = ee.FeatureCollection("users/dipankarIITB/Wheat_MB2016"),
    soybeanMB = ee.FeatureCollection("users/dipankarIITB/Soybean_MB2016"),
    oatsMB = ee.FeatureCollection("users/dipankarIITB/Oats_MB2016"),
    cornMB = ee.FeatureCollection("users/dipankarIITB/Corn_MB2016"),
    canolaMB = ee.FeatureCollection("users/dipankarIITB/Canola_MB2016");

//Locate bounding cordinates of interested test site
var MB = /* color: #d63000 */ee.Geometry.Polygon(
        [[[-98.14276611583256, 49.81478289632037],
          [-98.14276611583256, 49.34106317137146],
          [-97.76373779552006, 49.34106317137146],
          [-97.76373779552006, 49.81478289632037]]]);
          

// Load Sentinel-1 C-band SAR Ground Range collection 
var collection = ee.ImageCollection('COPERNICUS/S1_GRD').filterBounds(MB)
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))  
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  // Filter to get images collected in interferometric wide swath mode.
    .filter(ee.Filter.eq('instrumentMode', 'IW'));

//print avaliable list
//print('All collection VV+VH IW mode all dates', collection);

//Select VV and VH bands and store in variables
var vv=collection.select('VV');
var vh=collection.select('VH');


//Filter by date e.g. this following command selects 2016-07-19 image
var vv0=vv.filterDate('2016-07-18', '2016-07-20').mosaic();
var vh0=vh.filterDate('2016-07-18', '2016-07-20').mosaic();
//--------------------------------------------------------------------------------------
//List of avalible data ranges (You can also find the same at Sci-hub:  https://scihub.copernicus.eu/dhus/ )

//var vv0=vv.filterDate('2016-06-12', '2016-06-14').mosaic();
//var vh0=vh.filterDate('2016-06-12', '2016-06-14').mosaic();

//var vv1=vv.filterDate('2016-06-29', '2016-07-01').mosaic();
//var vh1=vh.filterDate('2016-06-29', '2016-07-01').mosaic();

//var vv2=vv.filterDate('2016-07-02', '2016-07-04').mosaic();
//var vh2=vh.filterDate('2016-07-02', '2016-07-04').mosaic();

//var vv2=vv.filterDate('2016-07-02', '2016-07-04').mosaic();
//var vh2=vh.filterDate('2016-07-02', '2016-07-04').mosaic();

//var vv3=vv.filterDate('2016-07-06', '2016-07-08').mosaic();
//var vh3=vh.filterDate('2016-07-06', '2016-07-08').mosaic();

//var vv4=vv.filterDate('2016-07-14', '2016-07-16').mosaic();
//var vh4=vh.filterDate('2016-07-14', '2016-07-16').mosaic();

//var vv5=vv.filterDate('2016-07-18', '2016-07-20').mosaic();
//var vh5=vh.filterDate('2016-07-18', '2016-07-20').mosaic();

//var vv6=vv.filterDate('2016-07-23', '2016-07-25').mosaic();
//var vh6=vh.filterDate('2016-07-23', '2016-07-25').mosaic();
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

// Display map
Map.centerObject(MB, 12);


//One the display check pixel values, it should be in range of -35 to 0; <0
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//For regression model we need to convert this dB values to linear power scale
//dB to Linear Scale Conversion
// Functions to convert from/to dB
function toNatural(img) {
  return ee.Image(10.0).pow(img.select(0).divide(10.0));
}
///Apply function on dB image
var vv0s = toNatural(vv0).rename('VV');
var vh0s = toNatural(vh0).rename('VH');
  
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

//Speckle filtering/smoothing
// Smooth the image by convolving with the boxcar kernel.
// Define a boxcar or low-pass kernel.
//A 3X3 Boxcar filter
var boxcar = ee.Kernel.square({radius: 1.5, units: 'pixels', normalize: true});
var vh0s1 = vh0s.convolve(boxcar);
var vv0s1 = vv0s.convolve(boxcar);
//similar for other images
//var vh1s = vh1.convolve(boxcar);
//var vv1s = vv1.convolve(boxcar);   
//Display individual data
Map.addLayer(vh0s1.clip(MB), {min:0,max:0.2}, 'VH_despeckled');
//--------------------------------------------------------------------------------------




//--------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------
//Vegetation model
//WCM model is calibrated and LUT is generated from forward WCM in a local system
//For inverse case, load the LUT then build the RFR model with VV-VH and PAI or BIOM

//Load ForwardWCM LUT as table
//This table is going to be used as training
var trainingfeatures = ee.FeatureCollection(wheatMB);
//--------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------
//INPUT: Image
var st_img = vv0s1.addBands(vh0s1);
var inputnames=['VV','VH'];
//--------------------------------------------------------------------------------------
//Let's train the RF model
var RFtrained = ee.Classifier.smileRandomForest({
 numberOfTrees: 74,
 minLeafPopulation: 3,
 bagFraction: 1, 
 seed: 123
})
.setOutputMode('REGRESSION') 
    .train({
    features: trainingfeatures,
    inputProperties: inputnames,
    classProperty: 'PAI'});   //Change classProperty to 'BIOM' if wetbiomass is required
//END of the training part of the code:

//Estimating the biophysical parameter
var estimation=st_img.classify(RFtrained);
//--------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------
//PAI/BIOM colour scale formation
var bioViz = {min: 0, max: 7, palette: ['1700FF','00ECFF','00FF4D', 'F0FF00', 'FF9E00','FF0000']};

//--------------------------------------------------------------------------------------
//Masking with crop inventory map
//Now we will mask the estimates with only Wheat fields

//Load crop mask or you can create one using classification in a seperate code
//Here, for the Carman region we have available crop mask from AAFC in GEE
///AAFC Annual crop inventory
var dataset = ee.ImageCollection('AAFC/ACI').filterBounds(MB);
var crop2016 = dataset
    .filter(ee.Filter.date('2016-01-01', '2016-12-31'))
    .first();
//print('AAFC',crop2016);

//Creating crop mask
var crop_mask = crop2016.eq(146); // create a mask for crop: Wheat==146

//--------------------------------------------------------------------------------------
//Resampling
//Getting projection
var band2 = vv0s1;
//print('CRS:', band2.projection().crs());

// Display a bilinear resampled image with 10m pixel spacing.
var wheatmask_10m = crop_mask.resample('bilinear').reproject({
  crs: band2.projection().crs(),
  scale: 10
});
//--------------------------------------------------------------------------------------


//Estimate wheat biom map
var wheatPAI = estimation.updateMask(wheatmask_10m); // mask it
//print('WheatPAI',wheatPAI);
//Map.addLayer(wheatPAI.clip(MB),bioViz,'WheatLAI Maps');
Map.addLayer(wheatPAI.clip(MB),bioViz,'WheatLAI Maps');

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////Soybean
//-------------------------------------------------------------------------
//https://fusiontables.google.com/data?docid=1kXt2c38y22zHQ4hkniJb__MOhTUtv0mxd_AzVKlo

//Training
//var trainingfeaturessoy=ee.FeatureCollection('ft:1kXt2c38y22zHQ4hkniJb__MOhTUtv0mxd_AzVKlo'); //C11 C22 PAI BIOM H SM (103 samples)

var trainingfeaturessoy = ee.FeatureCollection(soybeanMB);

//Let's train the classifier
var RFtrainedsoy = ee.Classifier.smileRandomForest({
 numberOfTrees: 74,
 minLeafPopulation: 3,
 bagFraction: 1, 
 //outOfBagMode: false,
 seed: 123
})
.setOutputMode('REGRESSION') 
    .train({
    features: trainingfeaturessoy,
    inputProperties: inputnames,
    classProperty: 'PAI'});   
    
//print('RFsoy',RFtrainedsoy);      
////END of the training part of the code:

var estimationsoy=st_img.classify(RFtrainedsoy);
//print('RFestimatesoy',estimationsoy);  

//var bioViz = {min: 0, max: 5, palette: ['00F12F', 'FF0000']};
//Map.addLayer(maps,{min:0,max:10},'LAI Maps');
//Map.addLayer(estimation, bioViz,'LAI Maps',false);


//Map.addLayer(estimation, {min: 0, max: 8, palette: ['blue', 'green', 'red']}, 'custom palette');

///AAFC Annual crop inventory

//Creating crop mask
var soy_mask = crop2016.eq(158); // create a mask for crop: Soybean==158
//print('Soybean',soy_mask);
//Map.addLayer(crop_mask);


//Resampling
//Getting projection
// Display a bilinear resampled image with 10m pixel spacing.
var soymask_10m = soy_mask.resample('bilinear').reproject({
  crs: band2.projection().crs(),
  scale: 10
});
//print('Soybean10m',soymask_10m);
//Map.addLayer(soymask_10m.clip(MB));


//Estimate Soybean biom map
var soyPAI = estimationsoy.updateMask(soymask_10m); // mask it
//print('SoybeanPAI',soyPAI);
Map.addLayer(soyPAI.clip(MB),bioViz,'SoybeanLAI Maps');

///---------------------------------------------------------------------------------------
///---------------------------------------------------------------------------------------
///---------------------------------------------------------------------------------------







//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////Canola
//-------------------------------------------------------------------------
//https://fusiontables.google.com/data?docid=1E1SUYmi09cznHWQoyQ9u088Kempiak1EUwgUPojx

//Training
//var trainingfeaturescan=ee.FeatureCollection('ft:1E1SUYmi09cznHWQoyQ9u088Kempiak1EUwgUPojx'); //C11 C22 PAI BIOM H SM (103 samples)
var trainingfeaturescan = ee.FeatureCollection(canolaMB);

//INPUT: Image
//var st_img = vv0sl.addBands(vh0sl);
//print('Stacked Img', st_img);

//Let's train the classifier
var RFtrainedcan = ee.Classifier.smileRandomForest({
 numberOfTrees: 74,
 minLeafPopulation: 3,
 bagFraction: 1, 
 //outOfBagMode: false,
 seed: 123
})
.setOutputMode('REGRESSION') 
    .train({
    features: trainingfeaturescan,
    inputProperties: inputnames,
    classProperty: 'PAI'});   
    
//print('RFcan',RFtrainedcan);      
////END of the training part of the code:

var estimationcan=st_img.classify(RFtrainedcan);
//print('RFestimatecan',estimationcan);  

//var bioViz = {min: 0, max: 5, palette: ['00F12F', 'FF0000']};
//Map.addLayer(maps,{min:0,max:10},'LAI Maps');
//Map.addLayer(estimation, bioViz,'LAI Maps',false);


//Map.addLayer(estimation, {min: 0, max: 8, palette: ['blue', 'green', 'red']}, 'custom palette');

///AAFC Annual crop inventory

//Creating crop mask
var can_mask = crop2016.eq(153); // create a mask for crop: Canola==153
//print('Canola',can_mask);
//Map.addLayer(crop_mask);


//Resampling
//Getting projection
// Display a bilinear resampled image with 10m pixel spacing.
var canmask_10m = can_mask.resample('bilinear').reproject({
  crs: band2.projection().crs(),
  scale: 10
});
//print('Canola10m',canmask_10m);
//Map.addLayer(canmask_10m.clip(MB));


//Estimate Soybean biom map
var canPAI = estimationcan.updateMask(canmask_10m); // mask it
//print('CanolaPAI',canPAI);
Map.addLayer(canPAI.clip(MB),bioViz,'CanolaLAI Maps');

///---------------------------------------------------------------------------------------
///---------------------------------------------------------------------------------------





//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////Oats
//-------------------------------------------------------------------------
//https://fusiontables.google.com/data?docid=14jOYODzEPR00Klx_1Xru9AhQzqbQ1uQyOT79FSKL

//Training
//var trainingfeaturesoat=ee.FeatureCollection('ft:14jOYODzEPR00Klx_1Xru9AhQzqbQ1uQyOT79FSKL'); //C11 C22 PAI BIOM H SM (103 samples)
var trainingfeaturesoat= ee.FeatureCollection(oatsMB);

//INPUT: Image
//var st_img = vv0sl.addBands(vh0sl);
//print('Stacked Img', st_img);

//Let's train the classifier
var RFtrainedoat = ee.Classifier.smileRandomForest({
 numberOfTrees: 74,
 minLeafPopulation: 3,
 bagFraction: 1, 
 seed: 123
})
.setOutputMode('REGRESSION') 
    .train({
    features: trainingfeaturesoat,
    inputProperties: inputnames,
    classProperty: 'PAI'});   
    
//print('RFoat',RFtrainedoat);      
////END of the training part of the code:

var estimationoat=st_img.classify(RFtrainedoat);
//print('RFestimateoat',estimationoat);  

//var bioViz = {min: 0, max: 5, palette: ['00F12F', 'FF0000']};
//Map.addLayer(maps,{min:0,max:10},'LAI Maps');
//Map.addLayer(estimation, bioViz,'LAI Maps',false);


//Map.addLayer(estimation, {min: 0, max: 8, palette: ['blue', 'green', 'red']}, 'custom palette');

///AAFC Annual crop inventory

//Creating crop mask
var oat_mask = crop2016.eq(136); // create a mask for crop: Oats==136
//print('Oats',oat_mask);
//Map.addLayer(crop_mask);


//Resampling
//Getting projection
// Display a bilinear resampled image with 10m pixel spacing.
var oatmask_10m = oat_mask.resample('bilinear').reproject({
  crs: band2.projection().crs(),
  scale: 10
});
//print('Oats10m',oatmask_10m);
//Map.addLayer(oatmask_10m.clip(MB));


//Estimate Soybean biom map
var oatPAI = estimationoat.updateMask(oatmask_10m); // mask it
//print('OatsPAI',oatPAI);
Map.addLayer(oatPAI.clip(MB),bioViz,'OatsLAI Maps');

///---------------------------------------------------------------------------------------
///---------------------------------------------------------------------------------------





//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////Corn
//-------------------------------------------------------------------------
//https://fusiontables.google.com/data?docid=1Rz4K1_xPFnDVgyyVtB_13I4JW-jSf_p-LrIACKH1

//Training
//var trainingfeaturescorn=ee.FeatureCollection('ft:1Rz4K1_xPFnDVgyyVtB_13I4JW-jSf_p-LrIACKH1'); //C11 C22 PAI BIOM H SM (103 samples)
var trainingfeaturescorn = ee.FeatureCollection(cornMB);

//INPUT: Image
//var st_img = vv0sl.addBands(vh0sl);
//print('Stacked Img', st_img);

//Let's train the classifier
var RFtrainedcorn = ee.Classifier.smileRandomForest({
 numberOfTrees: 74,
 minLeafPopulation: 3,
 bagFraction: 1, 
 seed: 123
})
.setOutputMode('REGRESSION') 
    .train({
    features: trainingfeaturescorn,
    inputProperties: inputnames,
    classProperty: 'PAI'});   
    
//print('RFcorn',RFtrainedcorn);      
////END of the training part of the code:

var estimationcorn=st_img.classify(RFtrainedcorn);
//print('RFestimatecorn',estimationcorn);  

//var bioViz = {min: 0, max: 5, palette: ['00F12F', 'FF0000']};
//Map.addLayer(maps,{min:0,max:10},'LAI Maps');
//Map.addLayer(estimation, bioViz,'LAI Maps',false);


//Map.addLayer(estimation, {min: 0, max: 8, palette: ['blue', 'green', 'red']}, 'custom palette');

///AAFC Annual crop inventory

//Creating crop mask
var corn_mask = crop2016.eq(147); // create a mask for crop: Corn==147
//print('corn',corn_mask);
//Map.addLayer(crop_mask);


//Resampling
//Getting projection
// Display a bilinear resampled image with 10m pixel spacing.
var cornmask_10m = corn_mask.resample('bilinear').reproject({
  crs: band2.projection().crs(),
  scale: 10
});
//print('corn10m',cornmask_10m);
//Map.addLayer(cornmask_10m.clip(MB));


//Estimate Corn biom map
var cornPAI = estimationcorn.updateMask(cornmask_10m); // mask it
//print('CornPAI',cornPAI);
Map.addLayer(cornPAI.clip(MB),bioViz,'CornLAI Maps');

///---------------------------------------------------------------------------------------
///---------------------------------------------------------------------------------------

var mergeredmap=wheatPAI.add(soyPAI).add(canPAI).add(oatPAI).add(cornPAI);
//var mergeredmap=wheatPAI.add(soyPAI.eq(null).eq(1));
//print('mergeredmap',mergeredmap);
var finalmap=mergeredmap.clip(MB);
//print('finalmap',finalmap);
//Map.addLayer(finalmap,bioViz,'LAI Maps');
////

//----------------------------------------------------------------------------------
// Export the image, specifying scale and region.
//Export.image.toDrive({
//  image: finalmap,
  //description: 'BioMaps',
  //scale: 10,
  //region: MB
//});

//Soybean PAI and Export
//var clipsoyPAI=cornPAI.clip(MB);

Export.image.toDrive({
  image: finalmap,
  description: 'BioMaps',
  scale: 20,
  region: MB
});



//--------------------------------------------------------------------------------------
//End




