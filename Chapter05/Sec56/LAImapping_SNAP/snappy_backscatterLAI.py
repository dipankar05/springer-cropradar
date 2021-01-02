import sys

import numpy
import numpy as np
from snappy import Product
from snappy import ProductData
from snappy import ProductIO
from snappy import ProductUtils
from snappy import FlagCoding
##############
import csv
###############MSVR
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler 
from sklearn.pipeline import make_pipeline
########################

if len(sys.argv) != 2:
    print("usage: %s <file>" % sys.argv[0])
    sys.exit(1)

file = sys.argv[1]

print("Reading...")
product = ProductIO.readProduct(file)
width = product.getSceneRasterWidth()
height = product.getSceneRasterHeight()
name = product.getName()
description = product.getDescription()
band_names = product.getBandNames()

print("Product:     %s, %s" % (name, description))
print("Raster size: %d x %d pixels" % (width, height))
print("Start time:  " + str(product.getStartTime()))
print("End time:    " + str(product.getEndTime()))
print("Bands:       %s" % (list(band_names)))
##---------------------------------------------------------------------------------
with open('rice_LUT.csv','r') as dest_f:
    data_iter = csv.reader(dest_f,
                           delimiter = ',',
                           quotechar = '"')
    data = [data for data in data_iter]
data_array = np.asarray(data, dtype = np.float32)

VV = data_array[:,1]
VH = data_array[:,2]
PAI = data_array[:,0]

X=np.column_stack((VV,VH))
Y = PAI

#SVR training
pipeline = make_pipeline(StandardScaler(),
    SVR(kernel='rbf', epsilon=0.105, C=250, gamma = 2.8),
)
SVRmodel=pipeline.fit(X,Y)

# Predictfor validation data
valX = X;
y_out = pipeline.predict(valX);



##---------------------------------------------------------------------------------
bandc11 = product.getBand('C11')
bandc22 = product.getBand('C22')

laiProduct = Product('LAI', 'LAI', width, height)
laiBand = laiProduct.addBand('lai', ProductData.TYPE_FLOAT32)
laiFlagsBand = laiProduct.addBand('lai_flags', ProductData.TYPE_UINT8)
writer = ProductIO.getProductWriter('BEAM-DIMAP')

ProductUtils.copyGeoCoding(product, laiProduct)
ProductUtils.copyMetadata(product, laiProduct)
ProductUtils.copyTiePointGrids(product, laiProduct)

laiFlagCoding = FlagCoding('lai_flags')
laiFlagCoding.addFlag("LAI_LOW", 1, "LAI below 0")
laiFlagCoding.addFlag("LAI_HIGH", 2, "LAI above 5")
group = laiProduct.getFlagCodingGroup()
#print(dir(group))
group.add(laiFlagCoding)

laiFlagsBand.setSampleCoding(laiFlagCoding)

laiProduct.setProductWriter(writer)
laiProduct.writeHeader('LAImap_output.dim')

c11 = numpy.zeros(width, dtype=numpy.float32)
c22 = numpy.zeros(width, dtype=numpy.float32)

print("Writing...")

for y in range(height):
    print("processing line ", y, " of ", height)
    c11 = bandc11.readPixels(0, y, width, 1, c11)
    c22 = bandc22.readPixels(0, y, width, 1, c22)
    
    Z=np.column_stack((c11,c22))

    #ndvi = (r10 - r7) / (r10 + r7)
    
    lai = pipeline.predict(Z);
    laiBand.writePixels(0, y, width, 1, lai)
    laiLow = lai < 0.0
    laiHigh = lai > 5.0
    laiFlags = numpy.array(laiLow + 2 * laiHigh, dtype=numpy.int32)
    laiFlagsBand.writePixels(0, y, width, 1, laiFlags)

laiProduct.closeIO()

print("Done.")