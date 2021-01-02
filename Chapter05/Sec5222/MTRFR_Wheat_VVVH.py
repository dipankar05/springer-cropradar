# -*- coding: utf-8 -*-

"""
Created on Fri Jun 23 15:04:18 2017

@author: Dipankar
"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.multioutput import MultiOutputRegressor
from sklearn.metrics import mean_absolute_error
import pandas as pd
########################

## Load LUT simulated backscatter intensities (VV and VH) and corresponding crop parameters

Y = pd.read_csv('calibrationcrop.csv',header=None)
X = pd.read_csv('LUT_sigma0.csv',header=None)

# Train-test split for RF      
Xtrain, Xtest, Ytrain, Ytest = train_test_split(X, Y, train_size=0.70)

# RF model
regr_multirf = MultiOutputRegressor(RandomForestRegressor(n_estimators=200, criterion='mse', max_depth=None,bootstrap=True,
                                                         min_impurity_decrease=0.0,
                                                          min_samples_leaf=2,random_state=0)) 

# Fitting RF with training data
regr_multirf.fit(Xtrain, Ytrain)
# Predict on test data
y_multirf = regr_multirf.predict(Xtest)


#rmse estimation
def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())


###-----------------------------------------------------------------------------
###----Testing on test data (splitted previously)--Commented
#bb=Ytest.values
###RMSE
#rmse_val = rmse(np.array(bb[:, 0]), np.array(y_multirf[:, 0]))
###Correlation coefficient 
#corrr=np.corrcoef(np.array(bb[:, 0]), np.array(y_multirf[:, 0]))
#rr= corrr[0,1]
###Plotting data
#plt.plot(bb[:, 0],y_multirf[:, 0], 'ro')
#plt.plot([0, 8], [0, 8], 'k:')
#plt.show()
#
#rmse_val3 = rmse(np.array(bb[:, 1]), np.array(y_multirf[:, 1]))
#corrr3=np.corrcoef(np.array(bb[:, 1]), np.array(y_multirf[:, 1]))
#rr3= corrr3[0,1]
###Plotting data
#plt.plot(bb[:, 1],y_multirf[:, 1], 'ro')
#plt.plot([0, 8], [0, 8], 'k:')
#plt.show()
###-----------------------------------------------------------------------------

###-----------------------------------------------------------------------------
#Final retrieval of PAI and wet biomass
# Read validation data
Yval = pd.read_csv('validation_cropparam.csv',header=None)
Xval = pd.read_csv('validation_sigma0.csv',header=None)

#Estimate with trained MTRFR model
y_out = regr_multirf.predict(Xval)
cc=Yval.values

## PAI estimates-------------------------------------------------------
#rmse estimation
rmse_value = rmse(np.array(cc[:, 0]), np.array(y_out[:, 0]))
#Correlation coefficient 
corrr_value=np.corrcoef(np.array(cc[:, 0]), np.array(y_out[:, 0]))
rr_value= corrr_value[0,1]
mae1=mean_absolute_error(cc[:, 0],y_out[:, 0])

#Plotting
plt.plot(cc[:, 0],y_out[:, 0], 'go')
plt.xlim([0, 8])
plt.ylim([0, 8])
plt.xlabel("Observed PAI ($m^2~m^{-2}$)")
plt.ylabel("Estimated PAI ($m^2~m^{-2}$)")
plt.plot([0, 8], [0, 8], 'k:')
plt.annotate('r = %.2f'%rr_value, xy=(0.5, 7.4))#round off upto 2decimals
plt.annotate('RMSE = %.3f'%rmse_value, xy=(0.5, 6.8))
plt.annotate('MAE = %.3f'%mae1, xy=(0.5, 6.2))
matplotlib.rcParams.update({'font.size': 18})
plt.gca().set_aspect('equal', adjustable='box')
#plt.savefig('PAI_MTRFR.png')
plt.show()

## Wet biomass estimates-------------------------------------------------------
#rmse estimation
rmse_value2 = rmse(np.array(cc[:, 1]), np.array(y_out[:, 1]))
#Correlation coefficient 
corrr_value2=np.corrcoef(np.array(cc[:, 1]), np.array(y_out[:, 1]))
rr_value2= corrr_value2[0,1]
mae2=mean_absolute_error(cc[:, 1],y_out[:, 1])

#Plotting
plt.plot(cc[:, 1],y_out[:, 1], 'go')
plt.xlim([0, 8])
plt.ylim([0, 8])
plt.xlabel("Observed WB ($kg~m^{-2}$)")
plt.ylabel("Estimated WB ($kg~m^{-2}$)")
plt.plot([0, 8], [0, 8], 'k:')
plt.annotate('r = %.2f'%rr_value2, xy=(0.5, 7.4))#round off upto 2decimals
plt.annotate('RMSE = %.3f'%rmse_value2, xy=(0.5, 6.8))
plt.annotate('MAE = %.3f'%mae2, xy=(0.5, 6.2))
matplotlib.rcParams.update({'font.size': 18})
plt.gca().set_aspect('equal', adjustable='box')
#plt.savefig('WB_MTRFR.png')
plt.show()


# End