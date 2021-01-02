# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 09:23:28 2018

@author: Administrator
"""
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error
import pandas as pd
########################

###Load and create LUT database first
Y1 = pd.read_csv('calibrationcrop.csv',header=None)
#Pandas dataframe to matrix conversion
Ym=Y1.values
#Incidence angle
thr=Ym[:,3]#Col3==Local Inc Angle in rad
thr = np.float64(thr)
Ym0=np.float64(Ym[:,0])#Col0==PAI_True
#Ym1=np.float64(Ym[:,1])#Col1==Wetbiomass;
Ym2=np.float64(Ym[:,2])#Col2==Soilmoisture

##-----------------------------------------------------------------------------
Y=np.column_stack((Ym0,Ym2))
#else if one parameter
#Y=Ym0


#XHH = pd.read_csv('HHsimulatedCanola.csv',header=None)
#XHH=np.float64(XHH.as_matrix(columns=None))
#XHV = pd.read_csv('HVsimulatedCanola.csv',header=None)
#XHV=np.float64(XHV.as_matrix(columns=None))
XVH = pd.read_csv('VHsimulatedCorn.csv',header=None)
XVH=np.float64(XVH.values)
XVV = pd.read_csv('VVsimulatedCorn.csv',header=None)
XVV=np.float64(XVV.values)
#X=np.column_stack((XHH,XHV,XVH,XVV))
X=np.column_stack((XVH,XVV))


#  MTRFR Model formation    
Xtrain, Xtest, Ytrain, Ytest = train_test_split(X, Y, train_size=0.9)

regr_multirf = RandomForestRegressor(n_estimators=150, criterion='mse', max_depth=None,bootstrap=True,
                                                         min_impurity_decrease=0.0,
                                                          min_samples_leaf=4,random_state=3)

#https://stackoverflow.com/questions/28064634/random-state-pseudo-random-numberin-scikit-learn
regr_multirf.fit(Xtrain, Ytrain)
# Predict on new data
y_multirf = regr_multirf.predict(Xtest)



########################################################################################################
###Load validation data
###Load validation data
#replacing 'no info' and '.' i.e. blank space or 'None' string with NaN
cornval = pd.read_excel('ValidationPoints.xlsx',na_values = ['no info', '.','None','#VALUE!', '#DIV/0!'],skiprows=[0],header=None);
vald=cornval.dropna(subset=[1])                                                              
                                                              
valdm=vald.values

#Incidence angle
valthr=(3.1415/180)*valdm[:,0]#Col15==Local Inc Angle
valthr = np.float64(valthr)
#valHH=np.float64(valdm[:,17])#Col16==HH; 17==HV; 18==VH; 19==VV
#valHV=np.float64(valdm[:,18])
valVH=np.float64(valdm[:,4])
valVV=np.float64(valdm[:,3])
vallai=np.float64(valdm[:,1])#Col6==PAI_True; Col7==LAI
valbiom=np.float64(valdm[:,1])#Col9==Wetbiomass; Col8==VWC; Col10==drybiomass
valsm=np.float64(valdm[:,2])#Col5==Soilmoisture
valY=np.column_stack((vallai,valsm))
#valY=vallai

#valX=np.column_stack((valHH,valHV,valVH,valVV))
valX=np.column_stack((valVH,valVV))

# Predictfor validation data
y_out = regr_multirf.predict(valX)

######################################3
##PAI estimation and error
#rmse estimation
def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())
rmselai = rmse(np.array(vallai), np.array(y_out[:, 0]))
#Correlation coefficient 
corrr_value=np.corrcoef(np.array(vallai), np.array(y_out[:, 0]))
rrlai= corrr_value[0,1]
maelai=mean_absolute_error(vallai,y_out[:, 0])
#Plotting
plt.plot(vallai,y_out[:, 0], 'go')
plt.xlim([0, 6])
plt.ylim([0, 6])
plt.xlabel("Observed LAI ($m^2 m^{-2}$)")
plt.ylabel("Estimated LAI ($m^2 m^{-2}$)")
plt.plot([0, 6], [0, 6], 'k:')
plt.annotate('r = %.2f'%rrlai, xy=(0.5, 5.5))#round off upto 3decimals
plt.annotate('RMSE = %.2f'%rmselai, xy=(0.5, 5.0))
plt.annotate('MAE = %.2f'%maelai, xy=(0.5, 4.5))
matplotlib.rcParams.update({'font.size': 20})
plt.yticks(np.arange(0, 7, 2))
plt.xticks(np.arange(0, 7, 2))
plt.gca().set_aspect('equal', adjustable='box')
#plt.savefig('PAIValidationLUT.png',bbox_inches="tight",dpi=100)
plt.show()
