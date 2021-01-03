# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 09:23:28 2018

@author: Administrator
"""
#import matplotlib.pyplot as plt
#import matplotlib
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error
import pandas as pd
########################
from sklearn.model_selection import GridSearchCV



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
Xtrain, Xtest, Ytrain, Ytest = train_test_split(X, Y, train_size=0.8)

regr_multirf = RandomForestRegressor(n_estimators=150, criterion='mse', max_depth=None,bootstrap=True,
                                                         min_impurity_decrease=0.0,
                                                          min_samples_leaf=4,random_state=3)

#https://stackoverflow.com/questions/28064634/random-state-pseudo-random-numberin-scikit-learn
regr_multirf.fit(Xtrain, Ytrain)
# Predict on new data
y_multirf = regr_multirf.predict(Xtest)


##--------------------------------------------------
# Find the best parameters for the RFR model and number of Tree selection
parameters = {
    #'max_depth': [70, 80, 90, 100],
    'n_estimators': [50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350]
}
gridforest = GridSearchCV(regr_multirf, parameters, cv = 3, n_jobs = -1, verbose = 1)
gridforest.fit(Xtrain, Ytrain)
print('Best no. of trees:',gridforest.best_params_)
##-----------------------------------------------------------------------------

MAEtree=mean_absolute_error(Ytest,y_multirf)
print('MAE=',MAEtree)