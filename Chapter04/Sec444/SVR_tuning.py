# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 09:23:28 2018

@author: Administrator
"""
#import matplotlib.pyplot as plt
#import matplotlib
import numpy as np
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler 
from sklearn.metrics import mean_absolute_error
import pandas as pd
from sklearn.pipeline import make_pipeline

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
#Y=np.column_stack((Ym0,Ym2))
Y=Ym0


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


###################MSVR########################################################################
###############MSVR
#svr = SVR(kernel='rbf', C=1000, gamma=0.05)
pipeline = make_pipeline(StandardScaler(),
    SVR(kernel='rbf', epsilon=0.5, C=4, gamma = 0.05),
)
MSVRmodel=pipeline.fit(X,Y)
y_out = pipeline.predict(X)


#svr = SVR(kernel='rbf', C=1000, gamma=1.1, epsilon=0.09)

#svr_rbf = SVR(kernel='rbf', C=1e3, gamma=0.1)
#svr_lin = SVR(kernel='linear', C=1e3)
#svr_poly = SVR(kernel='poly', C=1e3, degree=2)
#Intuitively, the gamma parameter defines how far the influence of a single training example reaches, with low values meaning ‘far’ and high values meaning ‘close’. The gamma parameters can be seen as the inverse of the radius of influence of samples selected by the model as support vectors.
#The C parameter trades off correct classification of training examples against maximization of the decision function’s margin. For larger values of C, a smaller margin will be accepted if the decision function is better at classifying all training points correctly. A lower C will encourage a larger margin, therefore a simpler decision function, at the cost of training accuracy. In other words``C`` behaves as a regularization parameter in the SVM.
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#regr = svr
#MSVRmodel=regr.fit(X,Y);


#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#Ypredtrain=regr.predict(X);

##--------------------------------------------------
# Find the best parameters for the RFR model and number of Tree selection
parameters = {
        'kernel': ['rbf','linear'],
        'epsilon' : [0.1, 0.3, 0.5],
        'gamma': [0.01, 0.05, 0.1],
        'C': [1, 5, 10]
}
gridSVR = GridSearchCV(SVR(), parameters, cv = 5, n_jobs = -1, verbose = 1)
gridSVR.fit(X,Y)
print('Best parameters:',gridSVR.best_params_)
##-----------------------------------------------------------------------------
MAEtree=mean_absolute_error(Y,y_out)
print('MAE=',MAEtree)