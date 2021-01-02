# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 19:36:25 2020

@author: Dipankar
"""

###Importing library
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit
import pandas as pd
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler 
from sklearn.metrics import mean_absolute_error
from sklearn.pipeline import make_pipeline


#Load CAL data
#With headers
riceallhead = pd.read_excel('DSR_Rice_SAR_CAL_combine.xlsx',na_values = ['no info', '.','None','#VALUE!', '#DIV/0!'],skiprows=[0],header=None);

#Pandas dataframe to matrix conversion
cald1m=riceallhead.values


#Incidence angle
thr=31*np.pi/180#Col15==Local Inc Angle
##################################################################################################
######
##PAI values
x1=np.float64(cald1m[:,11])#Col8==PAI
RH0=np.float64(cald1m[:,13])#Col16==HH; 17==HV; 18==VH; 19==VV
RV0=np.float64(cald1m[:,14])
##################################################################################################
##M-chi powers
Psm0=np.float64(cald1m[:,20])
Pdm0=np.float64(cald1m[:,21])
Pvm0=np.float64(cald1m[:,22])
##################################################################################################
##IS-Omega powers
Psi0=np.float64(cald1m[:,17])#Col6==Ps
Pdi0=np.float64(cald1m[:,16])#Col9==Pd
Pvi0=np.float64(cald1m[:,18])#Col5==Pv

##################################################################################################
##################################################################################################

### Forward modeling
#
###################################################################################################



##################################################################################################
##################################################################################################

## SVR model building and Inversion

##################################################################################################
##################################################################################################

# Read VAL data
riceallhead2 = pd.read_excel('DSR_Rice_SAR_VAL_combine.xlsx',na_values = ['no info', '.','None','#VALUE!', '#DIV/0!'],skiprows=[0],header=None);
cald1mv=riceallhead2.values
## RH-RV VAL data--------------------------------------------------------------
RHv=np.float64(cald1mv[:,13])#Col16==HH; 17==HV; 18==VH; 19==VV
RVv=np.float64(cald1mv[:,14])
## m-Chi VAL data--------------------------------------------------------------
Psmv=np.float64(cald1mv[:,20])
Pdmv=np.float64(cald1mv[:,21])
Pvmv=np.float64(cald1mv[:,22])
## iS_Omega VAL data--------------------------------------------------------------
Psiv=np.float64(cald1mv[:,17])#Col6==Ps
Pdiv=np.float64(cald1mv[:,16])#Col9==Pd
Pviv=np.float64(cald1mv[:,18])#Col5==Pv
## PAI VAL data--------------------------------------------------------------
x1v=np.float64(cald1mv[:,11])#Col8==PAI

##################################################################################################
## Building SVR models
## RH-RV SVR model-------------------------------------------------------------
Xr=np.column_stack((RH0,RV0))
valXr=np.column_stack((RHv,RVv))

pipelineRHRV = make_pipeline(StandardScaler(),SVR(kernel='rbf', epsilon=0.02, C=1000, gamma = 0.02),)
#pipelineRHRV = make_pipeline(StandardScaler(),SVR(kernel='rbf', epsilon=0.2, C=1000, gamma = 0.2),)
SVRmodel=pipelineRHRV.fit(Xr,x1)
# Predict for validation data
y_out = pipelineRHRV.predict(valXr);

## -------------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------------
# Read VAL data
riceallhead3 = pd.read_excel('PAI_estmc.xlsx',na_values = ['no info', '.','None','#VALUE!', '#DIV/0!'],skiprows=[0],header=None);
cald1mv3=riceallhead3.values
## RH-RV VAL data--------------------------------------------------------------
x1v=np.float64(cald1mv3[:,0])
y_out=np.float64(cald1mv3[:,1])
y_out1=np.float64(cald1mv3[:,2])
y_out2=np.float64(cald1mv3[:,3])


##PAI estimation and error-------------------------------------------------------------
def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())
rmselai = rmse(np.array(x1v), np.array(y_out))
#Correlation coefficient 
corrr_value=np.corrcoef(np.array(x1v), np.array(y_out))
rrlai= corrr_value[0,1]
maelai=mean_absolute_error(x1v,y_out)
#Plotting
fig, ax = plt.subplots(figsize=(5, 5))    
plt.plot(x1v,y_out, 'go')
plt.xlim([0, 8])
plt.ylim([0, 8])
plt.xlabel("Observed PAI ($m^{2}~m^{-2}$)")
plt.ylabel("Estimated PAI ($m^{2}~m^{-2}$)")
#plt.title("PAI plot")
plt.plot([0, 8], [0, 8], 'k:')
plt.annotate('r = %.2f'%rrlai, xy=(0.5, 7.4))#round off upto 3decimals
plt.annotate('RMSE = %.3f'%rmselai, xy=(0.5, 6.8))
plt.annotate('MAE = %.3f'%maelai, xy=(0.5, 6.2))
plt.xticks(np.arange(0, 8+1, 2.0))
matplotlib.rcParams.update({'font.size': 24})
plt.savefig('RHRV_PAI.png',bbox_inches="tight",dpi=100)
plt.show()
plt.close()



## ---------------------------------------------------------------------------------------------
## iS-Omega SVR model-------------------------------------------------------------
Xr=np.column_stack((Psi0,Pdi0,Pvi0))
valXr=np.column_stack((Psiv,Pdiv,Pviv))

#pipelineiso = make_pipeline(StandardScaler(),SVR(kernel='rbf', epsilon=0.7, C=100, gamma = 1.5),)
pipelineiso = make_pipeline(StandardScaler(),SVR(kernel='rbf', epsilon=0.1, C=100, gamma = 1.1),)
SVRmodel=pipelineiso.fit(Xr,x1)
# Predict for validation data
#y_out1 = pipelineiso.predict(valXr);


##PAI estimation and error-------------------------------------------------------------
def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())
rmselai = rmse(np.array(x1v), np.array(y_out1))
#Correlation coefficient 
corrr_value=np.corrcoef(np.array(x1v), np.array(y_out1))
rrlai= corrr_value[0,1]
maelai=mean_absolute_error(x1v,y_out1)
#Plotting
fig, ax = plt.subplots(figsize=(5, 5))    
plt.plot(x1v,y_out1, 'go')
plt.xlim([0, 8])
plt.ylim([0, 8])
plt.xlabel("Observed PAI ($m^{2}~m^{-2}$)")
plt.ylabel("Estimated PAI ($m^{2}~m^{-2}$)")
#plt.title("PAI plot")
plt.plot([0, 8], [0, 8], 'k:')
plt.annotate('r = %.2f'%rrlai, xy=(0.5, 7.4))#round off upto 3decimals
plt.annotate('RMSE = %.3f'%rmselai, xy=(0.5, 6.8))
plt.annotate('MAE = %.3f'%maelai, xy=(0.5, 6.2))
plt.xticks(np.arange(0, 8+1, 2.0))
matplotlib.rcParams.update({'font.size': 24})
plt.savefig('iSO_PAI.png',bbox_inches="tight",dpi=100)
plt.show()
plt.close()


## ---------------------------------------------------------------------------------------------
## m-chi SVR model-------------------------------------------------------------
Xr=np.column_stack((Psm0,Pdm0,Pvm0))
valXr=np.column_stack((Psmv,Pdmv,Pvmv))

#pipelinem = make_pipeline(StandardScaler(),SVR(kernel='rbf', epsilon=0.1, C=100, gamma = 0.1),)
pipelinem = make_pipeline(StandardScaler(),SVR(kernel='rbf', epsilon=0.2, C=100, gamma = 0.2),)
SVRmodel=pipelinem.fit(Xr,x1)
# Predict for validation data
#y_out2 = pipelinem.predict(valXr);


##PAI estimation and error-------------------------------------------------------------
def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())
rmselai = rmse(np.array(x1v), np.array(y_out2))
#Correlation coefficient 
corrr_value=np.corrcoef(np.array(x1v), np.array(y_out2))
rrlai= corrr_value[0,1]
maelai=mean_absolute_error(x1v,y_out2)
#Plotting
fig, ax = plt.subplots(figsize=(5, 5))    
plt.plot(x1v,y_out2, 'go')
plt.xlim([0, 8])
plt.ylim([0, 8])
plt.xlabel("Observed PAI ($m^{2}~m^{-2}$)")
plt.ylabel("Estimated PAI ($m^{2}~m^{-2}$)")
#plt.title("PAI plot")
plt.plot([0, 8], [0, 8], 'k:')
plt.annotate('r = %.2f'%rrlai, xy=(0.5, 7.4))#round off upto 3decimals
plt.annotate('RMSE = %.3f'%rmselai, xy=(0.5, 6.8))
plt.annotate('MAE = %.3f'%maelai, xy=(0.5, 6.2))
plt.xticks(np.arange(0, 8+1, 2.0))
matplotlib.rcParams.update({'font.size': 24})
plt.savefig('mChi_PAI.png',bbox_inches="tight",dpi=100)
plt.show()
plt.close()




