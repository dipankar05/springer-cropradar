# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 13:32:27 2020

@author: Dipankar
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
import pandas as pd
from sklearn.metrics import mean_absolute_error
##-----------------------------------------------------------------------------

# Incidence angle
th=30
thr=th*3.1415/180

##-----------------------------------------------------------------------------
## Linear scale function-WCM
## WCM parameters are estimated with calibration data for corn
## VV model
a=0.173
b=1.025
c=0.258
d=-0.745
e=0.037
f=0.32
def fitFuncVV(x1,x2):
    return (a*(np.power(x1,e))*np.cos(thr)*(1-np.exp((-2)*b*np.power((x1),f)/np.cos(thr))))+((d*np.exp(c*x2))*np.cos(thr)*np.exp((-2)*b*np.power((x1),f)/np.cos(thr)))


## VH model
a1=0.074
b1=0.551
c1=0.23
d1=-0.025
e1=0.221
f1=0.35
def fitFuncVH(x1,x2):
    return (a1*(np.power(x1,e1))*np.cos(thr)*(1-np.exp((-2)*b1*np.power((x1),f1)/np.cos(thr))))+((d1*np.exp(c1*x2))*np.cos(thr)*np.exp((-2)*b1*np.power((x1),f1)/np.cos(thr)))


## Single-channel algorithm
##--------------------------------------------
#def Error2(X):
#    x1,x2=X
#    return np.sqrt((y1-fitFuncVV(x1,x2))*(y1-fitFuncVV(x1,x2)))
#
## Multi-channel algorithm
##--------------------------------------------
#def Error1(X):
#    x1,x2=X
#    return np.sqrt(((y1-fitFuncVV(x1,x2))*(y1-fitFuncVV(x1,x2)))+((y2-fitFuncVH(x1,x2))*(y2-fitFuncVH(x1,x2)))) 


##--------------------------------------------
########################################################################################################
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
valsm=np.float64(valdm[:,2])#Col5==Soilmoisture
#valY=np.column_stack((vallai,valsm))
valY=vallai

#valX=np.column_stack((valHH,valHV,valVH,valVV))
valX=np.column_stack((valVH,valVV))

numrows = len(valX) 
laiout = np.zeros(numrows)
smout= np.zeros(numrows)

## Iterative optimization approch for inversion
##--------------------------------------------
## Initialization    
x0 = [0.5, 0.2] #LAI, Mv

for index, row in vald.iterrows():
    
    # Observed VV and HH
    y1 = row[3]
    y2 = row[4]

    ##--------------------------------------------
    ## Multi-channel algorithm
    ##--------------------------------------------
    def Error1(X):
        x1,x2=X
        return np.sqrt(((y1-fitFuncVV(x1,x2))*(y1-fitFuncVV(x1,x2)))+((y2-fitFuncVH(x1,x2))*(y2-fitFuncVH(x1,x2)))) 
    
    ### Constrained optimization
    sol1 = least_squares(Error1,x0, method='trf', bounds=([0, 0.05],[6.0, 0.5]), ftol=1e-06,)
    
    ### Unconstrained optimization
#    sol1 = least_squares(Error1,x0, method='trf', ftol=1e-06,)
    
    laiout[index] = sol1.x[0]
    smout[index] = sol1.x[1]
    
#    laiout2[index] = sol2.x[0]
#    smout2[index] = sol2.x[1]

### Ref
##-----------------------------------------------------------------------
##https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
### method{‘trf’, ‘dogbox’, ‘lm’}, optional
### ‘trf’ : Trust Region Reflective algorithm, particularly suitable for large sparse problems with bounds
###‘dogbox’ : dogleg algorithm with rectangular trust regions, typical use case is small problems with bounds. Not recommended for problems with rank-deficient Jacobian.
### ‘lm’ : Levenberg-Marquardt algorithm as implemented in MINPACK. Doesn’t handle bounds and sparse Jacobians. Usually the most efficient method for small unconstrained problems.
###ftol--Tolerance for termination by the change of the cost function. Default is 1e-8.

## End of for loop

##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
y_out = laiout
##PAI estimation and error
#rmse estimation
def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())
rmselai = rmse(np.array(vallai), np.array(y_out))
#Correlation coefficient 
corrr_value=np.corrcoef(np.array(vallai), np.array(y_out))
rrlai= corrr_value[0,1]
maelai=mean_absolute_error(vallai,y_out)
#Plotting
plt.plot(vallai,y_out, 'go')
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

## End