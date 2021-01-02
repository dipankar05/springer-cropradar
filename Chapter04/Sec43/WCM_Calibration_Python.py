#------------------------------------------------------------------------------
## Water Cloud Model calibration for Soybean crop
#------------------------------------------------------------------------------
"""
Version: Python 3.7.4
@author: Dipankar
"""

"""
Spyder Editor

"""
## import libraries
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import warnings
from scipy.optimize import differential_evolution

## import data for calibration
X1= pd.read_csv('calibrationcrop.csv',header=None); ## in-situ measured LAI and soil-moisture
Y1=pd.read_csv('calibrationHH.csv',header=None); ## SAR measured backscatter intensities in HH channel

#Pandas dataframe to matrix conversion
y=Y1.values;
x=X1.values;

#Incidence angle
th=35;
thr=th*3.1415/180;
y=y[:,0];
x1=x[:,0];
x2=x[:,1];


#------------------------------------------------------------------------------
## Linear scale function-WCM
def fitFunc(X,a,b,c,d,e):
    x1,x2=X
    return (a*(np.power(x1,e))*np.cos(thr)*(1-np.exp((-2)*b*np.power((x1),1)/np.cos(thr))))+((d*np.exp(c*x2))*np.cos(thr)*np.exp((-2)*b*np.power((x1),1)/np.cos(thr)))

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## function for genetic algorithm to minimize (RMSE error)
## bounds on parameters are set in generate_Initial_Parameters() below
## genetic algorithm for initial parameter estimation.
def Error(parameterTuple):
    warnings.filterwarnings("ignore") # do not print warnings by genetic algorithm
    return np.sqrt(np.sum((y - fitFunc((x1,x2), *parameterTuple)) ** 2).mean())


def generate_Initial_Parameters():
    ## min and max used for bounds
    parameterBounds = []
    parameterBounds.append([0,1.1]) # parameter bounds for a
    parameterBounds.append([0,0.5]) # parameter bounds for b
    parameterBounds.append([-0.5,1]) # parameter bounds for c
    parameterBounds.append([-0.5,1]) # parameter bounds for d
    parameterBounds.append([-1.5,1]) # parameter bounds for e
    ##parameterBounds.append([-100,100]) # parameter bounds for f
   

    ## "seed" the numpy random number generator for repeatable results
    ##https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.optimize.differential_evolution.html
    result = differential_evolution(Error, parameterBounds, strategy='best1bin',polish=True,seed=3,init='latinhypercube')
    return result.x

## generate initial parameter values
initialParameters = generate_Initial_Parameters()

##-----------------------------------------------------------------------------
## OR directly define initial parameters
#initialParameters=[0.2,1.357,2,4,-1.965]

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
# curve fit the test data
fitParams,fitCovariances = curve_fit(fitFunc,(x1,x2),y, initialParameters,method='lm',maxfev=6000,ftol=1e-8)


##-----------------------------------------------------------------------------

#predicting with fitted function
A=x.T
ypred=fitFunc(A,fitParams[0],fitParams[1],fitParams[2],fitParams[3],fitParams[4])

#rmse estimation
def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

rmse_val = rmse(ypred, y)
print('RMSE=',rmse_val)

#Correlation coefficient 
corrr=np.corrcoef(ypred,y)
rr= corrr[0,1]
print('R=',rr)
plt.scatter(y,ypred)
plt.xlim([0, 0.3])
plt.ylim([0, 0.3])
plt.xlabel("Observed $\sigma^0$")
plt.ylabel("Estimated $\sigma^0$")
plt.title("HH-Soybean")
plt.plot([0, 0.3], [0, 0.3], 'k:')
plt.annotate('r = %.3f'%rr, xy=(0.015,0.27))#round off upto 3decimals
plt.annotate('RMSE = %.3f'%rmse_val, xy=(0.015,0.24))
matplotlib.rcParams.update({'font.size': 20})
plt.show()
plt.savefig('HH_Soybean.png')

print('Fitted WCM coefficients =\n',fitParams)

## end