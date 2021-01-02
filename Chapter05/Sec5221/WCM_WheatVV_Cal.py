# -*- coding: utf-8 -*-
"""
Created on Sun Aug 20 10:53:47 2017

@author: Dipankar
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import warnings
from scipy.optimize import differential_evolution
from scipy.stats.distributions import  t

X1= pd.read_csv('calibrationcrop.csv',header=None)
Y1=pd.read_csv('calibrationC33.csv',header=None)


#Pandas dataframe to matrix conversion
y=Y1.values;
x=X1.values;

#Incidence angle
th=35
thr=th*3.1415/180
y=y[:,0]
x1=x[:,0]
x2=x[:,1]
x4=x[:,2] #SM

#def fitFunc(X,a,b,c,d,e,f):
#    x1,x2,x4=X
#    return (a*(np.power(x1,e))*(1-np.exp((-2)*b*np.power((x2),f)/np.cos(thr))))+((c+(d*x4))*np.exp((-2)*b*np.power((x2),f)/np.cos(thr)))

## Linear scale function-WCM
def fitFunc(X,a,b,c,d,e,f):
    x1,x2,x4=X
    return (a*(np.power(x1,e))*np.cos(thr)*(1-np.exp((-2)*b*np.power((x2),f)/np.cos(thr))))+((d*np.exp(c*x4))*np.cos(thr)*np.exp((-2)*b*np.power((x2),f)/np.cos(thr)))

def Error(parameterTuple):
    warnings.filterwarnings("ignore") # do not print warnings by genetic algorithm
    return np.sqrt(np.sum((y - fitFunc((x1,x2,x4), *parameterTuple)) ** 2).mean())


def generate_Initial_Parameters():
    ## min and max used for bounds
    parameterBounds = []
    parameterBounds.append([0,0.5]) # parameter bounds for a
    parameterBounds.append([-1.0,0.1]) # parameter bounds for b
    parameterBounds.append([1.0,2.0]) # parameter bounds for c
    parameterBounds.append([0,1.0]) # parameter bounds for d
    parameterBounds.append([-2.0,1.0]) # parameter bounds for e
    parameterBounds.append([-0.5,0.5]) # parameter bounds for f
   

    # "seed" the numpy random number generator for repeatable results
    #https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.optimize.differential_evolution.html
    result = differential_evolution(Error, parameterBounds, strategy='best1bin',polish=True,seed=3,init='latinhypercube')
    return result.x

###-----------------------------------------------------------------------------------------------
## generate initial parameter values
initialParameters = generate_Initial_Parameters()

###-----------------------------------------------------------------------------------------------

##-----------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------
#Final p0 best
#initialParameters=[0.45, -0.949, 1.833, 0.725, -1.727, -0.12]

##-----------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------

# curve fit the test data
fitParams,fitCovariances = curve_fit(fitFunc,(x1,x2,x4),y, initialParameters,method='lm',maxfev=4000,ftol=1e-6)

print('Fitted parameters=',fitParams)


##--------------t-test
alpha = 0.05 # 95% confidence interval = 100*(1-alpha)

n = len(y)    # number of data points
p = len(fitParams) # number of parameters

dof = max(0, n - p) # number of degrees of freedom

## student-t value for the dof and confidence level
##https://stackoverflow.com/questions/21494141/how-do-i-do-a-f-test-in-python
tval = t.ppf(1.0-alpha/2., dof) 
#
for i, p,var in zip(range(n), fitParams, np.diag(fitCovariances)):
   sigma = var**0.5
   print ('p{0}: {1} [{2}  {3}]'.format(i, p, p - sigma*tval, p + sigma*tval))


#Mutual information
# http://scikit-learn.org/stable/auto_examples/feature_selection/plot_f_test_vs_mi.html

#F-statistics tutorial
#http://www.statisticshowto.com/f-statistic/



#predicting with fitted function
A=x.T
ypred=fitFunc(A,fitParams[0],fitParams[1],fitParams[2],fitParams[3],fitParams[4],fitParams[5])

#rmse estimation
def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

#Correlation coefficient 
#corrr=np.corrcoef(ypred,y)
#rr= corrr[0,1]
#print('r=',rr)

ydb = 10.*np.log10(y)
ypreddb = 10.*np.log10(ypred)

rmse_val2 = rmse(ypreddb, ydb)
print('RMSE (dB)=',rmse_val2)

corrr=np.corrcoef(ypreddb,ydb)
rr= corrr[0,1]
print('r=',rr)

plt.scatter(ydb,ypreddb,marker='s',c='r',s=40)
plt.xlim([-30, 0])
plt.ylim([-30, 0])
plt.xlabel("Observed $\sigma^0$")
plt.ylabel("Estimated $\sigma^0$")
plt.plot([0, -30], [0, -30], 'k:')
plt.annotate('r = %.2f'%rr, xy=(-29, -3))#round off upto 3decimals
plt.annotate('RMSE = %.2fdB'%rmse_val2, xy=(-29, -6))
matplotlib.rcParams.update({'font.size': 18})
plt.yticks(np.arange(-30, 1, 10))
plt.xticks(np.arange(-30, 1, 10))
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

#end