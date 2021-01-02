# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 09:23:28 2018

@author: Administrator
"""
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
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
Y=np.column_stack((Ym0,Ym2))


XVH = pd.read_csv('VHsimulatedCorn.csv',header=None)
#XVH=np.float64(XVH.as_matrix(columns=None))
XVV = pd.read_csv('VVsimulatedCorn.csv',header=None)
#XVV=np.float64(XVV.as_matrix(columns=None))
#X=np.column_stack((XHH,XHV,XVH,XVV))
#X=np.column_stack((XHH,XHV,XVH,XVV,thr))

Y2=pd.concat([XVH,XVV,Y1],axis=1,ignore_index=True)
############################

df1=np.column_stack((Ym0,Ym2,XVH,XVV,thr))

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
valbiom=np.float64(valdm[:,1])#Col9==Wetbiomass; Col8==VWC; Col10==drybiomass
valsm=np.float64(valdm[:,2])#Col5==Soilmoisture
#valY=np.column_stack((vallai,valsm))
valY=vallai

#valX=np.column_stack((valHH,valHV,valVH,valVV))
valX=np.column_stack((valVH,valVV))

df2=vald;


numrows = len(valX) 
laiout = np.zeros(numrows)
smout= np.zeros(numrows)

for index, row in vald.iterrows():
    #m = 10000
    min = None
    for index1, row1 in Y2.iterrows():
        
        #RMSE
        rmse=np.sqrt((((row[4]-row1[0])**2)+
                   ((row[3]-row1[1])**2))/2)
#        
        #L1 estimate
#        rmse = (abs(row[4]-row1[0]) + abs(row[3]-row1[1]))
        
        #Bhattacharya distance
#        rmse = (-np.log(1 + (row[4]*row1[0])**0.5 - 0.5*(row[4] + row1[0]) +
#                (row[3]*row1[1])**0.5 - 0.5*(row[3] + row1[1])))
      
        if laiout[index] == 0 or rmse < min:
            min = rmse
            laiout[index] = row1[3]
            #smout[index] = row1[4]
        
    
 
#rmse value between datafrmae and list 
valrmselailut=((vallai - laiout) ** 2).mean() ** .5
valrrlailut=np.corrcoef(vallai, laiout)   
maelai=mean_absolute_error(vallai,laiout)



#df3 = pd.DataFrame(
#    {'PAIp': laiout
#    })
##write dataframe to excel
#writer = pd.ExcelWriter('LUTRetrievedPAIBiom.xlsx', engine='xlsxwriter')
## Convert the dataframe to an XlsxWriter Excel object.
#df3.to_excel(writer, sheet_name='Sheet1')
#
## Close the Pandas Excel writer and output the Excel file.
#writer.save()




#####Plotting PAI
#Plotting
plt.plot(vallai,laiout, 'go')
plt.xlim([0, 6])
plt.ylim([0, 6])
plt.xlabel("Observed LAI ($m^2 m^{-2}$)")
plt.ylabel("Estimated LAI ($m^2 m^{-2}$)")
plt.plot([0, 6], [0, 6], 'k:')
plt.annotate('r = %.2f'%valrrlailut[0,1], xy=(0.5, 5.5))#round off upto 3decimals
plt.annotate('RMSE = %.2f'%valrmselailut, xy=(0.5, 5.0))
plt.annotate('MAE = %.2f'%maelai, xy=(0.5, 4.5))
matplotlib.rcParams.update({'font.size': 20})
plt.yticks(np.arange(0, 7, 2))
plt.xticks(np.arange(0, 7, 2))
plt.gca().set_aspect('equal', adjustable='box')
#plt.savefig('PAIValidationLUT.png',bbox_inches="tight",dpi=100)
plt.show()

