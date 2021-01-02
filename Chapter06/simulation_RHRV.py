# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 17:24:34 2020

@author: Dipankar
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# a = np.arange(0,9,0.1)
x1 = np.linspace(0.01,4.5,100)
thr = 31*np.pi/180;
         
## Model parameters for RV
#a = 0.211653;
#b = 0.53534114;
#c = -0.0463636;
#d = -0.0580817;
#e = 0.044755954;
#f = 0.543188;

a = 0.1159
b = 0.34974 
c = -0.03646
d = -0.051702
e = 0.0740385
f = 0.65813

##--------------------------------------------------------
## Model parameters for RH
#a1 = 0.146145538;
#b1 = 0.10866143;
#c1 = -0.05348105;
#d1 = -0.06477025;
#e1 = 0.1915623;
#f1 = 0.59271111;

a1 = 0.19588056 
b1= 0.074842796 
c1= -0.0153009 
d1 = -0.261887401 
e1 = 0.1837192  
f1 = 0.51333463
      


for i in x1:
    ## RV model
    y=(a*(np.power(x1,e))*np.cos(thr)*(1-np.exp((-2)*b*np.power((0.035*np.power(x1,2.824)),f)/np.cos(thr))))+(c*(d*np.power(10,1.0))*np.exp((-2)*b*np.power((0.035*np.power(x1,2.824)),f)/np.cos(thr)));
    yveg = (a*(np.power(x1,e))*np.cos(thr)*(1-np.exp((-2)*b*np.power((0.035*np.power(x1,2.824)),f)/np.cos(thr))));
    yvegdb = 10*np.log10(yveg);
    ysoil = (c*(d*np.power(10,1.0))*np.exp((-2)*b*np.power((0.035*np.power(x1,2.824)),f)/np.cos(thr)));
    ysoildb = 10*np.log10(ysoil)
    ydb = 10*np.log10(y)
    ## RH model
    yh=(a1*(np.power(x1,e1))*np.cos(thr)*(1-np.exp((-2)*b1*np.power((0.035*np.power(x1,2.824)),f1)/np.cos(thr))))+(c1*(d1*np.power(10,1.0))*np.exp((-2)*b1*np.power((0.035*np.power(x1,2.824)),f1)/np.cos(thr)));
    yhdb = 10*np.log10(yh)
    yhsoil = (c1*(d1*np.power(10,1.0))*np.exp((-2)*b1*np.power((0.035*np.power(x1,2.824)),f1)/np.cos(thr)));
    yhsoildb = 10*np.log10(yhsoil)
    yhveg=(a1*(np.power(x1,e1))*np.cos(thr)*(1-np.exp((-2)*b1*np.power((0.035*np.power(x1,2.824)),f1)/np.cos(thr))));
    yhvegdb = 10*np.log10(yhveg)

fig, ax = plt.subplots(figsize=(5, 5))    
plt.plot(x1,ydb,'-k',label='$\sigma^{0}_{total}$',linewidth=2.5)
plt.xlim([0, 4.5])
plt.ylim([-35, 0])   
plt.xticks(np.arange(0, 4.5+0.5, 1.5))
plt.plot(x1,ysoildb,'-r',label='$\sigma^{0}_{soil}$',linewidth=2.5)  
#plt.legend()
plt.plot(x1,yvegdb,'-g',label='$\sigma^{0}_{veg}$',linewidth=2.5) 
##
#
plt.xlabel("PAI ($m^{2}~m^{-2}$)")
plt.ylabel("$\sigma^{\circ}_{RV}$ (dB)")
#plt.gca().legend(ncol=3)
#plt.title("RH-Rice")
matplotlib.rcParams.update({'font.size': 24})
plt.savefig('RVCal.png',bbox_inches="tight",dpi=100)
plt.show()
plt.close()





fig, ax = plt.subplots(figsize=(5, 5))   
## RH simulation plots
plt.xlim([0, 4.5])
plt.ylim([-35, 0])   
plt.xticks(np.arange(0, 4.5+0.5, 1.5))
plt.plot(x1,yhdb,'-k',linewidth=2.5)  
plt.plot(x1,yhsoildb,'-r',linewidth=2.5)  
plt.plot(x1,yhvegdb,'-g',linewidth=2.5) 
## 
plt.xlabel("PAI ($m^{2}~m^{-2}$)")
plt.ylabel("$\sigma^{\circ}_{RH}$ (dB)")
#plt.title("RH-Rice")
matplotlib.rcParams.update({'font.size': 24})
plt.savefig('RHCal.png',bbox_inches="tight",dpi=100)
plt.show()
plt.close() 