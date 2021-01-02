# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 21:34:08 2020

@author: Dipankar
"""

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

## Model parameters for Pd iS-O model
#a = 0.81165971;       
#b = -0.22474;
#c = -0.17093;
#e = 1.677;
#f = 0.5110129;

a = 0.75294585  
b=-0.02804072 
c=-0.27982399  
e=1.78607574  
f=0.47209796




## Model parameters Ps iS-O model
#b1 = 0.2141
#d1 = -0.1393
#e1 = 0.021
#f1 = 0.65

b1 = 0.329;  
d1 = -0.0119;  
e1 = 0.025398; 
f1 = 0.523;



## Model parameters Pv iS-O model
a2 = 0.45217;  
b2 = 0.1439 #0.91169045;  
c2 = 0.1562; #0.06829466; 
d2 = 0.17981; #0.21571469;
      


for i in x1:
    ## Pd IS-O model
    yd=(a*1.0*np.exp((-2)*f*x1/np.cos(thr))*((b*(1-np.exp((-1)*c*x1)))+(e*(0.035*np.power(x1,2.824)))));
    yddb = 10*np.log10(yd)
    ## Ps iS-O model
    ys = (np.power(b1,x1)*(0.035*np.power(x1,2.824))*np.exp((-2)*d1*x1/np.cos(thr)))+(e1*1.0*np.exp((-2)*d1*x1/np.cos(thr))*np.exp((-2)*f1*(0.035*np.power(x1,2.824))/np.cos(thr)))
    ysdb = 10*np.log10(ys)
    ## Pv iS-O model
    yv = (1-a2)*b2*(1-np.exp((-1)*c2*x1/0.5))*np.cos(thr)*(1-np.exp((-2)*d2*x1/np.cos(thr)));
    yvdb = 10*np.log10(yv)


fig, ax = plt.subplots(figsize=(5, 5))    
plt.plot(x1,yddb,color='#FA7805',linewidth=2.5)
plt.xlim([0, 4.5])
plt.ylim([-35, 0])   
 

#fig, ax = plt.subplots(figsize=(5, 5))   
### RH simulation plots
#plt.xlim([0, 9])
#plt.ylim([-35, 0]) 
plt.plot(x1,ysdb,color="#3361FF",linewidth=2.5)  
plt.plot(x1,yvdb,color="#2BEC14",linewidth=2.5) 
#plt.plot(x1,yhsoildb)  
#plt.plot(x1,yhvegdb)  
plt.xlabel("PAI ($m^{2}~m^{-2}$)")
plt.ylabel("$iS-\Omega$ powers (dB) ")
#plt.title("RH-Rice")
plt.xticks(np.arange(0, 4.5+0.5, 1.5))
matplotlib.rcParams.update({'font.size': 24})
plt.savefig('iSOCal.png',bbox_inches="tight",dpi=100)
plt.show()
plt.close()