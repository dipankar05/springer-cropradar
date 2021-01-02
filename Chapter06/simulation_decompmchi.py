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
#a = 0.40122218 
#b = -0.43592192  
#c = -0.06592063  
#e = 0.84668574  
#f = 0.28204715

a = 0.40304;       
b = -0.835837;
c = -0.03657;
e = 0.50666;
f = 0.20754;


## Model parameters Ps iS-O model
b1 = 0.12403
d1 = -0.3952
e1 = 0.01388
f1 = 0.676
#b1 = 0.30580876;  
#d1 = -0.1251497149;  
#e1 = 0.012258803; 
#f1 = 0.615;

## Model parameters Pv iS-O model
#a2 = 0.61624;  
#b2 = 0.3143 #0.91169045;  
#c2 = 0.3343; #0.06829466; 
#d2 = 0.28659; #0.21571469;
      
a2=2.55455746 
b2=-0.3198832  
c2=0.04291485 
d2=0.3707592



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
plt.xlim([0,4.5])
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
plt.ylabel("$m-\chi$ powers (dB) ")
#plt.title("RH-Rice")
plt.xticks(np.arange(0, 4.5+0.5,1.5))
matplotlib.rcParams.update({'font.size': 24})
plt.savefig('mchiCal.png',bbox_inches="tight",dpi=100)
plt.show()
plt.close()