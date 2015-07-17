import pandas as pd
import math as mt
import numpy as np
import csv as csv
import scipy as sp
from scipy import integrate
from __future__ import print_function, division

pre=open("./predata2.dat", "w+")
print (" ","r"," ","g",file=pre)
def drange(start , stop , step):
    r = start
    while r < stop:
        print(" ",r,"  0.0","\n",file=pre)
        r +=step
i0=drange(0.0,1.32,.02) #change middle value to match r value in original rdf file
rdf = open("hw_hw_rdf", "r") 
for line in rdf:
    print (line,file=pre)
pre.close()
rdf.close()

data= pd.read_csv("predata2.dat", sep="  ")
rows =len(data)-1

count1=.001
count2=0
struct=open('hw_hw_struc.test', 'w').close()
while True:
    struct=open('hw_hw_struc.test', 'ab')
    preintfile=open ("hw_hw_preint.test", "wb+")
    while True:
        preint = 4*3.14159*766/count1/(19.5970437582**3)*(data.ix[count2]["r"])*((data.ix[count2][" g"])-1)*mt.sin(count1*(data.ix[count2]["r"]))
        print (data.ix[count2]["r"], preint,file=preintfile)
        if count2 ==rows:
            break
        count2 +=1 
    preintfile.close()
    
    data2 = np.loadtxt("hw_hw_preint.test")
    [r, integrand] = np.hsplit(data2, 2)
    r = r.flatten()
    integrand = integrand.flatten()
    Integral = sp.integrate.cumtrapz(integrand,r)
    cum_sum = Integral[len(Integral)-1]
    struct_fact = cum_sum
    print (count1, '\t', struct_fact, file=struct)
    
    count2=0
    if count1 == .005:
        break
    count1 +=.001
