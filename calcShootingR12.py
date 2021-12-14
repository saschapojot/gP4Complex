from shootingFuncsR12 import *
import numpy as np
import pandas as pd



num=500
startG=1e-4
stopG=1e1
gnIndAll=mpmath.linspace(mpmath.log10(startG), mpmath.log10(stopG), num)
gAll = [mpmath.power(10,elem) for elem in gnIndAll]
EMax=20

#convert to lambda and F
inDataAll=[] #contains [lambda, FEst]
dE=0.5
for g in gAll:
    for E in mpmath.arange(1,EMax,dE):
        lmd=mpmath.power(g,-4/7)#g**(-4/7)
        FEst=mpmath.power(E,-2/7)#E*g**(-2/7)
        inDataAll.append([lmd,FEst])

threadNum=24
tShootingStart=datetime.now()
pool1=Pool(threadNum)
retAll=pool1.map(computeOneSolution,inDataAll)
tShootingEnd=datetime.now()

print("shooting time: ",tShootingEnd-tShootingStart)
#data serialization
gShootingVals=[]
EShootingVals=[]
for itemTmp in retAll:
    if len(itemTmp)==0:
        continue

    lmd,F=itemTmp
    if F<0: continue
    if math.isnan(F):
        continue
    gTmp=mpmath.power(lmd,-7/4)#lmd**(-7/4)
    ETmp=F*mpmath.power(lmd,-1/2)#F*lmd**(-1/2)
    if mpmath.fabs(ETmp)>25:
        continue
    gShootingVals.append(gTmp)
    EShootingVals.append(ETmp)

dataDict={"g":gShootingVals,"E":EShootingVals}
dfShooting=pd.DataFrame(data=dataDict)

dfShooting.to_csv("startL"+str(LEst)+"shootingR12.csv",index=False)