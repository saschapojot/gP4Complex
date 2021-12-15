from funcsParalComplexR12 import *
import pandas as pd

num=60
startG=1e-4
stopG=1e1
gnIndAll = mpmath.linspace(mpmath.log10(startG), mpmath.log10(stopG), num)

gAll = [10 ** elem for elem in gnIndAll]


threadNum = 24
# energyLevelMax = 4
levelStart=7
levelEnd=levelStart
levelsAll = range(levelStart, levelEnd + 1)
inDataAll=[]

for nTmp in levelsAll:
    for gTmp in gAll:
        EEst=(nTmp+1/2)*mp.pi
        inDataAll.append([nTmp,gTmp,EEst])


############################################################################
######Adj computation
# ###########parallel computation  for adj, may be memory consuming
tWKBAdjStart=datetime.now()
pool1=Pool(threadNum)
retAllAdj=pool1.map(vecComputeOneSolutionWith5AdjacentPairs,inDataAll)
tWKBAdjEnd=datetime.now()
print("parallel WKB time for adj pairs: ",tWKBAdjEnd-tWKBAdjStart)

####################################end of parallel computation#########
###############################parallel computation for conjugate(E)
inDataAllConj=[]
for itemTmp in retAllAdj:
    if len(itemTmp)==0:
        continue
    n,g,E=itemTmp
    inDataAllConj.append([n,g,mpmath.conj(E)])
tConjWKBStart=datetime.now()
pool2=Pool(threadNum)
retAllAdjConj=pool2.map(vecComputeOneSolutionWith5AdjacentPairs,inDataAllConj)
tConjWKBEnd=datetime.now()
print("parallel WKB for adj conj: ",tConjWKBEnd-tConjWKBStart)

########################################################################
#data serialization for adj
nValsAdj=[]
gValsAdj=[]
ERealValsAdj=[]
EImagValsAdj=[]
for itemTmp in retAllAdj:
    if len(itemTmp)==0:
        continue
    n,g,E=itemTmp
    EReTmp=mpmath.re(E)
    EImTmp=mpmath.im(E)
    nValsAdj.append(n)
    gValsAdj.append(g)
    ERealValsAdj.append(EReTmp)
    EImagValsAdj.append(EImTmp)

############################
###########data serialization for adj conj
for itemTmp in retAllAdjConj:
    if len(itemTmp)==0:
        continue
    n,g,E=itemTmp
    EReTmp=mpmath.re(E)
    EImTmp=mpmath.im(E)
    nValsAdj.append(n)
    gValsAdj.append(g)
    ERealValsAdj.append(EReTmp)
    EImagValsAdj.append(EImTmp)
######################
#write data of adj to csv
adjDat=np.array([nValsAdj,gValsAdj,ERealValsAdj,EImagValsAdj]).T

adjDf=pd.DataFrame(adjDat,columns=["n","g","ERe","EIm"])
adjDf.to_csv("level"+str(levelStart)+"adjGComplexR12.csv",index=False)