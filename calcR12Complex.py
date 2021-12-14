from funcsParalComplexR12 import *
import pandas as pd

num=15
startG=1e-4
stopG=1e1
gnIndAll = mpmath.linspace(mpmath.log10(startG), mpmath.log10(stopG), num)

gAll = [10 ** elem for elem in gnIndAll]


threadNum = 24
# energyLevelMax = 4
levelStart=0
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
retAllAdj=pool1.map(computeOneSolutionWith5AdjacentPairs,inDataAll)
tWKBAdjEnd=datetime.now()
print("parallel WKB time for adj pairs: ",tWKBAdjEnd-tWKBAdjStart)

####################################end of parallel computation#########
########################################################################

print(retAllAdj)