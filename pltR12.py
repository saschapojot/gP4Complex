import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
#load shooting data
prefix="startL4"
shootingDf=pd.read_csv(prefix+"shootingR120.01deleted.csv")

gShooting=shootingDf["g"]
EShooting=shootingDf["E"]
################load adj data
levelNum=8
adjDf=pd.read_csv("level"+str(levelNum)+"adjGComplexR12.csv")
# adjDf=adjDf.drop(adjDf[np.abs(adjDf.EIm)>1].index)
gAdj=np.array(adjDf["g"])
EReAdj=np.array(adjDf["ERe"])
EImAdj=np.array(adjDf["EIm"])
#############
EMax=80
indLarge=np.logical_and(np.abs(EImAdj)<EMax,np.abs(EImAdj)>1)
indSmall=np.where(np.abs(EImAdj)<=1)[0]
gAdjLarge=gAdj[indLarge]
gAdjSmall=gAdj[indSmall]
EReAdjLarge=EReAdj[indLarge]
EReAdjSmall=EReAdj[indSmall]
EImAdjSmall=EImAdj[indSmall]
EImAdjLarge=EImAdj[indLarge]
########
fig, ax = plt.subplots(figsize=(20, 20))
ax.set_ylabel("E")
# plt.yscale('symlog')
ax.set_xscale("log")
ax.set_xlabel("g")
ax.set_title("Eigenvalues for potential $V(x)=\lambda x^{2}-ix^{5}$")






sctShooting=ax.scatter(gShooting,EShooting,color="blue",marker=".",s=40,label="Shooting")
sctAdjReSmall=ax.scatter(gAdjSmall,EReAdjSmall,color="fuchsia",marker="+",s=50,label="WKB adj real |im|<1")
sctAdjImSmall=ax.scatter(gAdjSmall,EImAdjSmall,color="green",marker="^",s=50,label="WKB adj imag  |im|<1")
sctAdjReLarge=ax.scatter(gAdjLarge,EReAdjLarge,color="brown",marker="+",s=50,label="WKB adj real |im|>1")
sctAdjImLarge=ax.scatter(gAdjLarge,EImAdjLarge,color="cyan",marker="^",label="WKB adj imag |im|>1")
ax.hlines(y=0, xmin=1e-4, xmax=10, color='r')
plt.legend()
plt.savefig("n="+str(levelNum)+"tmp12.png")