import matplotlib.pyplot as plt
import pandas as pd

#load shooting data
prefix="startL4"
shootingDf=pd.read_csv(prefix+"shootingR12.csv")

gShooting=shootingDf["g"]
EShooting=shootingDf["E"]
##############load adj data frames
levelFirst=0
levelLast=7
# adjDfList=[]
gAdjAll=[]
EReAdjAll=[]
for j in range(levelFirst,levelLast+1):
    adjDfTmp=pd.read_csv("level"+str(j)+"adjGComplexR12.csv")
    gAdjAll.extend(list(adjDfTmp["g"]))
    EReAdjAll.extend(list(adjDfTmp["ERe"]))


#############
fig, ax = plt.subplots(figsize=(20, 20))
ax.set_ylabel("E")
# plt.yscale('symlog')
ax.set_xscale("log")
ax.set_xlabel("g")
ax.set_title("Eigenvalues for potential $V(x)=\lambda x^{2}-ix^{5}$")


sctShooting=ax.scatter(gShooting,EShooting,color="blue",marker=".",s=40,label="Shooting")
sctAdjReAll=ax.scatter(gAdjAll,EReAdjAll,color="fuchsia",marker="+",s=50,label="WKB adj real")
plt.legend()
plt.savefig("all"+"tmp12.png")
plt.close()

fig1,ax1=plt.subplots(figsize=(20,20))
ax1.set_ylabel("E")
# plt.yscale('symlog')
ax1.set_xscale("log")
ax1.set_xlabel("g")
ax1.set_title("Eigenvalues for potential $V(x)=\lambda x^{2}-ix^{5}$")
sctAdjReAll1=ax1.scatter(gAdjAll,EReAdjAll,color="fuchsia",marker="+",s=50,label="WKB adj real")
plt.legend()
plt.savefig("real"+"tmp12.png")
plt.close()