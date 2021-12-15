import matplotlib.pyplot as plt
import pandas as pd

#load shooting data
prefix="startL4"
shootingDf=pd.read_csv(prefix+"shootingR12.csv")

gShooting=shootingDf["g"]
EShooting=shootingDf["E"]
################load adj data
levelNum=7
adjDf=pd.read_csv("level"+str(levelNum)+"adjGComplexR12.csv")
gAdj=adjDf["g"]
EReAdj=adjDf["ERe"]
EImAdj=adjDf["EIm"]
#############
fig, ax = plt.subplots(figsize=(20, 20))
ax.set_ylabel("E")
# plt.yscale('symlog')
ax.set_xscale("log")
ax.set_xlabel("g")
ax.set_title("Eigenvalues for potential $V(x)=\lambda x^{2}-ix^{5}$")






sctShooting=ax.scatter(gShooting,EShooting,color="blue",marker=".",s=40,label="Shooting")
# sctAdjRe=ax.scatter(gAdj,EReAdj,color="fuchsia",marker="+",s=50,label="WKB adj real")
# sctAdjIm=ax.scatter(gAdj,EImAdj,color="green",marker="^",s=50,label="WKB adj imag")

ax.hlines(y=0, xmin=1e-4, xmax=10, color='r')
plt.legend()
plt.savefig("n="+str(levelNum)+"tmp12.png")