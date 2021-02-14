import numpy as np 
import pandas as pd 
import copy
import sys


def savematnew(m, name):
    strs = ""
    with open("000mat/" + name[:-3] + "matnew", "w") as f:
        for r in m:
            for l in r:
                strs = strs + str(l) + ","
            strs = strs[:-1] + '\n'
        f.write(strs)

centre = []
labels = []
with open("results.txt", "r") as f:
    centre = f.readline()
    labels = f.readline()

centre = np.array(centre.split()).astype("int32")
labels = np.array(labels.split()).astype("int32")


dfsubmatrix = pd.read_csv("submatrix.txt", header=None)
dfprotein = pd.read_csv("proteinname.txt", header=None).iloc[:,:]

centres = []

for i in range(len(centre)):
    centres.append(dfsubmatrix.iloc[centre[i],1:].values)
centres = np.array(centres)


m = 10
cha = np.zeros((100, 1))
for n in range(len(dfprotein.values)):

    mat = pd.read_csv("000mat/" + dfprotein[0][n][:-3] + "mat", header=None).values
 
    indexs = dfsubmatrix[dfsubmatrix[0]==dfprotein[0][n]].index

    l = 0
    for i in range(mat.shape[0]-m+1):
        for j in range(mat.shape[0]-m+1):
            if i==0 and j==0:#完全不重复
                mat[i:(i+m),j:(j+m)] = centres[labels[indexs[l]]].reshape(m, m)
            elif i==0 and j>0:#最后一列不重复
                temp1 = copy.deepcopy(mat[i:(i+m),j:(j+m)][:, :-1])
                temp2 = copy.deepcopy(centres[labels[indexs[l]]].reshape(m, m)[:, :-1])
                mat[i:(i+m),j:(j+m)] = centres[labels[indexs[l]]].reshape(m, m)
                mat[i:(i+m),j:(j+m)][:, :-1] = (temp1 + temp2) / 2
            elif i>0 and j==0:#最后一行不重复
                temp1 = copy.deepcopy(mat[i:(i+m),j:(j+m)][:-1, :])
                temp2 = copy.deepcopy((centres[labels[indexs[l]]].reshape(m, m)[:-1, :]))
                mat[i:(i+m),j:(j+m)] = centres[labels[indexs[l]]].reshape(m, m)
                mat[i:(i+m),j:(j+m)][:-1, :] = (temp1 + temp2) / 2
            else:#最后一个元素不重复
                x = copy.deepcopy(centres[labels[indexs[l]]].reshape(m, m)[-1, -1])
                mat[i:(i+m),j:(j+m)] = (centres[labels[indexs[l]]].reshape(m, m) + mat[i:(i+m),j:(j+m)]) / 2
                mat[i:(i+m),j:(j+m)][-1, -1] = x
            l = l + 1
    print(l==len(indexs))
    dim = mat.shape[0]
    cha[n] = np.linalg.norm(mat.reshape(1, dim**2) - pd.read_csv("000mat/" + dfprotein[0][n][:-3] + "mat", header=None).values.reshape(1, dim**2))
    print(dfprotein[0][n], cha[n])
    savematnew(mat, dfprotein[0][n])



print("\nTotal difference for all protein = " + str(np.sum(cha)))
print("\nAverage difference for all protein = " + str(np.sum(cha) / 100))

print("\nDifference for maximal protein = " + str(np.max(cha)))
print("\nDifference for minimal protein = " + str(np.min(cha)))
    
