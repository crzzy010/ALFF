#-*- coding:utf-8 -*-

import sys
import os
import numpy as np
import pandas as pd 


mulu = "C:\\Users\\Crzzy\\Desktop\\data\\matix"


def euclidean_distance(a, b):
    d = 0                       
    if a.ndim > 1 or b.ndim >1:     #如果是多位数组，向量操作
        d = 1
    return np.sqrt(np.sum(np.power(a-b, 2), axis=d))


def creat_sub_matrix(matrix, m, filename, proteinname):  #分割出满足条件的mxm矩阵
    for i in range(matrix.shape[0]-m):
        for j in range(matrix.shape[0]-m):
            sub = matrix[i:(i+m),j:(j+m)].reshape(1,m*m)[0]
            line = proteinname+','
            for o in range(len(sub)):
                line += str(sub[o])
                if o < len(sub) - 1:
                    line += ',' 
            line += "\n"
            filename.write(line)
    return

def start(dataname, filename):
    df = pd.read_csv(dataname)
    filelist = []
    with open(filename) as f:
        filelist = f.readlines()  #逐行读文本
    with open("submatrix.txt", 'a') as fma:
        for p in filelist:
            if len(df[df.name==p.strip()]) > 0:
                protein = np.array(df[df.name==p.strip()][['x','y','z']])
                m = np.zeros((len(protein), len(protein)))
                for i in range(len(protein)):
                    m[i] = euclidean_distance(protein, protein[i])
                creat_sub_matrix(m, 10, fma, p.strip())
    return




if __name__ == "__main__":
    start("result_aa.txt","proteinname.txt")
