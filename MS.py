# -*- coding=utf-8 -*-
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import sys
import re
import os
import copy
import random
from sklearn.externals import joblib
import time

from sklearn import cluster, datasets, mixture
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler
from itertools import cycle, islice

mulu = "/home/crz/protein/protein_structure/data/"


def euclidean_distance(a, b):
                      
    if a.ndim > 1 or b.ndim >1:     #如果是多位数组，向量操作
        return np.linalg.norm(a-b, axis=1)
    else:
        return np.linalg.norm(a-b)


#################################################################################################################

def check_chinese(str):
    return u'\u4e00' <= str <= u'\u9fa5'


def is_float1(s):
    value = re.compile(r'^[-+]?[0-9]+\.[0-9]+$')  # 定义正则表达式
    result = value.match(s)
    return result
    
def is_float(s):
    if(len(s.strip().split('.'))==2)and is_float1(s):
      return len(s.strip().split('.')[1])
    else:
      return 0


def countCA(path, name):

    with open(path+"/"+name, "r") as f:
        Cnum = 0
        while True:
            line = f.readline()  #逐行读文本
            if not line:
                break
            else:
                strline = str(line).strip().split()

                if len(strline) > 0 and strline[0] == "ATOM" and strline[2] == "CA":
                    
                    for i in range(5):
                        # if i==4:
                        #     if num1==True:
                        #         #print("出现异常")
                        #         num2 = 1
                        #         #xun=False
                        #         break
                        if is_float(str(strline[i+4]))==3:
                            if is_float(str(strline[i+5]))==3:
                                if is_float(str(strline[i+6]))==3:
                                    Cnum = Cnum + 1

        return Cnum




def travel(path):

    if os.path.exists(os.getcwd()+"/result_aa.txt"):
       os.remove(os.getcwd()+"/result_aa.txt")
    if os.path.exists(os.getcwd()+"/proteinname.txt"):
       os.remove(os.getcwd()+"/proteinname.txt")

    filenames = []
    mo = joblib.load("filename.list")

    l = 0

    while len(filenames) < 100:
        ccnum = countCA(path, mo[l])
        if ccnum > 50 and ccnum < 150:

            filenames.append(mo[l])
        l += 1

    # for i in filenames:
    #     print(i, countCA(path, i))

    # sys.exit()
    with open("proteinname.txt",'w') as pn:
        with open("result_aa.txt",'w') as ra:
            ra.write("name,type,x,y,z\n")
            for f in filenames:
                getdata(path, f, ra)
                pn.write(f+"\n")



#生成数据之前先删除文件
def getdata(path, name, ra):

    with open(path + "/" + name, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            else:
                strline = str(line).strip().split()
                if len(strline) > 0 and strline[0] == "ATOM" and strline[2] == "CA":
                    for i in range(5):
                        if is_float(str(strline[i + 4])) == 3:
                            if is_float(str(strline[i + 5])) == 3:
                                if is_float(str(strline[i + 6])) == 3:
                                    ra.write(
                                        name + "," + strline[2] + "," + strline[i + 4] + "," +
                                        strline[i + 5] + "," +
                                        strline[i + 6] + "\n")

    # nos = nos + 1
    return

#################################################################################################################






def creat_sub_matrix(matrix, m, filename, proteinname):  #分割出满足条件的mxm矩阵
    for i in range(matrix.shape[0]-m+1):
        for j in range(i, matrix.shape[0]-m+1):
            sub = matrix[i:(i+m),j:(j+m)].reshape(1,m*m)[0]
            line = proteinname+','
            
            for o in range(len(sub)):
                line += str(sub[o])
                if o < len(sub) - 1:
                    line += ',' 
            line += "\n"
            filename.write(line)
    return

def savemat(savepath, m, name):
    strs = ""
    with open(savepath + "/000mat/" + name[:-3] + "mat", "w") as f:
        for r in m:
            for l in r:
                strs = strs + str(l) + ","
            strs = strs[:-1] + '\n'
        f.write(strs)

def start(dataname, filename, m):
    ans = 0
    df = pd.read_csv(dataname)
    filelist = []
    with open(filename) as f:
        filelist = f.readlines()  #逐行读文本
    
    X = []

    if os.path.exists(os.getcwd()+"/submatrix.txt"):
       os.remove(os.getcwd()+"/submatrix.txt")

    with open("submatrix.txt", 'a') as fma:
        for p in filelist:
            if len(df[df.name==p.strip()]) > 0:
                protein = np.array(df[df.name==p.strip()][['x','y','z']])
                # if len(protein) > 200 or len(protein) < 15:
                #     ans += 1
                #     print("######## exit ########", p, len(protein))
                #     continue
                mat = np.zeros((len(protein), len(protein)))
                for i in range(len(protein)):
                    mat[i] = euclidean_distance(protein, protein[i])
                # print(p.strip(), m.shape)
                creat_sub_matrix(mat, m, fma, p.strip())
                savemat(os.getcwd()+"/data", mat, p.strip())
                X.append(mat)
                #print(p.strip(), m.shape)
    # print("error ent:", ans)
    return np.array(X)




#################################################################################################################
def get_costs(medoids, data, fdist=euclidean_distance):
    costs = []
    for m in medoids:
        c = fdist(data, m)
        costs.append(c)
    costs = np.array(costs)

    return np.sum(np.min(costs, axis=0))

# data (m, n)
def k_medoids(data, k, fdist=euclidean_distance, maxiters=50):
    
    #Initialize
    medoids = []
    centre = []     #index of medoids
    
    centre = random.sample(range(data.shape[0]), k)
    # while len(centre) < k:
    #     r = np.random.randint(0, data.shape[0])
    #     if r not in centre:
    #         centre.append(r)
    #         medoids.append(data[r])

    centre = np.array(centre)
    medoids = data[centre]
    cluster = np.zeros(data.shape[0])


    #PAM
    convergence = False
    no = 0
    while not convergence and no < maxiters:
    #while not no > 3:
        no = no + 1
        convergence = True
        # Associate each data point to the closest medoid.
        for p in range(k):
            costs = get_costs(medoids, data, fdist)
            medoids_temp = medoids.copy()
            for m in range(data.shape[0]):
                if m not in centre:
                    medoids_temp[p] = data[m]
                    cost = get_costs(medoids_temp, data, fdist)
                    if cost < costs:        #如果交换之后代价比原来小，就交换
                        costs = cost
                        medoids[p] = data[m]
                        centre[p] = m
                        convergence = False #簇心改变了，说明算法没有收敛
                    else:                   #否则取消交换
                        medoids_temp[p] = medoids[p]

    for i in range(data.shape[0]):

        diss = fdist(data[i], medoids)
        cluster[i] = np.argmin(diss)
        # min_dist = np.inf
        # for j in range(k):
        #     dist = fdist(data[i], medoids[j])
        #     if dist < min_dist:
        #         min_dist = dist
        #         cluster[i] = j
    print("---kmedoids iters times:"+str(no)+"  ---")
    return medoids, centre.astype("int32"), cluster.astype("int32")






def clara(data, k, fdist=euclidean_distance, maxiters=50):
    
    nums = 1
    N = 40 + 2*k
    mincost = np.inf
    medoids, centre, cluster = [], [], np.zeros(data.shape[0])
    for i in range(nums):
        sub_data_index = np.array(random.sample(range(data.shape[0]), N))
        # print(sub_data_index)
        sub_medoids, sub_centre, sub_cluster = k_medoids(data[sub_data_index], k)
        curcost = get_costs(sub_medoids, data)
        if  curcost< mincost:
            mincost = curcost
            medoids, centre = sub_medoids.copy(), sub_data_index[sub_centre]
    

    for i in range(data.shape[0]):

        diss = fdist(data[i], medoids)
        cluster[i] = np.argmin(diss)

    return medoids, centre, cluster.astype("int32")






def test_kmedoids():

    r1 = np.random.normal(100,10,size = (1000,2))

    r2 = np.random.normal(80,10,size = (1000,2))

    r = np.append(r1, r2).reshape((2000, 2))

    medoids, centre, cluster = k_medoids(r, 2)

    plt.scatter(r[cluster==0][:, 0],r[cluster==0][:, 1])

    plt.scatter(r[cluster==1][:, 0],r[cluster==1][:, 1])

    print(r[cluster==0])

    print(r[cluster==1])

    plt.show()




    
def test_clara():

    r1 = np.random.normal(100,10,size = (10000,2))

    r2 = np.random.normal(80,10,size = (10000,2))

    r = np.append(r1, r2).reshape((20000, 2))

    medoids, centre, cluster = clara(r, 2)

    plt.scatter(r[cluster==0][:, 0],r[cluster==0][:, 1])

    plt.scatter(r[cluster==1][:, 0],r[cluster==1][:, 1])

    print(r[cluster==0])

    print(r[cluster==1])

    plt.show()





#################################################################################################################


def savematnew(m, name):
    strs = ""
    with open("./data/000mat/" + name[:-3] + "matnew", "w") as f:
        for r in m:
            for l in r:
                strs = strs + str(l) + ","
            strs = strs[:-1] + '\n'
        f.write(strs)

def restructure(dfsubmatrix, dfprotein, centre=[], labels=[], m = 10, k=100):

    # centre = []
    # labels = []
    # with open("results.txt", "r") as f:
    #     centre = f.readline()
    #     labels = f.readline()

    # centre = np.array(centre.split()).astype("int32")
    # labels = np.array(labels.split()).astype("int32")


    # dfsubmatrix = pd.read_csv("submatrix.txt", header=None)
    # dfprotein = pd.read_csv("proteinname.txt", header=None).iloc[:,:]

    centres = []

    for i in range(len(centre)):
        centres.append(dfsubmatrix.iloc[centre[i],1:].values)
    centres = np.array(centres)

    print("centres_shape={}".format(centres.shape))

    cha = np.zeros((100, 1))
    for n in range(len(dfprotein.values)):

        mat = pd.read_csv("./data/000mat/" + dfprotein[0][n][:-3] + "mat", header=None).values
    
        indexs = dfsubmatrix[dfsubmatrix[0]==dfprotein[0][n]].index

        l = 0
        for i in range(mat.shape[0]-m+1):
            for j in range(i, mat.shape[0]-m+1):
                if i==0 and j==0:#完全不重复
                    print("l={}, indexs[l]={}, labels[indexs[l]={}, centres[labels[indexs[l]]]={}]".format(l, indexs[l], labels[indexs[l]], centres[labels[indexs[l]]]))
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
        # print(l==len(indexs))
        dim = mat.shape[0]
        cha[n] = np.linalg.norm(mat.reshape(1, dim**2) - pd.read_csv("./data/000mat/" + dfprotein[0][n][:-3] + "mat", header=None).values.reshape(1, dim**2))
        # print(dfprotein[0][n], cha[n])
        savematnew(mat, dfprotein[0][n])



    # print("\nTotal difference for all protein = " + str(np.sum(cha)))
    # print("\nAverage difference for all protein = " + str(np.sum(cha) / k))

    # print("\nDifference for maximal protein = " + str(np.max(cha)))
    # print("\nDifference for minimal protein = " + str(np.min(cha)))

    return np.sum(cha)




#################################################################################################################


if __name__ == "__main__":

    # 2019-08-12 MeanShift
    #################################################################################################################
    
    t0 = time.time()
    np.random.seed(0)

    M = 10


    # X = start("result_aa.txt","proteinname.txt", M)

    dfprotein = pd.read_csv("proteinname.txt", header=None).iloc[:,:]


    dfsubmatrix = pd.read_csv("submatrix.txt", header=None)


    X = dfsubmatrix.iloc[:,1:].values[:]


    # normalize dataset for easier parameter selection
    # X = StandardScaler().fit_transform(X)

    
    bandwidth = cluster.estimate_bandwidth(X, quantile=0.2)

    
    print(bandwidth)

    with open("res.txt", 'w') as f:
        f.write("bandwidth="+str(bandwidth))

    sys.exit()
    MS = cluster.MeanShift(bandwidth=bandwidth).fit(X)



    print(len(MS.cluster_centers_))

    print("loss: ", restructure(dfsubmatrix, dfprotein, MS.cluster_centers_, MS.labels_, m=M, k=len(MS.cluster_centers_)))


    bandwidth = 2

    MS = cluster.MeanShift(bandwidth=bandwidth).fit(X)



    print(len(MS.cluster_centers_))

    print("loss: ", restructure(dfsubmatrix, dfprotein, MS.cluster_centers_, MS.labels_, m=M, k=len(MS.cluster_centers_)))
    # medoids, centre, cluster = clara(X, 100, euclidean_distance)
    # print("clara: ", restructure(dfsubmatrix, dfprotein, centre, cluster, m=M, k=10))


    bandwidth = 3

    MS = cluster.MeanShift(bandwidth=bandwidth).fit(X)



    print(len(MS.cluster_centers_))

    print("loss: ", restructure(dfsubmatrix, dfprotein, MS.cluster_centers_, MS.labels_, m=M, k=len(MS.cluster_centers_)))



    bandwidth = 4

    MS = cluster.MeanShift(bandwidth=bandwidth).fit(X)



    print(len(MS.cluster_centers_))

    print("loss: ", restructure(dfsubmatrix, dfprotein, MS.cluster_centers_, MS.labels_, m=M, k=len(MS.cluster_centers_)))



    print(time.time()-t0)
    sys.exit()    


    """
    bandwidth=2
    n=18679


---kmedoids iters times:4  ---
clara:  20877.354838513154
679.8745038509369


    """

    # 2019-07-02 DBscan
    #################################################################################################################
    
    t0 = time.time()
    np.random.seed(0)


    dfprotein = pd.read_csv("proteinname.txt", header=None).iloc[:,:]


    dfsubmatrix = pd.read_csv("submatrix.txt", header=None)


    X = dfsubmatrix.iloc[:,1:].values[:]


    # normalize dataset for easier parameter selection
    # X = StandardScaler().fit_transform(X)

    # myeps = 2

    # dbscan = cluster.DBSCAN(eps=myeps).fit(X)


    # print(dbscan.core_sample_indices_, len(dbscan.core_sample_indices_), set(dbscan.labels_))




    print(time.time()-t0)
    sys.exit()

    # 2019-06-17
    #################################################################################################################
    # travel(os.getcwd()+"/data")
    # sys.exit()
    res = []
    M = 10
    dfprotein = pd.read_csv("proteinname.txt", header=None).iloc[:,:]
    start("result_aa.txt","proteinname.txt", M)
    for K in range(50, 60):
        print("k = ", K)
        dfsubmatrix = pd.read_csv("submatrix.txt", header=None)
        
        X = dfsubmatrix.iloc[:,1:].values[:]

        medoids, centre, cluster = clara(X, K, euclidean_distance)

        # dunn
        #############        
        dist = 0
        for i in range(K-1):
            dist = np.inf
            for j in range(i+1, K):
                mydist = euclidean_distance(centre[i], centre[j])
                if dist > mydist:
                    dist = mydist
        fenzi = dist
        maxdist = -1
        for o in range(K):
            dd = X[cluster==o]
            dist = -1
            for i in range(dd.shape[0]-1):
                for j in range(i+1, dd.shape[0]):
                    mydist = euclidean_distance(dd[i], dd[j])
                    if mydist > dist:
                        dist = mydist
            if dist > maxdist:
                maxdist = dist
        fenmu = maxdist
        #############

        print("dunn:", fenzi/fenmu)
        res.append(fenzi/fenmu)
    joblib.dump(res, "dunn.m")
    sys.exit()


    #################################################################################################################
    
    
    
    #################################################################################################################
    # travel(os.getcwd()+"/data")
    # sys.exit()
    res = []
    dfprotein = pd.read_csv("proteinname.txt", header=None).iloc[:,:]
    for M in range(2, 40):
        print("m = ", M)
        start("result_aa.txt","proteinname.txt", M)
        print("start f")
        dfsubmatrix = pd.read_csv("submatrix.txt", header=None)
        
        X = dfsubmatrix.iloc[:,1:].values[:]

        medoids, centre, cluster = clara(X, 100, euclidean_distance)
        print("clara f")
        res.append(restructure(dfsubmatrix, dfprotein, centre, cluster, m=M))
        print("restructure f")
    joblib.dump(res, "res.m")
    print("res = ", res)
    sys.exit()


    #################################################################################################################
    # travel(os.getcwd()+"/data")

    M = 2
    

    # X = start("result_aa.txt","proteinname.txt", M)


    dfsubmatrix = pd.read_csv("submatrix.txt", header=None)
    dfprotein = pd.read_csv("proteinname.txt", header=None).iloc[:,:]

    X = dfsubmatrix.iloc[:,1:].values[:]

    medoids, centre, cluster = clara(X, 100, euclidean_distance)

    restructure(dfsubmatrix, dfprotein, centre, cluster, m=M)

    sys.exit()
    
    res = []

    for k in range(10, 240, 3):

        medoids, centre, cluster = clara(X, k, euclidean_distance)

        res.append(restructure(dfsubmatrix, dfprotein, centre, cluster, m=M))

    joblib.dump(res, "res.m")
    # print(centre.shape) 75643



    # test_clara()





# 18
# Total difference for all protein = 3263.45995603981

# Average difference for all protein = 32.6345995603981

# Difference for maximal protein = 672.3323995328003

# Difference for minimal protein = 0.0


# Total difference for all protein = 2464.048118381719

# Average difference for all protein = 24.64048118381719

# Difference for maximal protein = 489.1272104500558

# Difference for minimal protein = 0.0


# Total difference for all protein = 2088.7476460312328

# Average difference for all protein = 20.887476460312328

# Difference for maximal protein = 404.75068171718095

# Difference for minimal protein = 0.0



# 10
# Total difference for all protein = 2304.9187928932774

# Average difference for all protein = 23.049187928932774

# Difference for maximal protein = 483.54118711614757

# Difference for minimal protein = 0.0


# 9
# Total difference for all protein = 2121.6129994671364

# Average difference for all protein = 21.216129994671363

# Difference for maximal protein = 443.3881006229356

# Difference for minimal protein = 0.0



# k=200
# Total difference for all protein = 2196.489055823164

# Average difference for all protein = 21.96489055823164

# Difference for maximal protein = 421.49631130837963

# Difference for minimal protein = 0.0

#k=250
# Total difference for all protein = 2088.7476460312328

# Average difference for all protein = 20.887476460312328

# Difference for maximal protein = 404.75068171718095

# Difference for minimal protein = 0.0