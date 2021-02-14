# -*- coding=utf-8 -*-
import numpy as np 
import matplotlib.pyplot as plt 

'''

'''
def creat_data_auto(k, m, n):   #成成随机数据
    data = np.zeros((k,m,n))
    for i in range(k):          #生成K份分布不一样的数据
        kk = np.array([2.5*i for j in range(n)])
        #print(kk)
        data[i] = np.random.normal(loc=kk, size=(m, n))

    return data


def draw(data, fp):     #只能画2维的图像
    
    ax1 = fp.add_subplot(121)
    for i in range(data.shape[0]):
        ax1.scatter(data[i,:,0],data[i,:,1])
    #fp.show()


def euclidean_distance(a, b):
    d = 0                       
    if a.ndim > 1 or b.ndim >1:     #如果是多位数组，向量操作
        d = 1
    return np.sqrt(np.sum(np.power(a-b, 2), axis=d))



def get_costs(medoids, data, fdist=euclidean_distance):
    costs = []
    for m in medoids:
        c = fdist(data, m)
        costs.append(c)
    costs = np.array(costs)

    return np.sum(np.min(costs, axis=1))

# data (m, n)
def k_medoids(data, k, fdist=euclidean_distance):
    #Initialize
    medoids = []
    centre = []     #index of medoids
    while len(centre) < k:
        r = np.random.randint(0, data.shape[0])
        if r not in centre:
            centre.append(r)
            medoids.append(data[r])
    medoids = np.array(medoids) 
    cluster = np.zeros(data.shape[0])


    #PAM
    convergence = False
    no = 0
    while not convergence:
    #while not no > 3:
        no = no + 1
        convergence = True
        # Associate each data point to the closest medoid.


        for p in range(k):
            costs = get_costs(data, medoids, fdist)
            temp = medoids.copy()
            for m in range(data.shape[0]):
                if m not in centre:
                    temp[p] = data[m]
                    cost = get_costs(data, temp, fdist)
                    if cost < costs:        #如果交换之后代价比原来小，就交换
                        costs = cost
                        medoids[p] = data[m]
                        centre[p] = m
                        convergence = False #簇心改变了，说明算法没有收敛
                    else:                   #否则取消交换
                        temp[p] = medoids[p]

    for i in range(data.shape[0]):
        min_dist = np.inf
        for j in range(k):
            dist = fdist(data[i], medoids[j])
            if dist < min_dist:
                min_dist = dist
                cluster[i] = j
    print("---time:"+str(no)+"  ---")
    return medoids, centre, cluster




def choice(data,u):
    c = []          #存u序号
    for n in data:
        m = 0
        value = np.power(n-u[0],2).sum()
        #print(value)
        for i in range(len(u)):
            if value > np.power(n-u[i],2).sum():
                m = i
                value = np.power(n-u[i],2).sum()
        c.append(m)
    #print(c)
    return np.array(c)


def updateU(data,c,u):
    

    for i in range(len(u)):
        fenmu = 0
        x = np.zeros((1,data.shape[-1]))
        for j in range(len(data)):
            if c[j] == i:
                fenmu = fenmu +1
                x = x + data[j]
        u[i] = x / fenmu


def kmean(data,k=2):
    lenght = len(data)
    rr = np.random.randint(0,lenght-1,k)
    u = []
    for i in range(k):
        u.append(data[rr[i]])
    u = np.array(u)
    i = 0
    while i < 4:
        c = choice(data,u)
        updateU(data,c,u)
        i = i + 1
    return data,c,u



if __name__ == "__main__":
    mydata = creat_data_auto(2,200,2)
    print(mydata)
    print(mydata.shape)

    mydata = mydata.reshape(2*200, 2)

    medoids, centre, cluster = k_medoids(mydata, 2, euclidean_distance)
    
    data,c,u = kmean(mydata, 2)


    #print(medoids)

    #print(centre)

    #print(cluster)

    fig = plt.figure()
    #draw(mydata.reshape(2,200,2), fig)
    
    
    ax1 = fig.add_subplot(131)
    for i in range(mydata.reshape(2,200,2).shape[0]):
        ax1.scatter(mydata.reshape(2,200,2)[i,:,0],mydata.reshape(2,200,2)[i,:,1], c='g')

    #ax1.plot(medoids[:,0], medoids[:,1], color='r')
    #ax1.plot(u[:,0], u[:,1], color='k')
    plt.title('data')


    ax2 = fig.add_subplot(132)
    for i in range(len(centre)):
        ax2.scatter(mydata[np.nonzero(cluster[:]==i)][:,0],mydata[np.nonzero(cluster[:]==i)][:,1])

    ax2.plot(medoids[:,0], medoids[:,1], color='r')
    plt.title('k-medoids')
    
    ax3 = fig.add_subplot(133)
    for i in range(len(centre)):
        ax3.scatter(mydata[np.nonzero(c[:]==i)][:,0],mydata[np.nonzero(c[:]==i)][:,1])

    #ax3.plot(u[:,0], u[:,1], color='k')
    plt.title('k-means')
    plt.show()
   