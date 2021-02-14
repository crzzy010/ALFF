#-*- coding:utf-8 -*-

import sys
import os

mulu = "C:\\Users\\Crzzy\\Desktop\\data\\qwer"

#sys.setdefaultencoding('utf-8')

def check_chinese(str):
    return u'\u4e00' <= str <= u'\u9fa5'


def getdata(path, name, fre):
    print("The  "+path+name+"\\"+"  has start:"+"\n\n\n------")
    with open(path+"\\"+name, "r") as f:
        while True:
            line = f.readline()  #逐行读文本
            if not line:
                break
            else:
                strline = str(line).strip().split()
                if len(strline) == 12:
                    if strline[2] == "CA" or strline[2] == "CB":    #只读CA CB的数据
                        fre.write(name+"    "+strline[2]+"    "+strline[6]+"  "+strline[7]+"  "+strline[8]+"\n")
    return


def travel(path):
    if os.path.exists(sys.path[0]+"\\result.txt"):
        os.remove(sys.path[0]+"\\result.txt")
    #os.mknod("result.txt")
    queue = []
    queue.append(path)

    with open(path+"\\result.txt", "a") as fre:

        while len(queue) > 0:
            path = queue.pop(0)
            files = os.listdir(path)
            for f in files:
                ff = os.path.join(path, f)
                if os.path.isdir(ff):
                    queue.append(ff)
                   # print(f+"--------")
                else:
                    #print(ff)
                    if not f == 'getdata.py': 
                        getdata(path, f, fre)


if __name__ == "__main__":
   travel(mulu)
