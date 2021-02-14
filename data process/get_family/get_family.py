
"""
https://blog.csdn.net/wolinghuanyun/article/details/52474515

<tag attrib = > text </tag> tail 
例：<APP_KEY channel = 'CSDN'> hello123456789 </APP_KEY> 
- tag，即标签，用于标识该元素表示哪种数据，即APP_KEY 
- attrib，即属性，用Dictionary形式保存，即{‘channel’ = ‘CSDN’} 
- text，文本字符串，可以用来存储一些数据，即hello123456789 
- tail，尾字符串，并不是必须的，例子中没有包含。
"""


#coding:utf-8
from lxml import etree
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

import sys
import requests
from sklearn.externals import joblib


sys.setrecursionlimit(10000)

count = 0
attrib = {}
attr_list = []

lavel = ["Folds", "Superfamilies", "Families", "Protein", "Species", "domain", "subdomain", "PDB"]
mylav = {"Folds":1, "Superfamilies":2, "Families":3, "Protein":4, "Species":5, "domain":6, "subdomain":7, "PDB":8}

myre = []

flag = False

pre_family = ''

def dfs(myroot, depth, myattrib, lav):
	
	global count, myre, flag

	if depth == 5:
		myre.append(myroot.attrib["name"])
	if depth == lav:
		
		""""""
		try:
			myattrib["PDB_name"] = myroot.attrib["name"]
			myattrib["download"] = myroot.attrib["download"]

			joblib.dump(myattrib, myroot.attrib["name"].split()[0]+".info")

			with open(myroot.attrib["name"].split()[0]+".ent", "w") as f:
				text =  requests.get(myattrib["download"])
				f.write(text.text)
			with open("result.txt", "a") as f:
				strs = ""
				for i in myattrib.keys():
					
					strs +=  i + " : " + myattrib[i] + "\n"
				f.write(strs+"\n\n\n")
		except:
			print("shit, download error!", "\n", myroot.attrib)
			return

		if not flag:
			count = count + 1
			flag = True
		
		return
	
	for (i, t) in enumerate(myroot.getchildren()):
		#print(t.tag, ":", t.attrib)
		myattrib[lavel[depth]] = ''
		try:
			myattrib[lavel[depth]] = t.attrib["name"]
		except:
			pass


		if depth == 2:
			flag = False
			pre_family = t.attrib["name"]

		dfs(t, depth+1, myattrib, lav) 



		# if not flag:
		# 	dfs(t, depth+1, myattrib, lav) 

		# else:
		# 	if depth == 3:
		# 		print("bb")
		# 		flag = False
		# 	else:
		# 		return



if '__main__' == __name__:

	alpha = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"]
	
	for i in alpha:

		tree = ET.ElementTree(file=str(i)+'.xml')

		root = tree.getroot()

		dfs(root.getchildren()[0], 0, attrib, 8)
		# for (i, child) in enumerate(tree.iter(tag='PDB')):
		# 	count = count + 1
		# 	print(child.attrib, child.find("").attrib)
		# 	sys.exit()
			
		print("family",":",count)
	
	# for i in root.find("./Classes/Folds").iter():
	# 	print(i.attrib)