
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

sys.setrecursionlimit(10000)




#print(type(tree.iter()))
#getchildren()
#root.iter()

count = 0
attrib = {}
lavel = ["Folds", "Superfamilies", "Families", "Protein", "Species", "domain", "subdomain", "PDB"]
mylav = {"Folds":1, "Superfamilies":2, "Families":3, "Protein":4, "Species":5, "domain":6, "subdomain":7, "PDB":8}

myre = []

def dfs(myroot, depth, myattrib, lav):
	
	global count, myre

	if depth == 5:
		myre.append(myroot.attrib["name"])
	if depth == lav:
		# for (i, t) in enumerate(myroot.getchildren()):
		# 	print(t.tag, ":", t.attrib)
		
		count = count + 1

		"""
		if count % 1710 == 0:
			
			#try:
			myattrib["PDB_name"] = myroot.attrib["name"]
			myattrib["download"] = myroot.attrib["download"]

			with open(myroot.attrib["name"].split()[0]+".ent", "w") as f:
				text =  requests.get(myattrib["download"])
				f.write(text.text)
			with open("result.txt", "a") as f:
				strs = ""
				for i in myattrib.keys():
					
					strs +=  i + " : " + myattrib[i] + "\n"
				f.write(strs+"\n\n\n")
			# except:
			# 	print("shit, error!")
		"""
		#print(myattrib)

		return 
	
	for (i, t) in enumerate(myroot.getchildren()):
		#print(t.tag, ":", t.attrib)
		myattrib[lavel[depth]] = ''
		try:
			myattrib[lavel[depth]] = t.attrib["name"]
		except:
			pass
			#print("fail")
		# if count%10000==0:
		# 	print(myattrib)
			#sys.exit()
		
		if lav >= 5:
			if depth == 4:
				if str(t.attrib["name"]) == "Human (Homo sapiens)":		# 只要人类的
					dfs(t, depth+1, myattrib, lav) 
			else:
				dfs(t, depth+1, myattrib, lav) 

		else:
			dfs(t, depth+1, myattrib, lav)

		'''	
		if lav >= 5:
			if depth == 4:
				if str(t.attrib["name"]) != "Human (Homo sapiens)":		# 只要人类的
					return
				else:
					if str(t.attrib["name"]) == "Human (Homo sapiens)":
						dfs(t, depth+1, myattrib, lav) 
			else:
				dfs(t, depth+1, myattrib, lav) 

		else:
			dfs(t, depth+1, myattrib, lav)
		'''



		# if lav >= 5:
		# 	if depth == 4:
		# 		if t.attrib["name"] != "Human (Homo sapiens)":		# 只要人类的
		# 			return


		# dfs(t, depth+1, myattrib, lav)





# for (i, t) in enumerate(root.getchildren()):
# 	print(t.tag, ":", t.attrib)
	
# 	for (j, x) in enumerate(t.getchildren()):
# 		print(x.tag, ":", x.attrib)
# 		if j==10:
# 			sys.exit()


if '__main__' == __name__:

	alpha = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"]
	
	for l in lavel:
		count = 0
		for i in alpha:

			tree = ET.ElementTree(file=str(i)+'.xml')

			root = tree.getroot()

			dfs(root.getchildren()[0], 0, attrib, mylav[l])
			# for (i, child) in enumerate(tree.iter(tag='PDB')):
			# 	count = count + 1
			# 	print(child.attrib, child.find("").attrib)
			# 	sys.exit()
			
		print(l,":",count)
		
	
	print(set(myre))

	# for i in root.find("./Classes/Folds").iter():
	# 	print(i.attrib)