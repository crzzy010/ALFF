
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


flag = 0

def dfs(myroot, depth, myattrib, lav):
	
	global count, flag


	if depth == 6:
		if (flag+1) % 100 == 0:
			print("finish")
			print(count)
			sys.exit()
		else:
			flag = flag + 1
	if depth == 8:
		# for (i, t) in enumerate(myroot.getchildren()):
		# 	print(t.tag, ":", t.attrib)
		
		

		"""
		if count % 1710 == 0:
		"""	
		if flag == 88:
			try:
				count = count + 1
				myattrib["PDB_name"] = myroot.attrib["name"]
				myattrib["download"] = myroot.attrib["download"]

				with open(myroot.attrib["name"].split()[0]+".ent", "w") as f:
					text =  requests.get(myattrib["download"])
					f.write(text.text)
				with open("dddr.txt", "a") as f:
					strs = ""
					for i in myattrib.keys():
						
						strs +=  i + " : " + myattrib[i] + "\n"
					f.write(strs+"\n\n\n")
			except:
				print("shit, error!")
		
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
		dfs(t, depth+1, myattrib, lav)





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
		
	


	# for i in root.find("./Classes/Folds").iter():
	# 	print(i.attrib)