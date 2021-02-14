# -*- coding: utf-8 -*-
from bs4 import BeautifulSoup
import requests
import sys
import time
import random
import datetime


# import httplib
# httplib.HTTPConnection._http_vsn = 10
# httplib.HTTPConnection._http_vsn_str = 'HTTP/1.0'

import http.client
http.client.HTTPConnection._http_vsn = 10
http.client.HTTPConnection._http_vsn_str = 'HTTP/1.0'

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

sys.setrecursionlimit(10000)

acc_num = 10

count = 0

#starttime = time.clock()

USER_AGENTS = """
				Mozilla/5.0(Windows;U;WindowsNT6.1;en-us)AppleWebKit/534.50(KHTML,likeGecko)Version/5.1Safari/534.50

				Mozilla/5.0(Macintosh;U;IntelMacOSX10_6_8;en-us)AppleWebKit/534.50(KHTML,likeGecko)Version/5.1Safari/534.50

				Mozilla/5.0(Macintosh;IntelMacOSX10.6;rv:2.0.1)Gecko/20100101Firefox/4.0.1

				Mozilla/5.0(WindowsNT6.1;rv:2.0.1)Gecko/20100101Firefox/4.0.1

				Opera/9.80(Macintosh;IntelMacOSX10.6.8;U;en)Presto/2.8.131Version/11.11

				Opera/9.80(WindowsNT6.1;U;en)Presto/2.8.131Version/11.11

				Mozilla/5.0(Macintosh;IntelMacOSX10_7_0)AppleWebKit/535.11(KHTML,likeGecko)Chrome/17.0.963.56Safari/535.11

				Mozilla/5.0(iPhone;U;CPUiPhoneOS4_3_3likeMacOSX;en-us)AppleWebKit/533.17.9(KHTML,likeGecko)Version/5.0.2Mobile/8J2Safari/6533.18.5

				Mozilla/5.0(iPod;U;CPUiPhoneOS4_3_3likeMacOSX;en-us)AppleWebKit/533.17.9(KHTML,likeGecko)Version/5.0.2Mobile/8J2Safari/6533.18.5

				Mozilla/5.0(iPad;U;CPUOS4_3_3likeMacOSX;en-us)AppleWebKit/533.17.9(KHTML,likeGecko)Version/5.0.2Mobile/8J2Safari/6533.18.5

				Mozilla/5.0(Linux;U;Android2.3.7;en-us;NexusOneBuild/FRF91)AppleWebKit/533.1(KHTML,likeGecko)Version/4.0MobileSafari/533.1

				MQQBrowser/26Mozilla/5.0(Linux;U;Android2.3.7;zh-cn;MB200Build/GRJ22;CyanogenMod-7)AppleWebKit/533.1(KHTML,likeGecko)Version/4.0MobileSafari/533.1

				Opera/9.80(Android2.3.4;Linux;OperaMobi/build-1107180945;U;en-GB)Presto/2.8.149Version/11.10

				Mozilla/5.0(Linux;U;Android3.0;en-us;XoomBuild/HRI39)AppleWebKit/534.13(KHTML,likeGecko)Version/4.0Safari/534.13
			"""

IP_proxies = """
				101.236.22.141:8866
				114.215.95.188:3128
				139.224.80.139:3128
				218.60.8.98:3129
				121.42.167.160:3128

				49.70.209.159:9000
				139.196.51.201:8118
				110.72.242.203:53281

				218.66.253.145:8800
				101.93.200.150:9000
				101.37.79.125:3128
				124.207.82.166:8008
				114.113.126.83:80
				120.83.48.153:9000
				218.66.253.145:8800
				121.13.54.251:808
				203.130.46.108:9090
				171.88.52.125:9999
				123.172.68.67:53281
				14.20.235.87:9797
				120.34.73.236:53281
				39.108.76.176:3128
				119.57.108.65:53281
			 """

headers = {"Accept-Encoding":""}

proxies = {}

lavel = ["Folds", "Superfamilies", "Families", "Protein", "Species", "domain", "subdomain", "PDB"]

# link = [] #判断蛋白质文件是否重复



def save(filename, tree):
	tree = ET.ElementTree(tree)
	tree.write('test.xml', encoding="us-ascii", xml_declaration=True, default_namespace=None, method="xml", short_empty_elements=False)


def travel_get(url, depth, parent):

	print(url, ": depth=", depth)

	global count, acc_num, headers, proxies, lavel, USER_AGENTS, IP_proxies

	acc_num = acc_num + 1

	while True:
		try:
			res = requests.get(url, headers=headers, proxies=proxies)
			break
		except requests.exceptions.ConnectionError:
			print("connect refused!\nurl = ", url)
			print("headers:\n", headers, "\nproxies:\n", proxies)
			headers["user-agent"] = random.choice(USER_AGENTS)
			#proxies["http"] = random.choice(IP_proxies)
			time.sleep(5)
			continue

	soup = BeautifulSoup(res.text, "lxml")

	mylist = soup.find_all(class_="browse")

	ispagination = soup.find_all(class_="pagination")

	if len(ispagination) > 0:	#ispagination
		print("====="+url+"=====")
		sys.exit()
		with open("pagination_e.txt", "a") as f:
			f.white(str(url)+"\n")
	#name =  soup.find_all("h3")[1].contents[0][:-1]

	name = lavel[depth-1]

	if len(mylist) < 2:
		print("******error!******\n url = ", url, " \nuser-agent = ", headers["user-agent"])
		return

	if depth == 4 or depth == 5:	#protein #domian
		for i in mylist[1].find_all(class_="browse"):	#移除掉下下层的节点
			#print(i.prettify())
			i.clear()	

	
	attr = {}

	if depth == 6:
		for li in mylist[1].children:
			
			# if li.a.text == "1dlw":
			# 	print(url)
			# 	for g, xx in enumerate(mylist[1].children):
			# 		print(g)
			# 		print(xx)
			# 	sys.exit()

			domain_attr = {}
			domain_attr["name"] = li.a.text
			domain_attr["href"] = li.a["href"]
			
			domain = ET.SubElement(parent, "domain", domain_attr)
			
			lis = li.ul.children if li.ol == None else li.ol.children
			for	subli in lis: 
				if type(subli) != type(mylist[1]):
					continue
				
				subdomain_attr = {}
				try:
					if subli.table.a["title"] != None or str(subli.table.a["title"]) != "":
						subdomain_attr["name"] = subli.table.a["title"].split(':')[1]
						subdomain_attr["scope"] = subli.table.a["title"].split(':')[2].split()[0]

					aa = subli.table.find_all("a", "sunid")[1]

					subdomain_attr["name"] = aa.contents[0].split(':')[1]

					#subdomain_attr["scope"] = subli.table.a["title"].split(':')[2].split()[0]
					#print(subli.table.a.contents[0])
					scopess = aa.contents[0].split(':')[2].split()

					if len(scopess) > 0:
						subdomain_attr["scope"] = scopess[0]
					else:
						subdomain_attr["scope"] = "None"
				except:
					print("subli.table.find_all(\"a\", \"sunid\") except \nurl:", url)
				
				subdomain = ET.SubElement(domain, "subdomain", subdomain_attr)
				mylinks = subli.find_all("a", "sunid")
				for a in mylinks:
					if a.img != None:
						continue
					attr = {}
					attr["data-sunid"] = a["data-sunid"]

					############# debug info ############# 
					
					if count % 100 == 0:
						print(url, ": depth =", depth, "  count =", count, "  time =", int(time.time()-startTime), "s")
					############# debug info ############# 
					#time.sleep(random.uniform(1, 2))

					attr["name"] = a.contents[0].split(':')[0]


					while True:
						try:
							#res = requests.get(url, headers=headers, proxies=proxies)
							mypdb = requests.get(a["href"], headers=headers, proxies=proxies)
							break
						except requests.exceptions.ConnectionError:
							print("connect refused!\n", url)
							print("headers:=>\n", headers, "\nproxies:=>\n", proxies)
							headers["user-agent"] = random.choice(USER_AGENTS)
							#proxies["http"] = random.choice(IP_proxies)
							time.sleep(5)
							continue

					
					pdbsoup = BeautifulSoup(mypdb.text, "lxml")

					mydiv = pdbsoup.find_all("div", "indented")

					if len(mydiv) < 2:
						print("pdb download error! url is ", a["href"], " user-agent = ", headers["user-agent"])
					else:
						attr["download"] = mydiv[1].find("a")["href"]	#获取pdb文件的下载地址
						node = ET.SubElement(subdomain, "PDB", attr)
						count = count + 1


	########################################################

	else:
		mylinks = mylist[1].find_all("a", "sunid")
		for a in mylinks:

			if a.img != None:
				continue
			
			# if a.contents[0] == "Ciliate (Paramecium caudatum)":
			# 	print(url)
			# 	for x in mylinks:
			# 		print(x)
			# 	sys.exit()
			attr["data-sunid"] = a["data-sunid"]

			# if a["href"] in link:
			# 	continue
			# link.append(a["href"])
			attr["href"] = a["href"]

			attr["name"] = a.contents[0]

			node = ET.SubElement(parent, name, attr)
			
			travel_get(a["href"], depth + 1, node)		#递归创建


if '__main__' == __name__:
	
	# global acc_num, count

	requests.adapters.DEFAULT_RETRIES = 5

	s = requests.session()
	s.keep_alive = False

	USER_AGENTS = USER_AGENTS.split()

	IP_proxies = IP_proxies.split()

	headers["user-agent"] = random.choice(USER_AGENTS)
	
	#proxies["http"] = random.choice(IP_proxies)


	root = ET.Element('root')


	# class depth=1
	sub = ET.SubElement(root, 'Classes', {"sunid":"56572","href":"http://scop.berkeley.edu/sunid=56572", "name":"Multi-domain proteins (alpha and beta)"})

	#startTime = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

	startTime = time.time()

	travel_get("http://scop.berkeley.edu/sunid=56572", 1, sub)
	
	# travel_get("http://scop.berkeley.edu/sunid=88552", 6, sub)

	tree = ET.ElementTree(root)
	tree.write('e.xml', encoding="us-ascii", xml_declaration=True, default_namespace=None, method="xml", short_empty_elements=False)
	
	with open("results_e.txt", 'w') as f:
		strs = "acc_num=" + str(acc_num) + "\n" + "pdb num = " + str(count)
		f.write(strs)

