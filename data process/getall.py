from bs4 import BeautifulSoup
import requests
import sys
import time
import random

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

sys.setrecursionlimit(10000)

acc_num = 0

count = 0

starttime = time.clock()

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
				101.236.19.165:8866
				118.190.95.35:9001
				60.216.177.152:8118
				101.236.22.141:8866
				118.31.220.3:8080
				118.190.95.43:9001
				118.212.137.135:31288
				114.215.95.188:3128
				139.224.80.139:3128
				218.60.8.98:3129
				119.188.162.165:8081
				121.42.167.160:3128
				123.138.89.132:9999
				110.185.227.237:9999
			"""
headers = {}

proxies = {}

# link = [] #判断蛋白质文件是否重复



def save(filename, tree):
	tree = ET.ElementTree(tree)
	tree.write('test.xml', encoding="us-ascii", xml_declaration=True, default_namespace=None, method="xml", short_empty_elements=False)


def travel_get(url, depth, parent):

	print(url, ": depth=", depth)

	global count

	global acc_num

	acc_num = acc_num + 1

	#if acc_num % 4 == 0:
	time.sleep(random.uniform(1, 4))

	if count % 6 == 0:
		#random.choice(USER_AGENTS)
		headers["user-agent"] = random.choice(USER_AGENTS)
		proxies["http"] = random.choice(IP_proxies)
	# if count % 60 == 0:
	# 	time.sleep(random.uniform(1,2)) 
	# if count % 666 == 0:
	# 	time.sleep(random.uniform(5, 10))

	if (count + 1) % 3000 == 0:
		time.sleep(random.uniform(1000,2400))

	try:
		res = requests.get(url, headers=headers, proxies=proxies)
	except:
		print("headers:=>\n", headers, "\nproxies:=>\n", proxies)
		res = requests.get(url)

	soup = BeautifulSoup(res.text, "html.parser")

	mylist = soup.find_all(class_="browse")

	name =  soup.find_all("h3")[1].contents[0][:-1]

	if len(mylist) < 2:
		print("error! url is ", url, " user-agent = ", headers["user-agent"])
		return

	if depth == 4 or depth == 5:	#protein #domian
		for i in mylist[1].find_all(class_="browse"):	#移除掉下下层的节点
			#print(i.prettify())
			i.clear()	

	mylinks = mylist[1].find_all("a", "sunid")
	attr = {}

	#global link
	for a in mylinks:
		#time.sleep(random.uniform(1, 4))

		if a.img != None:
			continue
		attr["data-sunid"] = a["data-sunid"]

		# if a["href"] in link:
		# 	continue
		# link.append(a["href"])
		attr["href"] = a["href"]
		if depth == 6:

			############# debug info ############# 
			
			if count % 100 == 0:
				print(url, ": depth =", depth, "  count =", count, "  time =", int(time.clock()-starttime), "s")
			############# debug info ############# 
			time.sleep(random.uniform(1, 4))

			attr["name"] = a.contents[0].split(':')[0]

			mypdb = requests.get(a["href"], headers=headers, proxies=proxies)

			pdbsoup = BeautifulSoup(mypdb.text, "html.parser")

			mydiv = pdbsoup.find_all("div", "indented")

			if len(mydiv) < 2:
				print("pdb download error! url is ", a["href"], " user-agent = ", headers["user-agent"])
			else:
				attr["download"] = mydiv[1].find("a")["href"]	#获取pdb文件的下载地址
				node = ET.SubElement(parent, "PDB", attr)
				count = count + 1

		else:
			attr["name"] = a.contents[0]

			node = ET.SubElement(parent, name, attr)
			
			travel_get(a["href"], depth + 1, node)		#递归创建


if '__main__' == __name__:
	
	#global acc_num

	requests.adapters.DEFAULT_RETRIES = 5

	s = requests.session()
	s.keep_alive = False

	USER_AGENTS = USER_AGENTS.split()

	IP_proxies = IP_proxies.split()

	# print(random.choice(IP_proxies))

	# sys.exit()

	# rooturl = "https://scop.berkeley.edu/sunid=0"

	# res = requests.get(rooturl)

	# soup = BeautifulSoup(res.text, "html.parser")

	# myclass = soup.find_all(class_="browse")

	# mylinks = myclass[0].find_all("a", "sunid")

	# for i in mylinks:
	# 	print(i.prettify())

	root = ET.Element('root')



	# # protein depth=5
	# sub = ET.SubElement(root, 'Protein', {"sunid":"46456","href":"https://scop.berkeley.edu/sunid=46460", "name":"Protozoan/bacterial hemoglobin"})

	# travel_get("https://scop.berkeley.edu/sunid=46460", 5, sub)


	# class depth=1
	sub = ET.SubElement(root, 'Protein', {"sunid":"46456","href":"https://scop.berkeley.edu/sunid=46456", "name":"All alpha proteins"})

	travel_get("https://scop.berkeley.edu/sunid=46456", 1, sub)
	

	tree = ET.ElementTree(root)
	tree.write('test.xml', encoding="us-ascii", xml_declaration=True, default_namespace=None, method="xml", short_empty_elements=False)
	print("\nacc_num=", acc_num)	

