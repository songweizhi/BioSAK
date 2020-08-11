from bs4 import BeautifulSoup


# html_doc = """
# <html><head><title>The Dormouse's story</title></head>
# <body>
# <p class="title"><b>The Dormouse's story</b></p>
#
# <p class="story">Once upon a time there were three little sisters; and their names were
# <a href="http://example.com/elsie" class="sister" id="link1">Elsie</a>,
# <a href="http://example.com/lacie" class="sister" id="link2">Lacie</a> and
# <a href="http://example.com/tillie" class="sister" id="link3">Tillie</a>;
# and they lived at the bottom of a well.</p>
#
# <p class="story">...</p>
# """
#
# soup = BeautifulSoup(html_doc, 'html.parser')
#
# print(soup.title)
# print(soup.title.string)
# print(soup.title.text)
# print(soup.p)
# print(soup.p['class'])
# print(soup.a)
# print(soup.find_all('a'))
# print(soup.find(id="link1"))
# print(soup.find(id="link2"))
# print(soup.find(id="link3"))
# print(soup.get_text())
#
#
# soup = BeautifulSoup('<b class="boldest" id="DemoID" > Extremely bold </b>', features='lxml')
# print(soup.b)
# print(soup.b.name)
# print(soup.b['class'])
# print(soup.b.get('class'))
# print(soup.b['id'])
#
#
# html_doc = """
# <html><head><title>The Dormouse's story</title></head>
# <body>
# <p class="title"><b>The Dormouse's story</b></p>
#
# <p class="story">Once upon a time there were three little sisters; and their names were
# <a href="http://example.com/elsie" class="sister" id="link1">Elsie</a>,
# <a href="http://example.com/lacie" class="sister" id="link2">Lacie</a> and
# <a href="http://example.com/tillie" class="sister" id="link3">Tillie</a>;
# and they lived at the bottom of a well.</p>
#
# <p class="story">...</p>
# """
#
# from bs4 import BeautifulSoup
# soup = BeautifulSoup(html_doc, 'html.parser')
#
# print(soup.head)
# print(soup.title)
# print(soup.title.text)
# print(soup.body)
# print(soup.body.b)
# print(soup.a)
# print()




test_html = '''
<tr><th class="th40" align="left" valign="top" style="border-color:#000; border-width: 1px 0px 0px 1px; border-style: solid"><nobr>Other DBs</nobr></th>
<td class="td40" style="border-color:#000; border-width: 1px 1px 0px 1px; border-style: solid"><table style="border:0;border-collapse:collapse;"><tr><td valign="top"><nobr>COG:&nbsp;</nobr></td><td><a href="https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=COG5499">COG5499</a></td></tr></table></td></tr>
'''

soup = BeautifulSoup(test_html, features="lxml")


print()
print(soup.th.text)
print(soup.td.text)

print(soup.find_all('th'))

print(soup.th)



