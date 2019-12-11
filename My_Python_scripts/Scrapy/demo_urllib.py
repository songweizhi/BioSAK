import re
from urllib.request import urlopen


# if has Chinese, apply decode()
html = urlopen("https://morvanzhou.github.io/static/scraping/basic-structure.html").read().decode('utf-8')


# get Page title
Page_title = re.findall(r"<title>(.+?)</title>", html)[0]
print('\n------Page_title------')
print(Page_title)


# get Page paragraph
Page_paragraph = re.findall(r"<p>(.*?)</p>", html, flags=re.DOTALL)[0]   # re.DOTALL if multi line
print('\n------Page_paragraph------')
print(Page_paragraph)


# select links
Page_links = re.findall(r'href="(.*?)"', html)
print('\n------Page_links------')
print(Page_links)

