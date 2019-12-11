import re
import random
from bs4 import BeautifulSoup
from urllib.request import urlopen


base_url = 'https://baike.baidu.com/'
url = 'https://baike.baidu.com/item/%E7%BD%91%E7%BB%9C%E7%88%AC%E8%99%AB/5162711'


html = urlopen(url).read().decode('utf-8')
html_soup = BeautifulSoup(html, features='lxml')

print(html_soup.h1.get_text())
print(html_soup.find('h1').get_text())
