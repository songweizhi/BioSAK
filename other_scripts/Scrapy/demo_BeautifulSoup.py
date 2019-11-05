from bs4 import BeautifulSoup
from urllib.request import urlopen


html = urlopen("https://morvanzhou.github.io/static/scraping/basic-structure.html").read().decode('utf-8')
html_soup = BeautifulSoup(html, features='lxml')


print('\nhtml_soup')
print(html_soup)


print('\n\nhtml_soup.h1')
print(html_soup.h1)


print('\n\nhtml_soup.title')
print(html_soup.title)

print("\n\nhtml_soup.find_all('a')")
print(html_soup.find_all('a'))


for each in html_soup.find_all('a'):
    print()
    print(each)
    print(each['href'])

