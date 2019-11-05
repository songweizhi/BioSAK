import requests
from bs4 import BeautifulSoup


# https://www.dataquest.io/blog/web-scraping-tutorial-python/






page = requests.get("http://dataquestio.github.io/web-scraping-pages/simple.html")
print(page)
soup = BeautifulSoup(page.content, 'html.parser')

html = list(soup.children)[2]


soup = BeautifulSoup(page.content, 'html.parser')
soup.find_all('p')