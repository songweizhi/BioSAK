from requests import get
from bs4 import BeautifulSoup


# https://www.dataquest.io/blog/web-scraping-beautifulsoup/

url = 'http://www.imdb.com/search/title?release_date=2017&sort=num_votes,desc&page=1'
html_soup = BeautifulSoup(get(url).text, 'html.parser')
print(html_soup)

movie_containers = html_soup.find_all('div', class_='lister-item mode-advanced')


# Extract data
print('Year\tIMDB\tM_Score\tVote\tName')
for movie_container in movie_containers:

    # If the movie has Metascore, then extract:
    if movie_container.find('div', class_='ratings-metascore') is not None:
        name =    movie_container.h3.a.text
        year =    movie_container.h3.find('span', class_='lister-item-year').text
        imdb =    movie_container.strong.text
        m_score = movie_container.find('span', class_='metascore').text
        vote =    movie_container.find('span', attrs={'name': 'nv'})['data-value']

        print('%s\t%s\t%s\t%s\t%s' % (year, imdb, m_score, vote, name))

