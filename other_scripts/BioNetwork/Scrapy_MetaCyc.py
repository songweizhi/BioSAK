from requests import get
from bs4 import BeautifulSoup

MetaCyc_url = 'https://metacyc.org/META/class-instances?object=Reactions'
MetaCyc_url = 'https://metacyc.org/META/search-query?type=REACTION&name=1.1.1.-'

MetaCyc_soup = BeautifulSoup(get(MetaCyc_url).text, 'html.parser')

reaction_containers = MetaCyc_soup.find_all('tr', class_={'trDark', 'trLight'})

for reaction_container in reaction_containers:

    reaction = reaction_container.td.a.text
    reaction_id = reaction_container.td.a['href'].split('object=')[-1]

    print('%s\t%s' % (reaction_id, reaction))
