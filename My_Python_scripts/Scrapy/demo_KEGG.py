from requests import get
from bs4 import BeautifulSoup


web_address_K18831 = 'https://www.genome.jp/dbget-bin/www_bget?ko:K18831'
output = '/Users/songweizhi/Desktop/K18831.fasta'


web_address_K18831_soup = BeautifulSoup(get(web_address_K18831).text, 'html.parser')


#print(web_address_K18831_soup.find_all('a'))

# get total number of genes
total_num = 0
for gene in web_address_K18831_soup.find_all('a'):

    try:
        if '/dbget-bin/www_bget?' in gene['href']:
            total_num += 1
    except:
        pass


# download sequences
output_handle = open(output, 'w')
processed_seq = 1
for gene in web_address_K18831_soup.find_all('a'):

    try:
        if '/dbget-bin/www_bget?' in gene['href']:

            full_url = 'https://www.genome.jp/%s' % gene['href']
            full_url_soup = BeautifulSoup(get(full_url).text, 'html.parser')

            for each_row in full_url_soup.find_all('tr'):
                if each_row.text.find("AA seq") == 0:

                    print('Downloading the %s/%sth gene: %s' % (processed_seq, total_num, gene['href'].split('/dbget-bin/www_bget?')[1]))
                    output_handle.write('>%s %s\n' % (gene['href'].split('/dbget-bin/www_bget?')[1], each_row.td.text))

                    processed_seq += 1

    except:
        pass

output_handle.close()
