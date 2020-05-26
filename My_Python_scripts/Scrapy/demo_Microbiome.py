import requests
import json
from lxml import etree
import time
import datetime


#获取网页相应内容
def get_one_page(url):

    headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/65.0.3325.181 Safari/537.36'}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.text
    else:
        return None



n = 1
out_html = '/Users/songweizhi/Desktop/out_html.txt'

# 查找关键字网址
url = 'https://microbiomejournal.biomedcentral.com/articles?searchType=journalSearch&sort=PubDate&page={num}'.format(num=n)
print(url)

# 获取相应内容
headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/65.0.3325.181 Safari/537.36'}
response = requests.get(url, headers=headers)
html = response.text
print(html)
out_html_handle = open(out_html, 'w')
out_html_handle.write(html)
out_html_handle.close()

# xpath解析文献href
html = etree.HTML(html)
href = html.xpath('//h3/a[@data-test="title-link"]/@href')
print(href)




#
# # xpath解析文献href
# html = etree.HTML(html)
# href = html.xpath('//h3/a[@data-test="title-link"]/@href')
#
# # 对每个href进行关键字筛选
# for i in href:
#         # 单篇文献url
#         url = 'https://microbiomejournal.biomedcentral.com' + str(i)
#         doi = i.replace('/articles', '')
#         url_pdf = 'https://microbiomejournal.biomedcentral.com/track/pdf' + doi + '#page=1'
#         print(url_pdf)
#         html = get_one_page(url)
#         name = []
#         # 判断关键字是否存在于文献中
#         for item in parse_one_page(html):
#                 with open('data.csv', 'a', encoding='utf-8') as f:
#                         f.write(json.dumps(item['Title'], ensure_ascii=False) + ',')
#                         f.write(json.dumps(item['DOI'], ensure_ascii=False) + ',')
#                         f.write(json.dumps(item['Accesses'], ensure_ascii=False) + ',')
#                         f.write(json.dumps(item['Citations'], ensure_ascii=False) + ',')
#                         f.write(json.dumps(item['Altmetric'], ensure_ascii=False) + ',')
#                         f.write(json.dumps(item['published'], ensure_ascii=False) + ',')
#                         f.write(json.dumps(item['author'], ensure_ascii=False) + '\n')
#                         time.sleep(2)
#                         f.close()
#                 print('======读取pdf中=====')
#                 r = get_one_page_pdf(url_pdf)
#                 print('======保存pdf中=====')
#                 with open('data_pdf/' + item['Title'][0] + '.pdf', 'wb') as p:
#                         p.write(r.content)
#                         print('======保存完毕=====')
#                         p.close()
#                 time.sleep(2)
#

