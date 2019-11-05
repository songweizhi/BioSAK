#! /usr/bin/env python
# -*- coding: utf-8 -*-

import requests
from requests.auth import HTTPBasicAuth


class BacdiveClient(object):
#     This class provides method to get taxonomy informations corresponding
#     users entry(speciesname)
#
#     The attributes of the class are:
#     * ``headers`` -- sets the content-type of the HTTP request-header to json
#     * ``credentials`` -- attaches the username and the password to HTTPBasicAuth for using with the `requests` library

    headers = {'Accept': 'application/json'}
    USERNAME = 'songwz03@gmail.com'
    PASSWORD = '52shan88A'
    credentials = HTTPBasicAuth(USERNAME, PASSWORD)

#   _____________________________________________________________________________________________________________________________________
    def getAllLinks(self):
#         To get a list with links to the Detail View for BacDive ID, but only for the first page.
#         To consider the other pages, you habve to change the url to 'http://bacdive.dsmz.de/api/bacdive/bacdive_id/?page=2' etc.
#
#         INPUT:
        response = requests.get(
            'http://bacdive.dsmz.de/api/bacdive/bacdive_id/',
            headers=self.headers,
            auth=self.credentials
        )

        if response.status_code == 200:
            results = response.json()

#             OUTPUT:
#             object of type 'dict' with the fields 'count', 'previous', 'results', 'next'
#             the entries in field 'results' contain the URLs to correspondong strains separated by ',' e.g. {u'url:url1},{u'url:url2},{u'url:url3},etc
            return results

#   _____________________________________________________________________________________________________________________________________
    def getLinkByCultureno(self,culturecolnumber):
#         To get an URL to a strain corresponding the culture collection number
#
#         INPUT:
        response = requests.get(
            'http://bacdive.dsmz.de/api/bacdive/culturecollectionno/%s/' % (culturecolnumber),
            headers=self.headers,
            auth=self.credentials
        )

        if response.status_code == 200:
            results = response.json()

#             OUTPUT :
#             object of type 'dict' with the URL correspondong the culture collection number  e.g. {u'url:url1}
            return results

#   _____________________________________________________________________________________________________________________________________
    def getLinksByGenus(self, genus):
#         To get a list with links to the Detail View for BacDive ID, but only for the first page.
#         To consider the other pages, you habve to change the url to 'http://bacdive.dsmz.de/api/bacdive/taxon/{genus}/?page=2' etc.
#
#         INPUT:
        response = requests.get(
            'http://bacdive.dsmz.de/api/bacdive/taxon/%s/' % (genus ),
            headers=self.headers,
            auth=self.credentials
        )

        if response.status_code == 200:
            results = response.json()

#             OUTPUT:
#             object of type 'dict' with the fields 'count', 'previous', 'results', 'next'
#             the entries in field 'results' contain the URLs to correspondong strains separated by ',' e.g. {u'url:url1},{u'url:url2},{u'url:url3},etc
            return results

#   _____________________________________________________________________________________________________________________________________
    def getLinksBySpecies(self, genus, species_epithet):
#         To get a list with links to the Detail View for BacDive ID, but only for the first page.
#         To consider the other pages, you habve to change the url to 'http://bacdive.dsmz.de/api/bacdive/taxon/{genus}/{species_epithet}/?page=2' etc.
#
#         INPUT:
        response = requests.get(
            'http://bacdive.dsmz.de/api/bacdive/taxon/%s/%s/' % (genus,species_epithet),
            headers=self.headers,
            auth=self.credentials
        )

        if response.status_code == 200:
            results = response.json()

#             OUTPUT:
#             object of type 'dict' with the fields 'count', 'previous', 'results', 'next'
#             the entries in field 'results' contain the URLs to correspondong strains separated by ',' e.g. {u'url:url1},{u'url:url2},{u'url:url3},etc
            return results

#   _____________________________________________________________________________________________________________________________________
    def getLinksBySubspecies(self, genus, species_epithet, subspecies_epithet):
#         To get a list with links to the Detail View for BacDive ID, but only for the first page.
#         To consider the other pages, you habve to change the url to 'http://bacdive.dsmz.de/api/bacdive/taxon/{genus}/{species_epithet}/{subspecies_epithet}/?page=2' etc.
#
#         INPUT:
        response = requests.get(
            'http://bacdive.dsmz.de/api/bacdive/taxon/%s/%s/%s/' % (genus,species_epithet,subspecies_epithet),
            headers=self.headers,
            auth=self.credentials
        )

        if response.status_code == 200:
            results = response.json()

#             OUTPUT:
#             object of type 'dict' with the fields 'count', 'previous', 'results', 'next'
#             the entries in field 'results' contain the URLs to correspondong strains separated by ',' e.g. {u'url:url1},{u'url:url2},{u'url:url3},etc
            return results

#   _____________________________________________________________________________________________________________________________________
    def getLinksBySeqAccNum(self, seq_acc_num):
#         To get a list with links to the Detail View for BacDive ID. Usually it is a single entry and rarely few, mostly 2 or 3.
#
#         INPUT:
        response = requests.get('http://bacdive.dsmz.de/api/bacdive/sequence/%s/' % (seq_acc_num),
            headers=self.headers,
            auth=self.credentials
        )

        if response.status_code == 200:
            results = response.json()
#             OUTPUT:
#             object of type 'dict' with the fields 'count', 'previous', 'results', 'next'
#             the entries in field 'results' contain the URLs to correspondong strains separated by ',' e.g. {u'url:url1},{u'url:url2},{u'url:url3},etc
            return results
#   _____________________________________________________________________________________________________________________________________
    def run(self):
#         Run the client
        allLinks = self.getAllLinks()
        print(allLinks)
        cultureno = self.getLinkByCultureno('DSM 1')
        print(cultureno)
        genus = self.getLinksByGenus('acetobacter')
        print(genus)
        species = self.getLinksBySpecies('acetobacter','aceti')
        print(species)
        subspecies = self.getLinksBySubspecies('Bacillus','subtilis', 'subtilis')
        print(subspecies)
        sec_acc_num = self.getLinksBySeqAccNum('ALAS01000001')
        print(sec_acc_num)

if __name__ == '__main__':
    BacdiveClient().run()
