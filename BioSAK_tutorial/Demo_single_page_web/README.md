
Script for generating a single-page website to visualize biodiversity data
---


### Produce a webpage without bar chart
  
    python3 get_single_page_web.py -m doc/metadata.txt -w doc/template.html -o single_page_web.html

![webpage_without_barchart.png](doc/webpage_without_barchart.png)


### Produce a webpage with bar chart

    python3 get_single_page_web.py -m doc/metadata.txt -w doc/template_with_barchart.html -o single_page_web_with_barchart.html

![webpage_with_barchart.png](doc/webpage_with_barchart.png)
