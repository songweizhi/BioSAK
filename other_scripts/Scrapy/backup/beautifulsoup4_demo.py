import urllib.request
req = urllib.request.Request(url='https://localhost/cgi-bin/test.cgi',
                      data=b'This data is passed to stdin of the CGI')
with urllib.request.urlopen(req) as f:
    print(f.read().decode('utf-8'))
    