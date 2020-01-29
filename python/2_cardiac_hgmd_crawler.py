"""
modified from <http://stockrt.github.com/p/emulating-a-browser-in-python-with-mechanize/>
"""
import mechanize
import http.cookiejar
import html2text
import re
import time
import sys
from bs4 import BeautifulSoup


hgmd_login_url = "http://www.hgmd.cf.ac.uk/docs/login.html"
email_address = "ssubedi@houstonmethodist.org"
password = "HGMD881656"
if email_address == "":
    sys.exit("define your email address and password")



def initialize_browser():

    br = mechanize.Browser()
    # Cookie Jar
    cj = http.cookiejar.LWPCookieJar()
    br.set_cookiejar(cj)

    # Browser options
    br.set_handle_equiv(True)
#    br.set_handle_gzip(True)
    br.set_handle_redirect(True)
    br.set_handle_referer(True)
    br.set_handle_robots(False)

    # Follows refresh 0 but not hangs on refresh > 0
    br.set_handle_refresh(mechanize._http.HTTPRefreshProcessor(), max_time=1)

    # Want debugging messages?
#    br.set_debug_http(True)
#    br.set_debug_redirects(True)
#    br.set_debug_responses(True)

    # User-Agent (this is cheating, ok?)
    br.addheaders = [('User-agent', 'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.9.0.1) Gecko/2008071615 Fedora/3.0.1-1.fc9 Firefox/3.0.1')]
    br.addheaders.append(('email', email_address))

    return br


def login_hgmd(br):
    """
    login to HGMD

    After calling this function, you will be able to search in HGMD programmatically
    """
    response = br.open(hgmd_login_url)
    html = response.read()

    # print ( response to STDOUT for debugging purposes
    # the html2text library is used for formatting the output in a more readable form
    # print (html)

    # print all the forms in the current page
    # print ([f for f in br.forms()])

    # select login form
    br.select_form(nr=0)
    # print  (br.form)

    # print ( ( all controls in the current form, for debugging purposes
    # print ([c.name for c in br.form.controls])

    # set username and password
    br.form['email'] = email_address
    br.form['password'] = password

    # submit form
    response_form = br.submit()

    # Now, you should have successfully logged in. The contents of the page will be changed. Check the contents of br.read()
    # html_response = response_form.read()
    # print (html_response)

    # wait 2 seconds to not overload the server
    time.sleep(2)

    return br

def browse_get(br, gene):
    query =  br.open("http://www.hgmd.cf.ac.uk/ac/gene.php?gene=%s"%gene)
    query_result = br.response().read()
    soup = BeautifulSoup(query_result)
    try:
        return soup.find("input", {"name":"refcore"})['value']
    except:
        return "ERROR"



if __name__ == '__main__':
    br = initialize_browser()
    br = login_hgmd(br)

    cardiac = pd.read_csv("1_cardiac_genes_unq.csv",header=None)
    cardiac.columns = ['gene']

    refseq = []
    for g in cardiac['gene'].values:
        ref = browse_get(br, g)
        print(g,ref)
        refseq.append([g,ref])
        time.sleep(2)
