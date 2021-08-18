#!/usr/bin/python3
import requests,urllib3,sys
from xml.etree.ElementTree import XML
from xml.etree import ElementTree
http=urllib3.PoolManager()
import mygene, pandas as pd, numpy as np
mg = mygene.MyGeneInfo()
import os

import pyodbc as db
import mysql.connector
from mysql.connector import Error

#function will be called later and take the gene symbol the best score for identity of the gene


search_DB = ""

def fetch(field,term):
    url = 'http://rest.genenames.org/fetch/'+field+'/'+term 
    requestedTerm = requests.get(url)
    termDictionary = {} 
    if requestedTerm.status_code == 200:
        result = requestedTerm.text
        root=ElementTree.fromstring(requestedTerm.text)
        if len(root)>1 and len(root[1])>0:
            for a in root[1][0]:
                termDictionary[a.attrib['name']]=[]
                if len(a)>0:
                    for e in a:
                        termDictionary[a.attrib['name']].append(e.text)
                else:
                    termDictionary[a.attrib['name']].append(a.text)
    return termDictionary



def search(field,term):
    if field == '':
        url = 'http://rest.genenames.org/search/'+term
    else:
        url = 'http://rest.genenames.org/search/'+field+'/'+term
    a = ''
    r = requests.get(url,stream=False)
    d = {}
    d['symbol']=''
    d['score']=0.0
    if r.status_code == 200:
        result = r.text
        if result.find('hgnc_id')>-1:
            result = result[result.find('hgnc_id')+14:]
            d['hgnc_id'] = result[:result.find('<')]
            result = result[result.find('symbol')+8:]
            d['symbol'] = result[:result.find('<')]
            result = result[result.find('score')+7:]
            d['score'] =result[:result.find('<')]
    return d

#Searches all HUGO fields for result
def do_search(term):
    to_search = ['symbol','prev_symbol','alias_symbol', 'alias_name', 'prev_name']
    best_result = ''
    for s in to_search:
        if best_result == '':
            a = search(s, term)
            score = float(a['score'])
            symbol = a['symbol']
            if term.find('-') > -1:
                a = search(s, term.replace('-', ''))
                score_nohyphen = float(a['score'])
                symbol_nohyphen = a['symbol']
                if score_nohyphen > score:
                    best_result = symbol_nohyphen
                elif score > 0.0:
                    best_result = symbol
            else:
                if score > 0:
                    best_result = symbol
    result = fetch('symbol',best_result)
   
    return result

def format_output(orig_term,d):
    symbol = ''
    name = ''
    prev_symbol = ''
    prev_name = ''
    alias_symbol = ''
    alias_name = ''
    to_write = orig_term+'\t'
    if 'symbol' in d:
        symbol = ','.join(d['symbol'])
    if 'name' in d:
        name = ','.join(d['name'])
    if 'prev_symbol' in d:
        prev_symbol = ','.join(d['prev_symbol'])
    if 'prev_name' in d:
        prev_name = ','.join(d['prev_name'])
    if 'alias_symbol' in d:
        alias_symbol = ','.join(d['alias_symbol'])
    if 'alias_name' in d:
        alias_name = ','.join(d['alias_name'])
    return symbol, name, prev_symbol, prev_name, alias_symbol, alias_name


#Get gene identities
to_query = []
f=open(sys.argv[1],'r',encoding="UTF-8")


initial_query = []
to_query = []
#Import all genes for query
for line in f:
    line=line.strip()
    if line !='' and line!=' ' and line not in to_query:
        orig_line = ''.join(line)
        if line.find('/')>-1 and line.find(' ')==-1:
            line = line[:line.find('/')]
        if line.find('(')>-1:
            line = line[:line.find('(')]
        to_query.append(line.strip())
        initial_query.append(orig_line)
f.close()
query_set = set(initial_query)


#Query all from MyGene
df_names = mg.querymany(query_set, scopes='symbol', species = 'Human',as_dataframe=True)
print(df_names)
#Check to make sure there actually was a result
if "_score" in df_names:
    df_names = df_names.reset_index()
    df_names = df_names[df_names.index.isin(df_names.groupby('query')._score.idxmax())]

df_names.to_csv('result_csv.csv', index=False)



customer_gene_dic = {}


for idx in range(len(to_query)):
    initial = initial_query[idx]
    query = to_query[idx]
    if "_score" not in df_names or query not in df_names['query'].values or df_names.loc[df_names['query']== query,'symbol'].iloc[0] != query:
        d = do_search(query)
        if d:
       
            symbol, name, prev_symbol, prev_name, alias_symbol, alias_name = format_output(orig_line,d)
            search_DB = search_DB + ("'" + symbol + "',")
            customer_gene_dic[symbol] = query 
        else:
            print("Could not find symbol for: " + query)
            search_DB = search_DB + ("'" + query + "',")
    else:
        continue



for each in (df_names['symbol']):

    search_DB = search_DB + ("'"+ each + "',")

search_DB = search_DB.rstrip(search_DB[-1])

top_hits = "select distinct gene, accession_nover, probepair_nsid, TVHitCount, usecount, source, design_remarks from( select * from ( (select  a.gene, a.accession_nover, a.probepair_nsid, count(distinct b.tv_id) + 1 as TVHitCount, count(distinct c.library_part_number) as usecount,'INV' as source, a.design_remarks from v_probepair_data a left join v_library_probes c on a.probepair_id = c.probepair_id left join v_transvar b on a.probepair_nsid = b.probepair_nsid where a.active = 1 and a.design_remarks not like ('%mismatch%') and a.design_remarks not like('%no longer%') and (a.probepair_nsid like 'NM_%' or a.probepair_nsid like 'NR_%' or a.probepair_nsid like 'XR_%' or a.probepair_nsid like 'XM_%') and organism_id = 1 and  a.gene in (" + search_DB +") group by a.probepair_nsid order by gene ASC ) UNION (select  a.gene, a.accession_nover, a.probepair_nsid, count(distinct b.tv_id) + 1 as TVHitCount, 0 as usecount, 'DTA' as source, '' as design_remarks from v_DTA_probepair_data a left join v_DTA_transvar b on a.probepair_nsid = b.probepair_nsid where a.active = 1 and (a.probepair_nsid like 'NM_%' or a.probepair_nsid like 'NR_%' or a.probepair_nsid like 'XR_%' or a.probepair_nsid like 'XM_%') and organism_id = 1 and  a.gene in (" + search_DB + ") group by a.probepair_nsid order by gene ASC ) )as x group by probepair_nsid order by gene ASC, TVHitCount DESC, usecount DESC ) as y group by gene"
probes = []



try:
    connection = mysql.connector.connect(host='Loki', database="BISProd", user='nanostring', password='user')
    if connection.is_connected():
        mycursor = connection.cursor()
        mycursor.execute(top_hits)

        probes = pd.DataFrame(mycursor.fetchall())
        print(probes)
       
        

except Error as e:
    print("Error while connecting to MySQL", e)
finally:
    if connection.is_connected():
        mycursor.close()
        connection.close()


def mapper(x):
    value = customer_gene_dic.get(x)
    if value:
        return value
    else:
        return x


#request_file = sys.argv[1]
#library_dir = os.path.abspath(request_file)
#directName = os.path.basename(os.path.dirname(library_dir))

libName = sys.argv[2]
sqNumber = sys.argv[3] 

inter = probes[[2,0]].copy()

inter[0] = inter[0].map(mapper)

NSID = inter.rename(columns={2:'params=' + libName, 0: sqNumber + ';'})

NSID.to_csv('NISD-' + libName + '-'+ sqNumber + '.csv', index=False)
print(customer_gene_dic)
