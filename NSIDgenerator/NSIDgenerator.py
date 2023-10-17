#!/usr/bin/python3
import requests
import urllib3
import sys
from xml.etree.ElementTree import XML
from xml.etree import ElementTree
import mygene
import pandas as pd
import numpy as np
import os
import pyodbc as db
import mysql.connector
from mysql.connector import Error

# Initialize an HTTP pool manager for making HTTP requests.
http = urllib3.PoolManager()

# Initialize the MyGeneInfo object.
mg = mygene.MyGeneInfo()

# Initialize an empty string to store search results.
search_DB = ""

# Function to fetch information for a given field and term from the Genenames API.
def fetch(field, term):
    url = 'http://rest.genenames.org/fetch/' + field + '/' + term
    requestedTerm = requests.get(url)
    termDictionary = {}

    if requestedTerm.status_code == 200:
        result = requestedTerm.text
        root = ElementTree.fromstring(requestedTerm.text)

        if len(root) > 1 and len(root[1]) > 0:
            for each in root[1][0]:
                termDictionary[each.attrib['name']] = []
                if len(each) > 0:
                    for hit in each:
                        termDictionary[each.attrib['name']].append(hit.text)
                else:
                    termDictionary[a.attrib['name']].append(a.text)

    return termDictionary

# Function to search for a term in Genenames and retrieve relevant information.
def search(field, term):
    if field == '':
        url = 'http://rest.genenames.org/search/' + term
    else:
        url = 'http://rest.genenames.org/search/' + field + '/' + term

    request = requests.get(url, stream=False)
    scores = {'symbol': '', 'score': 0.0}

    if request.status_code == 200:
        result = request.text
        if result.find('hgnc_id') > -1:
            result = result[result.find('hgnc_id') + 14:]
            scores['hgnc_id'] = result[:result.find('<')]
            result = result[result.find('symbol') + 8:]
            scores['symbol'] = result[:result.find('<')]
            result = result[result.find('score') + 7:]
            scores['score'] = float(result[:result.find('<')])

    return scores

# Searches all HUGO fields for the best result.
def do_search(term):
    to_search = ['symbol', 'prev_symbol', 'alias_symbol', 'alias_name', 'prev_name']
    best_result = ''

    for field in to_search:
        if best_result == '':
            search_result = search(field, term)
            score = search_result['score']
            symbol = search_result['symbol']

            if term.find('-') > -1:
                search_result = search(field, term.replace('-', ''))
                score_nohyphen = search_result['score']
                symbol_nohyphen = search_result['symbol']

                if score_nohyphen > score:
                    best_result = symbol_nohyphen
                elif score > 0.0:
                    best_result = symbol
            else:
                if score > 0:
                    best_result = symbol

    result = fetch('symbol', best_result)
    return result

# Format the output for a given term and dictionary of information.
def format_output(orig_term, dictionary):
    symbol, name, prev_symbol, prev_name, alias_symbol, alias_name = '', '', '', '', '', ''
    to_write = orig_term + '\t'

    if 'symbol' in dictionary:
        symbol = ','.join(dictionary['symbol'])
    if 'name' in dictionary:
        name = ','.join(dictionary['name'])
    if 'prev_symbol' in dictionary:
        prev_symbol = ','.join(dictionary['prev_symbol'])
    if 'prev_name' in dictionary:
        prev_name = ','.join(dictionary['prev_name'])
    if 'alias_symbol' in dictionary:
        alias_symbol = ','.join(dictionary['alias_symbol'])
    if 'alias_name' in dictionary:
        alias_name = ','.join(dictionary['alias_name'])

    return symbol, name, prev_symbol, prev_name, alias_symbol, alias_name

# Get gene identities
to_query = []
file = open(sys.argv[1], 'r', encoding="UTF-8")

initial_query = []
to_query = []

# Import all genes for query
for line in file:
    line = line.strip()
    if line != '' and line != ' ' and line not in to_query:
        orig_line = ''.join(line)
        if line.find('/') > -1 and line.find(' ') == -1:
            line = line[:line.find('/')]
        if line.find('(') > -1:
            line = line[:line.find('(')]
        to_query.append(line.strip())
        initial_query.append(orig_line)

file.close()
query_set = set(initial_query)

# Query all genes using MyGene
df_names = mg.querymany(query_set, scopes='symbol', species='Human', as_dataframe=True)
print(df_names)

# Check to make sure there are actual results
if "_score" in df_names:
    df_names = df_names.reset_index()
    df_names = df_names[df_names.index.isin(df_names.groupby('query')._score.idxmax())]

df_names.to_csv('result_csv.csv', index=False)

# Dictionary of customer preferred alias and official gene symbol pairs
customer_gene_dic = {}

# Process each query
for idx in range(len(to_query)):
    initial = initial_query[idx]
    query = to_query[idx]

    if "_score" not in df_names or query not in df_names['query'].values or df_names.loc[
        df_names['query'] == query, 'symbol'].iloc[0] != query:
        result_dict = do_search(query)
        if result_dict:
            symbol, name, prev_symbol, prev_name, alias_symbol, alias_name = format_output(orig_line, result_dict)
            search_DB = search_DB + ("'" + symbol + "',")
            customer_gene_dic[symbol] = query
        else:
            print("Could not find symbol for: " + query)
            search_DB = search_DB + ("'" + query + "',")
    else:
        continue

# Fotmat SQL query
for each in (df_names['symbol']):
    search_DB = search_DB + ("'" + each + "',")
search_DB = search_DB.rstrip(search_DB[-1])

# Build SQL query to select optimal probe
top_hits = "SELECT DISTINCT gene, accession_nover, probepair_nsid, TVHitCount, usecount, source, design_remarks FROM (SELECT * FROM ((SELECT  a.gene, a.accession_nover, a.probepair_nsid, COUNT(DISTINCT b.tv_id) + 1 AS TVHitCount, COUNT(DISTINCT c.library_part_number) AS usecount,'INV' AS source, a.design_remarks FROM v_probepair_data a LEFT JOIN v_library_probes c ON a.probepair_id = c.probepair_id LEFT JOIN v_transvar b ON a.probepair_nsid = b.probepair_nsid WHERE a.active = 1 AND a.design_remarks NOT LIKE ('%mismatch%') AND a.design_remarks NOT LIKE ('%no longer%') AND (a.probepair_nsid LIKE 'NM_%' OR a.probepair_nsid LIKE 'NR_%' OR a.probepair_nsid LIKE 'XR_%' OR a.probepair_nsid LIKE 'XM_%') AND organism_id = 1 AND  a.gene IN (" + search_DB + ") GROUP BY a.probepair_nsid ORDER BY gene ASC) UNION (SELECT  a.gene, a.accession_nover, a.probepair_nsid, COUNT(DISTINCT b.tv_id) + 1 AS TVHitCount, 0 AS usecount, 'DTA' AS source, '' AS design_remarks FROM v_DTA_probepair_data a LEFT JOIN v_DTA_transvar b ON a.probepair_nsid = b.probepair_nsid WHERE a.active = 1 AND (a.probepair_nsid LIKE 'NM_%' OR a.probepair_nsid LIKE 'NR_%' OR a.probepair_nsid LIKE 'XR_%' OR a.probepair_nsid LIKE 'XM_%') AND organism_id = 1 AND  a.gene IN (" + search_DB + ") GROUP BY a.probepair_nsid ORDER BY gene ASC))) AS x GROUP BY probepair_nsid ORDER BY gene ASC, TVHitCount DESC, usecount DESC) AS y GROUP BY gene"

# Initialize empty list to add selected probes
probes = []

# Query internal database with try catch to handle errors
try:
    connection = mysql.connector.connect(host='{HOST NAME HERE}', database="{DATABASE ID HERE}", user='{USERNAME}', password='{PASSWORD}')
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

# Function to map gene symbols to customer-defined symbols.
def mapper(x):
    value = customer_gene_dic.get(x)
    if value:
        return value
    else:
        return x

# Write file needed for pipeline input
libName = sys.argv[2]
sqNumber = sys.argv[3] 
inter = probes[[2,0]].copy()
inter[0] = inter[0].map(mapper)
NSID = inter.rename(columns={2:'params=' + libName, 0: sqNumber + ';'})
NSID.to_csv('NISD-' + libName + '-'+ sqNumber + '.csv', index=False)
print(customer_gene_dic)
