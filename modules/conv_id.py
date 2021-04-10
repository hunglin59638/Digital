# -*- coding: utf-8 -*-
"""
id convert and annotation
"""

import pandas as pd
import sys,re, os, platform
import json
import requests
import urllib.parse
import tempfile
tmp_dir = tempfile.gettempdir()
from tqdm import tqdm
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET


class GO_api:
    def __init__(self):
        self.url= "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
    
    def id2complete(self, ids):
        
        if type(ids) == list or set:
            ids = list(ids)
        elif type(ids) == str:
            ids= [ids]
        else:
            raise TypeError("GO id format is required for str, set or list")
            
        batch= 500
        i=0
        results= []
        pbar= tqdm(total= len(ids), desc= "convert GO ids")
        while True:
            chunk= ids[i*batch:(i+1)*batch]
            i +=1

            if not chunk:
                break
            
            idStr = ",".join(chunk)
            requestURL = self.url + urllib.parse.quote(idStr) + "/complete"
            r = requests.get(requestURL, headers={ "Accept" : "application/json"})
            
            if not r.ok:
              print(r.raise_for_status(), file=f"{tmp_dir}/quickgo.log")
              sys.exit()
            
            responseBody = r.text
            json_dict_go = json.loads(responseBody)['results']
            results.extend(json_dict_go)
            pbar.update(len(chunk))
        pbar.close()
        
        final = {}
        for row in results:
            if row["id"] in ids:
                final[row["id"]] = row
            if "secondaryIds" in row.keys():
                for x in row["secondaryIds"]:
                    if x in ids:
                        final[x] = row
            
        return final

class KEGG_api:
    from Bio.KEGG import REST

    def entry(self, entries):
        # process query ids
        if type(entries) == str and "+" not in entries:
            entries = [entries]
        elif type(entries) == str and "+" in entries:
            entries= [x for x in entries.split("+") if x != ""]
        elif type(entries) == list:
            pass
        return entries
    
    def loop(self, func, entries,target_db = None ,chunk= 100):
        
        def request(entries,func,target_db = None, url= "http://rest.kegg.jp/"):
            import urllib.request
            from time import sleep
            url= url+"/"+func+"/"
            if func in ["list","get"]:
                url = url + entries
            else:
                url = url + target_db +"/" + entries
            
            while True:
                try:
                    with urllib.request.urlopen(url) as response:
                        html = response.read()
                    return html
                except:
                    sleep(10)
            
        
        if chunk > 100:
            return 1
        
        result = ""
        loops = int(len(entries)/chunk) + 1
        pbar = tqdm(total=loops, desc="convert kegg ids: ")
        for l in range(loops):
            entry = entries[l*chunk: l*chunk+chunk]
            entry = "+".join(entry)
            if not entry:
                break
            if func == "list":
                try:
                    x = self.REST.kegg_list(database= entry).read()
                except:
                    x = request(entries= entry, func= func)
            elif func == "conv":
                try:
                    x= self.REST.kegg_conv(target_db= target_db, source_db=entry).read()
                except:
                    x = request(entries= entry, target_db= target_db)
            elif func == "link":
                try:
                    x = self.REST.kegg_link(target_db, source_db= entry).read()
                except:
                    x = request(entries= entry, target_db= target_db, func= func)
            elif func == "get":
                try:
                    x = self.REST.kegg_get(dbentries= entry)
                except:
                    x= request(entries= entry, func= func)
            result= result + x
            pbar.update(1)
        pbar.close()
        return result
    
    def conv_format(self, query, outfmt= "dict"):
        if outfmt == "dict":
            result = {}
            for row in [x for x in query.split("\n") if x!= ""]:
                row= row.split("\t")
                entry= row[0]
                target= row[1]
                if entry not in result.keys():
                    result[entry]= [target]
                elif target not in result[entry]:
                    result[entry].append(target)
        elif outfmt== "tab":
            result = query
        return result
        
        
    def list_kegg(self, entries):
        entries= self.entry(entries= entries)
        result= self.loop(func= "list", entries= entries)
        result= self.conv_format(query= result, outfmt = "dict")
        return result
    
    def conv(self, target_db, entries):
        try:
            result = self.REST.kegg_conv(target_db, entries).read()
        except:
            pass
        
        
if platform.system() == "Windows":
    hostname = os.environ['COMPUTERNAME']
elif platform.system() == "Linux":
    hostname = os.uname()[1]

def uniqList(List):  
    tmp = []
    for x in List:
        if x in tmp:
            pass
        else:
            tmp.append(x)
    return tmp

class idMapKit:
    import requests
    def __init__(self):
        self.email = "hunglin59638@gmail.com"
        self.batch_size = 100
        self.uniprot_regex = r'([O,P,Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9]|([A-N,R-Z][0-9][A-Z][A-Z,0-9]{2}[0-9])([A-Z][A-Z,0-9]{2}[0-9]){0,1})\.{0,1}[1-9]{0,1}'
        self.cog_regex = r'[C,K,P]OG[0-9]{4}'
    """  
    with open(jsonPath, 'r', encoding='utf-8') as f:
        merge_fileds = json.load(f)
    """  
    def __fieldTrans(self, field = "uniprot", dict_type = "mygene"):
        try:
            return self.merge_fileds[field][dict_type]
        except KeyError:
            print(field, " is not found.")
            return field
    def sort_by_list(self, df, List):
        tmp = pd.DataFrame(columns= df.columns)
        for i in List:
            if i in list(df["query"]):
                tmp = tmp.append(df[df["query"] == i], ignore_index= True)
        return tmp
        
    def mygene(self, query = [], fr = "accession_prot", to = ["uniprot", "name"], as_dataframe = False,returnall = False):
        from biothings_client import get_client
        mg = get_client('gene')
        fr = self.__fieldTrans(field = fr, dict_type = "mygene")
        to = [self.__fieldTrans(field = x, dict_type = "mygene") for x in to ]
        return mg.querymany(query , scopes = fr, fields = to,
                     as_dataframe = as_dataframe, returnall= returnall)
    

    def uniprotMap(self, fr= "accession_prot", to= ["uniprot","name","ec_number"], * , 
                   query = [],as_dataframe = False, returnall = False):
        import urllib.parse
        import urllib.request
        import urllib.error
        if query == []:
            return
        url = 'https://www.uniprot.org/uploadlists/'
        list_str = " ".join(query)
        fr = self.__fieldTrans(fr, "uniprot_map")
        params = {
        'from': fr,
        'to': "ACC",
        'format': 'tab',
        'query': list_str
        }
        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        try:
            with urllib.request.urlopen(req) as f:
                yourlist = f.geturl()[32:-4]
                f.close()
        except urllib.error.HTTPError or urllib.error.URLError:
            print("URL or Server error")
            return None
            
        columns = ",".join([self.__fieldTrans(x, "uniprot_col") for x in to])
        url = "https://www.uniprot.org/uniprot/?query=yourlist%3A" + yourlist + \
             "&format=tab&columns=yourlist(" + yourlist + ")," + columns
        print("yourlist url: ",url)
        tab_result = self.requests.get(url).text.split("\n")
        df = pd.DataFrame(columns = tab_result[0].split("\t"))
        for raw in tab_result:
            if "yourlist" not in raw and raw != "":
                df_add = pd.DataFrame([raw.split("\t")], columns = tab_result[0].split("\t"))
                df = df.append(df_add, ignore_index=True)
        df = df.rename(columns = {df.columns[0]: "query","Entry":"Uniprot"})
        
        muti_query = [i for i in df["query"] if len(i.split(",")) > 1 ]
        for x in muti_query:
            tmp = df[df["query"]== x]
            for i in range(len(x.split(","))-1):
                tmp = tmp.append(tmp, ignore_index = True)
            tmp["query"] = x.split(",")
            df = df.append(tmp, ignore_index= True)
            df = df[df["query"] != x]
        
        try:
            missing = [x for x in query if x not in ",".join(list(df["query"])).split(",")]
        except:
            missing = []
        if returnall == True:
            df["not found"] = [False] * len(df)
            mis_df = pd.DataFrame({"query": missing,"not found": [True] * len(missing)},
                                  columns = list(df.columns))
            df = df.append(mis_df, ignore_index=True)
        else:
            pass
        
        if as_dataframe == True:
            return self.sort_by_list(df,query)
        else:
            df = self.sort_by_list(df,query)
            to_dict = df.to_dict("index")
            if returnall == True:
                tmp = {"out":{}, "missing": []}
                for key, value in to_dict.items():
                    if value["not found"] == False:
                        del value["not found"]
                        if value["query"] in tmp.keys():
                            tmp["out"][value["query"]].update(value)
                        else:
                            tmp["out"].update([(value["query"],value)])
                    else:
                        del value["not found"]
                        tmp["missing"].append(value["query"])
            else:
                tmp = {}
                for key,value in to_dict.items():
                    if value["query"] in tmp.keys():
                        tmp[value["query"]].update(value)
                    else:
                        tmp.update([(value["query"],value)])                    
            return tmp
        
    
    def ncbiSearch(self, query = [], fr = ["accession_prot" , "uniprot"] , 
                   to= ["name" ,"ec_number", "uniprot"],db = "protein"):
        from Bio import Entrez
        import urllib.error
        if type(fr) == str:
            fr = fr.split()
        Entrez.email = self.email
        #queryList = query.copy()
        count = len(query)
        query = ",".join(query)
        search_results = Entrez.read(Entrez.epost(db, id= query ))
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]
        print("webenv: ",webenv+"\n")
        print("qerykey: ",query_key+"\n")
        batch_size = self.batch_size
        out_handle = open("tmp.xml", "w")
        def makeXmlRoot():
            fetch_handle = Entrez.efetch(db= db, rettype="gp",
                            retmode="xml", retstart= start, retmax= batch_size,
                             webenv= webenv, query_key= query_key)
            data = fetch_handle.read()
            root = ET.fromstring(data)
            return fetch_handle, data, root
        #for start in range(0, count, batch_size):
        start = 0
        result = {}
        pbar = tqdm(total=100)
        while True:
            end = start + batch_size
            #print("Going to download record %i to %i" % (start+1, end))
            if end <= count:
                try:
                    pbar.update(count/batch_size)
                except:
                    pass            
            try:
                fetch_handle = makeXmlRoot()[0]
                data = makeXmlRoot()[1]
                root = makeXmlRoot()[2]
            except urllib.error.HTTPError:
                print("finished")
                break
            for subroot in root:
                """
                subroot represents a result
                """
                tmp = {}
                for i in to:
                    tmp[i] = []
                if "uniprot" in fr or "accession_prot" in fr:
                    for elem in subroot.iter("GBSeq_accession-version"):
                        queryid = elem.text
                    tmp['query'] = queryid
                if "accession_prot" in to or "uniprot" in to:
                    for elem in subroot.iter("GBSeqid"):
                        value = elem.text
                        if "gb|" in value or "ref|" in value:
                            tmp["accession_prot"].append(value.split('|')[1])
                        if "sp|" in value:
                            tmp["uniprot"].append(value.split('|')[1])
                for elem in subroot.iter("GBQualifier"):
                    for subelem in elem:
                        if subelem.tag == "GBQualifier_name":
                            name = subelem.text
                        if subelem.tag == "GBQualifier_value":
                            value = subelem.text
                    if "entrezid" in to and name == "db_xref" and "GeneID" in value:
                        tmp["entrezid"].append(re.search(r'[0-9]+',value)[0])
                    if "name" in to and name == "product":
                        tmp["name"].append(value)
                    if "ec_number" in to and name == "EC_number":
                        tmp["ec_number"].append(value)
                    if "uniprot" in to and "UniProtKB" in value:
                        try:
                            uniprot = re.search(self.uniprot_regex, value)[0]
                        except:
                            print(value)
                        #print(uniprot)
                        if uniprot not in tmp["uniprot"]:
                            tmp["uniprot"].append(uniprot)
                    if "taxid" in to and "taxon" in value:
                        tmp["taxid"].append(re.search(r'[0-9]+',value)[0])
                    if "org_name" in to and name == "organism":
                        tmp['org_name'].append(value)
                    if "symbol" in to and name == "gene":
                        tmp['symbol'].append(value)
                    if "pfam" in to and name == "note" and "pfam" in value:
                        if value not in tmp['pfam']:
                            tmp['pfam'].append(value)
                        
                if "sequence_prot" in to:
                    for elem in subroot.iter("GBSeq_sequence"):
                        seq = elem.text
                        tmp["sequence_prot"].append(seq)
                        
                if "go" in to or "accession_prot" in to or "kegg" in to or \
                "pfam" in to:
                    for elem in subroot.iter("GBSeq_source-db"):
                        value = elem.text
                        if "go" in to and "GO:" in value:
                            tmp["go"] = re.findall(r'GO:[0-9]+',value)
                        if "accession_prot" in to and "xrefs:" in value:
                            accession_prot = re.search(r'xrefs:.*?;',value)[0].strip('xrefs:').strip(";").replace(" ","").split(",")
                            if tmp['accession_prot'] == []:
                                tmp['accession_prot'].append(accession_prot)
                        if "kegg" in to and "KEGG" in value:
                            kegg = re.findall(r'KEGG:[a-z]{3}:[0-9]+',value)
                            kegg = [x[5:] for x in kegg]
                            tmp["kegg"] = kegg;del kegg
                        if "pfam" in to and "Pfam:" in value:
                            pfam = re.findall(r'(PF\d{5}|CL\d{4})',value)
                            tmp['pfam'] = pfam;del pfam
                        if "string" in to and "STRING:" in value:
                            string = re.findall(r'STRING:.+, ',value)
                            string = [x[7:-1] for x in string]
                        if "interpro" in to and "InterPro":
                            interpro = re.findall(r'IPR[0-9]{6}',value)
                            tmp['interpro'].append(interpro)
                if "uniprot" in fr:
                    for elem in subroot.iter("GBSeq_source-db"):
                        pass
                    #Pfam, eggNog...
                    
                if queryid not in result.keys():
                    "this query id input to result dict first"
                    result[queryid] = tmp 
                else:
                    if type(result[queryid]) == dict:
                        "this query id only had one result in dict"
                        result[queryid] = [result[queryid]]
                        result[queryid].append(tmp)
                    elif type(result[queryid]) ==  list:
                        "more than one result in dict"
                        result[queryid].append(tmp)
                
            fetch_handle.close()
            out_handle.write(str(data))
            start = start + batch_size            
        out_handle.close()
        pbar.close()
        #os.remove("tmp.xml")
        return result