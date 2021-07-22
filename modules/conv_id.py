# -*- coding: utf-8 -*-

import re
from modules.web_api import go, kegg

def ko2name(ko_ids):
    kegg_api = kegg()
    ko2name_d = {}
    response = kegg_api.loop(func="list", entries=list(ko_ids))
    for row in [i for i in response.split("\n") if i ]:
        ko_id, names = row.split("\t")
        name = re.sub(r'\[EC:.{1,3}\..{1,3}\..{1,3}\..{1,3}]', "", 
                      ';'.join(names.split(";")[1:])).strip()
        ko2name_d[ko_id] = name
    return ko2name_d
        
def ko_path2name(ko_paths):
    kegg_api = kegg()
    ko_path2name_d = {}
    response = kegg_api.loop(func="list", entries=list(ko_paths))
    for row in [i for i in response.split("\n") if i ]:
        path, name = row.split("\t")
        path = re.sub(r'^path:', '', path)
        ko_path2name_d[path] = name
    return ko_path2name_d

def ko_reaction2name(ko_reactions):
    kegg_api = kegg()
    ko_reaction2name_d = {}
    response = kegg_api.loop(func="list", entries=list(ko_reactions))
    for row in [i for i in response.split("\n") if i ]:
        reaction, name = row.split("\t")
        reaction = re.sub(r'^rn:', "", reaction)
        ko_reaction2name_d[reaction] = name
    return ko_reaction2name_d

def ko_rclass2name(ko_rclasses):
    kegg_api = kegg()
    ko_rclass2name_d = {}
    response = kegg_api.loop(func="list", entries=list(ko_rclasses))
    for row in [i for i in response.split("\n") if i ]:
        rclass, name = row.split("\t")
        ko_rclass2name_d[rclass] = name
    return ko_rclass2name_d

def ko_module2name(ko_modules):
    kegg_api = kegg()
    ko_module2name_d = {}
    response = kegg_api.loop(func="list", entries=list(ko_modules))
    for row in [i for i in response.split("\n") if i ]:
        module, name = row.split("\t")
        module = re.sub(r'^md:', "", module)
        ko_module2name_d[module] = name
    return ko_module2name_d

def brite2name(ko_brs):
    kegg_api = kegg()
    brite2name_d = {}
    response = kegg_api.loop(func="list", entries=[f"br:{br}"for br in list(ko_brs)])
    for row in [i for i in response.split("\n") if i ]:
        br, name = row.split("\t")
        br = re.sub(r"^br:", "", br)
        brite2name_d[br] = name
    return brite2name_d

def convert_cazy(cazy_ids):
    pass

def go2name(go_ids):
    go_api = go()
    go_d = {}
    go_result = go_api.id2complete(ids=list(go_ids))
    for key, value in go_result.items():
        go_d[key] = {"name": value["name"], "ontology": value["aspect"]}
    return go_d