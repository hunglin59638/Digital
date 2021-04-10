#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 20:36:54 2020

@author: hunglin
"""
import subprocess
import os
import re
from termcolor import colored
from pathlib import Path

class diamond:
    def __init__(self):
        self.diamond = Path(self.which_diamond())
    
    def which_diamond(self, path = None):
        if path == None:
            try:
                diamond = subprocess.check_output(["which","diamond"]).decode().\
                    rstrip("\n")
            except:
                diamond = None
        else:
            diamond = path
        
        self.diamond = Path(diamond)
        return self.diamond
            
    def makedb(self, faa, db, taxonmap = None, taxonnodes= None, taxonnames = None):
        "Create a DIAMOND formatted reference database from a FASTA input file."
        if self.diamond == None:
            print("diamond is not found")
            return 1
        cmd = [self.diamond, "makedb", "--in", faa, "--db", db] # basic command
        
        if taxonmap != None:
            cmd.extend(["--taxonmap", taxonmap])
        if taxonnodes != None:
            cmd.extend(["--taxonnodes", taxonnodes])
        if taxonnames != None:
            cmd.extend(["--taxonnames", taxonnames])
        
        check = subprocess.check_call(cmd)
        if not check:
            print("make diamond database failed")
            
    def download_taxdb(self, save_path,db ,decompress = True):
        if db == "accession2protein":
            ftp = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"
        elif db == "taxdmp":
            ftp = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
        md5 = ftp + ".md5"
        md5_value = subprocess.check_output(["curl", md5]).decode()[0]

        old_md5 = save_path+os.path.basename(md5)
        print(os.path.exists(old_md5))
        if os.path.exists(old_md5):
            old_md5 = subprocess.check_output(["cat", old_md5]).decode().split(" ")[0]
            check_lastest = True if old_md5 != md5_value else False
        else:
            check_lastest = False
        if not check_lastest:
            db_file = save_path + os.path.basename(ftp)
            cmd = ["wget","-c", ftp, "-O", db_file]
            success = subprocess.check_output(cmd)
            cmd[-1] = save_path+ os.path.basename(md5)
            subprocess.call(cmd)
            if decompress and success == 0:
                if db == "accession2protein":
                    cmd = ["pigz","-p", str(os.cpu_count()),"-d", db_file]
                elif db == "taxdmp":
                    cmd = ["unzip", db_file]
                subprocess.call(cmd)
                
    def blastp(self, query, db, out=None ,outfmt= 100, threads= os.cpu_count(), 
                **kwargs):
        outfmt = [str(i) for i in outfmt] if isinstance(outfmt, (list, tuple)) else [str(outfmt)]
        cmd = [self.diamond, "blastp", "--db", db, "--query", query, "--out", 
               out, "--outfmt"]
        cmd.extend(outfmt)
        if out:
            if not Path(out).parent.is_dir():
                Path(out).parent.mkdir()
        for arg, value in kwargs.items():
            arg= arg.replace("_","-")
            arg= "--{}".format(arg) if len(arg) >1 else "-{}".format(arg)
            if isinstance(value, bool):
                if value:
                    cmd.append(arg)
            else:
                cmd.extend([arg, str(value)])
        if not set(["--block-size", "-b"]).intersection(set(cmd)):
            ram_size = subprocess.check_output(["cat", "/proc/meminfo"]).decode().split("\n")[0]
            ram_size = int(re.search(r'[0-9]+', ram_size)[0]) / (1024*1024)
            if ram_size < 16 and ("block_size" not in kwargs.keys() or "b" not in kwargs.keys()):
                cmd.extend(["--block-size", "1"])
        cmd = [str(x) for x in cmd]
        print(colored("CMD: {}".format(" ".join(cmd)), "yellow", attrs=["bold"]))
        success = subprocess.check_call(cmd)
        
        if success:
            with open(out, "r") as f:
                result = f.read()
            if "6" in str(outfmt).split(" "):
                if str(outfmt) == "6":
                    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
                blast_result = {}
                fields = [x for x in str(outfmt).split(" ") if x != "6"]
                count = 0
                for row in [x for x in result.split("\n") if x != ""]:
                    row= row.split("\t")
                    blast_result[count] = {}
                    for i in range(len(fields)):
                        blast_result[count].update({fields[i]:row[i]})
                    count +=1
                tmp = {}
                if "qseqid" in fields and "pident" in fields:
                    for i in range(len(blast_result)):
                        query_id = blast_result[i]["qseqid"]
                        if query_id not in tmp.keys():
                            tmp.update({query_id: blast_result[i]})
                        else:
                            if blast_result[i]['pident'] > tmp[query_id]["pident"]:
                                tmp.update({query_id: blast_result[i]})
                result = tmp
            elif str(outfmt) == "102":
                tmp = {}
                for row in [x for x in result.split("\n") if x != ""]:
                    row = row.split("\t")
                    tmp[row[0]] = {
                        "query_id": row[0],
                        "taxid" : row[1],
                        "e_value": row[2]
                        }
                result = tmp
            
            return result
        
    def view(self, daa, out= None, outfmt=6, forwardonly= False, **kwargs):
        if isinstance(outfmt, (int, str)):
            outfmt = str(outfmt).split(" ")
        elif isinstance(outfmt, list):
            pass
        
        cmd= [self.diamond, "view", "--daa", daa, "--outfmt"]
        cmd.extend(outfmt)
        if out:
            cmd.extend(["--out", out])
        for arg, value in kwargs.items():
            arg= "--{}".format(arg) if len(arg) > 1 else "-{}".format(arg) 
            if isinstance(value, bool):
                if value:
                    cmd.append(arg)
            else:
                cmd.extend([arg, value])
        if not out:
            result= subprocess.check_output(cmd).decode().strip("\n")
            return result
        elif os.path.exist(out):
            return 0
        else:
            raise FileNotFoundError(out)
            
class hmmer:
    def __init__(self, path = None):
        self.hmmer_path = self.which_hmmer(path)
    
    def which_hmmer(self, path):
        if path== None:
            path = subprocess.check_output(["which","hmmbuild"]).decode().replace("/hmmbuild\n","")
        elif os.path.isfile(path):
            pass
        else:
            path = None
        return path
        
        

    
    
    
    
