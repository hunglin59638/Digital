#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 21:47:36 2021

@author: hunglin
"""
import os
import json
import tempfile
tmp_dir = tempfile.gettempdir()
from Bio import SeqIO
from pathlib import Path
import subprocess
import tempfile
from modules.ident_prot import diamond
from modules.conv_id import GO_api, KEGG_api

"""
root_dir = Path(__file__).parent.parent
fasta = root_dir / "test/test.faa"
threads = str(os.cpu_count())
prefix = "out"
out_dir = root_dir / "test/output"
db = Path().home() / "databases/diamond/nr.dmnd"
force = True
"""
def main():
    def get_argument():
        import argparse
        def check_file(path):
            try:
                path = Path(path).resolve()
                if not path.is_file():
                    raise FileNotFoundError("No such a file")
                else:
                    return path
            except:
                raise TypeError("Please input path")
                
        parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("--fasta", help="input protein fasta", 
                            required=True, type=check_file)
        parser.add_argument("--nr_dmnd", help="NR database with diamond format", 
                            required=True, type=check_file)
        parser.add_argument("--out_dir","-o", help="output directory", 
                            default=os.getcwd())
        parser.add_argument("--prefix", help="prefix name of output files")
        parser.add_argument("--force", "-f", help="overwrite the exisiting files",
                            action="store_true")
        parser.add_argument("--emapper_path", help="path of eggnog-emapper",
                            default="emapper.py")
        parser.add_argument("--data_dir", help="eggnog-emapper database")
        parser.add_argument("--diamond_path", help="path of diamond",
                            default="diamond")
        args = parser.parse_args()
        if not args.prefix:
            args.prefix = os.path.basename(args.fasta).split(".")[0]
        args.threads = str(args.threads)
        args.out_dir = Path(args.out_dir)
        return args
        
    args = get_argument()
    fasta = Path(args.fasta)
    out_dir = Path(args.out_dir)
    threads = args.threads
    prefix = args.prefix
    force = args.force
    nr_dmnd = args.nr_dmnd
    data_dir = args.data_dir
    
    # identification
    total_ids = list(SeqIO.index(str(fasta), "fasta").keys())
    d = diamond()
    daa_file = out_dir / "diamond.txt"
    if daa_file.is_file() and force:
        daa_file.unlink()
    outfmt= ["6", "qseqid", "stitle", "scovhsp", "pident", "bitscore"]
    d.blastp(query=fasta, db=nr_dmnd, outfmt=outfmt, threads=threads,
             out=daa_file, more_sensitive=True,
             block_size=4, tmpdir="/dev/shm", quiet=True,
             id=0.8, query_cover=0.8, max_target_seqs=1)
    
    ident_dict = {}
    with open(f"{out_dir}/diamond.txt", "r") as f:
        while True:
            row = f.readline()
            if not row:
                break
            qseqid, stitle, scovhsp, pident, bitscore = row.rstrip("\n").split("\t")
            sseqid = stitle.split(' ')[0]
            product = ' '.join(stitle.split(" ")[1:])
            bitscore = float(bitscore)
            value = {"stitle": stitle, "stitle": stitle, "sseqid": sseqid,
                     "scovhsp": scovhsp, "product":product, "pident": pident,
                     "bitscore": bitscore}
            if qseqid not in ident_dict.keys():
                ident_dict[qseqid] = value    
            else:
                if ident_dict[qseqid]["bitscore"] < bitscore:
                    ident_dict[qseqid].update(value)

    seqid2name = {}
    for k, v in ident_dict.items():
        seqid2name[k]= v['product']
    mapped_ids = list(seqid2name.keys())
    for _id in total_ids:
        if _id not in mapped_ids:
            seqid2name[_id] = "hypothetical protein"
            
    # annotation
    annot_dict = emapper(fa=fasta,
                         output=f"{out_dir}/{prefix}",
                         threads=threads, mode="diamond",
                         dmnd_db=data_dir / "eggnog_proteins.dmnd")
    for key in annot_dict.keys():
        annot_dict[key].update({"NR": seqid2name[key]})
    go_list = set()
    ko_list = set()
    for key, value in annot_dict.items():
        go_ids = value["GOs"]
        go_list = go_list.union(go_ids)
        ko_ids = value["KEGG_ko"]
        ko_list = ko_list.union(ko_ids)
    if "-" in go_list:
        go_list.remove("-")
    if "-" in ko_list:
        ko_list.remove("-")
    go_api = GO_api()
    go_result = go_api.id2complete(ids=list(go_list))
    kegg_api = KEGG_api()
    ko2name = kegg_api.loop(func="list", entries=list(ko_list))
    ko2name= kegg_api.conv_format(query=ko2name, outfmt="dict")
    with open(f"{out_dir}/ko2name.json", "w") as f:
        json.dump(ko2name, f)
    with open(f"{out_dir}/go_result.json", "w") as f:
        json.dump(go_result, f)
    with open(f"{out_dir}/{prefix}_annotation.json", "w") as f:
        json.dump(annot_dict, f)
    emapper_dir =  (out_dir / "emapper")
    if not emapper_dir.is_dir():
        emapper_dir.mkdir()
    os.system(f"mv {out_dir}/*emapper.* {emapper_dir}")
        
    return 0

def emapper(fa, output, threads=os.cpu_count(), mode="diamond",
            temp_dir="/dev/shm",
            dmnd_db=None, emapper_path="emapper.py", **kwargs):

    try:
        emapper_path = subprocess.check_output(["which", emapper_path]).decode().\
            rstrip("\n")
    except:
        if not Path(emapper_path).is_file():
            raise FileNotFoundError(f"No such a file: {emapper_path}")
    dmnd_db = Path(dmnd_db) if dmnd_db else \
        Path(emapper_path).parent / "data/eggnog_proteins.dmnd"
    if not dmnd_db.is_file():
        raise FileNotFoundError(f"No such a file: {dmnd_db}")
    
    cmd = [emapper_path, "-i", fa, "-o", output, "--cpu",
           str(threads), "-m", mode, "--no_annot", 
           "--temp_dir", temp_dir]
    if dmnd_db != None and mode == "diamond":
        cmd.extend(["--dmnd_db", dmnd_db])

    for arg, value in kwargs.items():
        arg = "--" + arg if len(arg) > 1 else "-" + arg
        if type(value) == bool:
            if value:
                cmd.append(arg)
        elif type(value) == str:
            cmd.extend([arg, value])
    subprocess.call(cmd)

    cmd = [emapper_path, "-i", fa, "--annotate_hits_table",
           str(output) + ".emapper.seed_orthologs",
           "-o", output, "--cpu", str(threads),"--pfam_realign","denovo",
           "--override"]
    subprocess.call(cmd)
    result = {}
    with open(output + ".emapper.annotations", "r") as f:
        content = [x for x in f.read().rstrip("\n").split("\n")]

    for row in content:
        if row.startswith("# "):
            pass
        elif row.startswith("#"):
            cols = row.strip("#").split("\t")
        else:
            row = row.split("\t")
            result[row[0]] = {}
            for col, value in zip(cols, row):
                if "," in value:
                    value = value.split(",")
                result[row[0]].update({col: value})
    return result
