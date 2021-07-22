#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Jan  6 21:47:36 2021

@author: hunglin
"""
import os
import json
import shutil
from Bio import SeqIO
from copy import deepcopy
from pathlib import Path
import subprocess
import argparse
from modules.ident_prot import diamond
from modules.conv_id import *

def main():
    
    def get_argument():
        
        def check_file(path):
            if not path:
                raise TypeError("Please input path")
            else:
                path = Path(path)
                if not path.exists():
                    raise argparse.ArgumentTypeError("No such as a file or directory")
                    # raise FileNotFoundError("No such as a file or directory")
                else:
                    return path
            raise TypeError("Please input path")
                
        parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
                                         description="A pipeline to annotate protein sequences")
        parser.add_argument("--fasta", help="input protein fasta", 
                            required=True, type=check_file)
        parser.add_argument("--nr_dmnd", help="NR database with diamond format", 
                            required=False)
        parser.add_argument("--out_dir","-o", help="output directory", 
                            default=os.getcwd())
        parser.add_argument("--prefix", help="prefix name of output files")
        parser.add_argument("--threads","-t", help="cores of cpu", 
                            default=os.cpu_count())
        parser.add_argument("--force", "-f", help="overwrite the exisiting files",
                            action="store_true")
        parser.add_argument("--emapper_path", help="path of eggnog-emapper",
                            default="emapper.py")
        parser.add_argument("--data_dir", help="eggnog-emapper database",
                            type=check_file, required=True)
        parser.add_argument("--diamond_path", help="path of diamond",
                            default="diamond")
        args = parser.parse_args()
        if not args.prefix:
            args.prefix = os.path.basename(args.fasta).split(".")[0]
        args.threads = str(args.threads)
        args.out_dir = Path(args.out_dir)
        args.data_dir = check_file(args.data_dir)
        return args
        
    args = get_argument()
    fasta = Path(args.fasta)
    out_dir = Path(args.out_dir)
    threads = args.threads
    prefix = args.prefix
    force = args.force
    nr_dmnd = None if not args.nr_dmnd else Path(args.nr_dmnd)
    data_dir = args.data_dir
    emapper_path = args.emapper_path
    data_dir = data_dir if data_dir else Path(subprocess.check_output(["which", emapper_path]).decode().\
            rstrip("\n")).parent / "data"

    # identification
    total_ids = list(SeqIO.index(str(fasta), "fasta").keys())
    diamond_f = out_dir / "diamond.txt"
    if force:
        diamond_f.unlink(missing_ok=True)
    if (not diamond_f.is_file() or not diamond_f.stat().st_size) and nr_dmnd:
        diamond_f.unlink(missing_ok=True)
        print("Protein identification using diamond against nr database")
        d = diamond()
        daa_file = out_dir / "diamond.txt"
        if daa_file.is_file() and force:
            daa_file.unlink()
        outfmt= ["6", "qseqid", "stitle", "scovhsp", "pident", "bitscore"]
        d.blastp(query=fasta, db=nr_dmnd, outfmt=outfmt, threads=threads,
                 out=daa_file, more_sensitive=True,
                 block_size=4, tmpdir="/dev/shm", quiet=True,
                 id=0.8, query_cover=0.8, max_target_seqs=1)
    else:
        print(f"{diamond_f} is found and use it")
        
    if nr_dmnd:
        ident_dict = {}
        with open(diamond_f, "r") as f:
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
    if force:
        for file in out_dir.iterdir():
            if "emapper" in file.name:
                if file.is_file():
                    file.unlink()
                elif file.is_dir():
                    shutil.rmtree(file, ignore_error=True)
    print("Annotation using eggnog-emapper")
    annot_dict = emapper(fa=fasta,
                         output=f"{out_dir}/{prefix}",
                         threads=threads, mode="diamond",
                         dmnd_db=data_dir / "eggnog_proteins.dmnd")
    if nr_dmnd:
        for key in annot_dict.keys():
            annot_dict[key].update({"NR": seqid2name[key]})
    print("Converting GO and KEGG ids")
    go_ids = set()
    ko_ids = set()
    ko_paths = set()
    ko_brs = set()
    ko_modules = set()
    ko_rclasses = set()
    ko_reactions = set()
    for key, value in annot_dict.items():
        go_ids = go_ids.union(value["GOs"])
        ko_ids = ko_ids.union(value["KEGG_ko"])
        ko_paths = ko_paths.union(value["KEGG_Pathway"])
        ko_brs = ko_brs.union(value["BRITE"])
        ko_modules = ko_modules.union(value["KEGG_Module"])
        ko_rclasses = ko_rclasses.union(value["KEGG_rclass"])
        ko_reactions = ko_reactions.union(value["KEGG_Reaction"])
    go_d = go2name(go_ids)
    ko2name_d = ko2name(ko_ids)
    ko_path2name_d = ko_path2name(ko_paths)
    brite2name_d = brite2name(ko_brs)
    ko_module2name_d = ko_module2name(ko_modules)
    ko_rclass2name_d = ko_rclass2name(ko_rclasses)
    ko_reaction2name_d = ko_reaction2name(ko_reactions)
    
    for key, value in deepcopy(annot_dict).items():
        d = {}
        for go_id in value["GOs"]:
            go_d.setdefault(go_id, {})
            d[go_id] = go_d[go_id]
        annot_dict[key]["GOs"] = d
        d.clear()
        for ko_id in value["KEGG_ko"]:
            ko2name_d.setdefault(ko_id)
            d[ko_id] = ko2name_d[ko_id]
        annot_dict[key]["KEGG_ko"] = d
        d.clear()
        for ko_path in value["KEGG_Pathway"]:
            ko_path2name_d.setdefault(ko_path)
            d[ko_path] = ko_path2name_d[ko_path]
        annot_dict[key]["KEGG_Pathway"] = d
        d.clear()
        for br in value["BRITE"]:
            brite2name_d.setdefault(br)
            d[br] = brite2name_d[br]
        annot_dict[key]["BRITE"] = d
        d.clear()
        for md in value["KEGG_Module"]:
            ko_module2name_d.setdefault(md)
            d[md] = ko_module2name_d[md]
        annot_dict[key]["KEGG_Module"] = d
        d.clear()
        
        for rclass in value["KEGG_rclass"]:
            ko_rclass2name_d.setdefault(rclass)
            d[rclass] = ko_rclass2name_d[rclass]
        annot_dict[key]["KEGG_rclass"] = d
        d.clear()
        
        for reaction in value["KEGG_Reaction"]:
            ko_reaction2name_d.setdefault(reaction)
            d[reaction] = ko_reaction2name_d[reaction]
        annot_dict[key]["KEGG_Reaction"] = d
        d.clear()
    
    out_f = out_dir / f"{prefix}_annotation.json"
    with open(out_f, "w") as f:
        json.dump(annot_dict, f)
    emapper_dir =  (out_dir / "emapper")
    if not emapper_dir.is_dir():
        emapper_dir.mkdir()
    os.system(f"mv {out_dir}/*emapper.* {emapper_dir}")
    print("The pipeline is finished")
    print(f"output file: {out_f}")
    return out_f

def emapper(fa, output, threads=os.cpu_count(), mode="diamond",
            temp_dir="/dev/shm",
            dmnd_db=None, emapper_path="emapper.py", **kwargs):
    ema_f = Path(f"{output}.emapper.annotations")
    if not ema_f.is_file():
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
        
    with open(ema_f, "r") as f:
        content = [x for x in f.read().rstrip("\n").split("\n")]
        
    result = {}
    for row in content:
        if row.startswith("# "):
            pass
        elif row.startswith("#"):
            cols = row.strip("#").split("\t")
        else:
            row = row.split("\t")
            result[row[0]] = {}
            for col, value in zip(cols, row):
                if col in ("query", "COG_category", "Description",
                           "evalue", "score", "seed_ortholog", "Preferred_name",
                           "max_annot_lvl"):
                    try:
                        value  = float(value)
                    except:
                        value = str(value)
                elif value == "-":
                    value = []
                else:
                    value = value.split(",")
                result[row[0]].update({col: value})
    return result

