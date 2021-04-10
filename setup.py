#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 19 09:32:07 2020

@author: hunglin
"""
import subprocess
import os
from pathlib import Path
lib_path = (Path(__file__).parent / "lib")
if not lib_path.is_dir():
    lib_path.mkdir()
print("Setting eggNog...\n")
if not subprocess.call(["which","emapper.py"]):
    print("Already install eggnog-emapper\n")
else:
    os.chdir(lib_path)
    os.system("git clone https://github.com/jhcepas/eggnog-mapper.git")
    os.chdir("eggnog-mapper")
    os.system("./download_eggnog_data.py -y")
    
print("Setting diamond...\n")
if not subprocess.call(["which","diamond"]):
    print("Already install diamond\n")
else:
    os.chdir(lib_path)
    os.system("wget -qO- http://github.com/bbuchfink/diamond/releases/download/v2.0.8/diamond-linux64.tar.gz | tar xzf -")


print("Setting python packages\n")
requirements = Path(__file__).parent / "requirements.txt"
if subprocess.call(["pip", "install", "-r", requirements]):
    raise Exception("fail to install packages in requirements.txt")
    
print("Finished !")