#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 19 09:32:07 2020

@author: hunglin
"""
import subprocess
import os
from pathlib import Path
root_dir = Path(__file__).parent
print(f"root_dir: {root_dir.resolve()}")
lib_path = (root_dir / "lib")
lib_path.mkdir(exist_ok=True)
print("Setting eggNog...\n")
if not subprocess.call(["which","emapper.py"]):
    print("Already install eggnog-emapper\n")
else:
    os.chdir(lib_path)
    os.system("git clone https://github.com/jhcepas/eggnog-mapper.git")
    os.chdir("eggnog-mapper")
    os.system("python3 setup.py install")
    # os.system("python download_eggnog_data.py -y")
    
print("Setting diamond...\n")
if not subprocess.call(["which","diamond"]):
    print("Already install diamond\n")
else:
    lib_path.mkdir(exist_ok=True)
    os.chdir(lib_path)
    os.system("wget -qO- http://github.com/bbuchfink/diamond/releases/download/v2.0.8/diamond-linux64.tar.gz | tar xzf -")


print("Setting python packages\n")
require_f = root_dir / "requirements.txt"
if require_f.is_file():
    os.chdir(root_dir)
    if subprocess.call(["pip", "install", "-r", require_f]):
        print("fail to install packages in requirements.txt")
else:
    print(f"{require_f} is not found")
    for pkg in ("biopython", "tqdm", "termcolor", "fake_useragent","setuptools",
                "requests", "urllib3", "pandas"):
        os.system(f"pip install {pkg}")
print("Finished !")