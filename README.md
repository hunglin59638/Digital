# Digitama: pipeline to annotate protein sequences  
## Introduction
Digitama is a pipeline for annotation of proteins sequence. It could annoate with eggnog-mapper and then alingn NR database with diamond.
## Requirements
+ linux
+ python >= 3.6

## Installation
Running the `setup.py` will install [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper) and [diamond](https://github.com/bbuchfink/diamond) to `lib` directory.
If the software exists, the installation will be skipped.

```
git clone https://github.com/hunglin59638/digitama.git
cd digitama
python3 setup.py
digitama -h
```
Downloading NR database by diamond format or you can make by yourself
paste the link (`https://docs.google.com/u/0/uc?export=download&confirm=aPno&id=1TnXJVYOiJ7vuCxSFdqfTwXKe6ov_jVDX`) to browser and download file `nr.dmnd.gz`
```
gzip -d nr.dmnd.gz
```
## Usage
```
digitama --fasta [protein.fasta] --nr_dmnd [data_dir]/nr.dmnd --out_dir [out_dir]
```
the results will output to [out_dir]
[out_dir]/output
+ [prefix]_annotation.json  
annotation result
+ ko2name.json  
ko number in annotation result to its definition
+ go_result.json  
GO ids in annotation result to its complete imformation  
