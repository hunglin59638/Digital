# Digital: A pipeline to annotate protein sequences  
## Introduction
Digital is a pipeline for annotation of proteins sequence. It could annoate with eggnog-mapper and then alingn NR database with diamond.
## Requirements
+ linux
+ python >= 3.6
+ eggnog-mapper >= v2.1.2
+ diamond

## Installation
Running the `setup.py` will install [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper) and [diamond](https://github.com/bbuchfink/diamond) to `lib` directory.
If the software exists, the installation will be skipped.

```
git clone https://github.com/hunglin59638/digitama.git
cd digitama
./setup.py
./digital -h
```
### Databases
+ diamond  
Downloading NR database by diamond format or you can build by yourself  
paste the link (`https://docs.google.com/u/0/uc?export=download&confirm=aPno&id=1TnXJVYOiJ7vuCxSFdqfTwXKe6ov_jVDX`) to web browser and download file `nr.dmnd.gz`
```
gzip -d nr.dmnd.gz
```
+ eggnog-mapper  
```
 download_eggnog_data.py -y --data_dir [output_dir]
```
## Test
```
digital --fasta test/test.faa --out_dir test/out --prefix test --data_dir [eggnog_data_dir]
```
## Usage
+ use diamond and emapper
```
digital --fasta [protein.fasta] --nr_dmnd [data_dir]/nr.dmnd --out_dir [out_dir]
```
+ skip diamond
  ```
  digital --fasta [protein.fasta] --out_dir [out_dir]
  ```
the results will output to [out_dir]
[out_dir]/output
+ [prefix]_annotation.json  
annotation result
+ diamond.txt  
protein alingment file from diamond  
+ emmaper  
the directory contain files from emapper 

