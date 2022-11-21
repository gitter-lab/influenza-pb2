# FloPro

FloPro is a @flo@w @pro@gram. This repository is a supplement to our manuscript:

[Alternative splicing liberates a cryptic cytoplasmic isoform of mitochondrial MECR that antagonizes influenza virus](https://doi.org/10.1101/2020.11.09.355982)
Steven F Baker, Helene Meistermann, Manuel Tzouros, Aaron Baker, Sabrina Golling, Juliane Siebourg Polster, Mitchell P Ledwith, Anthony Gitter, Angelique Augustin, Hassan Javanbakht, Andrew Mehle.
_bioRxiv_ 2020 doi:10.1101/2020.11.09.355982

## Contents
- [`code`](code): a Python library and scripts for network flow and enrichment analysis
- [`condor`](condor): scripts for running the analyses in a HTCondor cluster
- [`data`](data): source and target gene lists for the network flow analysis

## Install python dependencies
```
# setup virtual environment for flopro dependencies
python -m venv env
source env/bin/activate

# install dependencies
pip install -r requirements.txt

# install flopro
cd code/python
python setup.py install
```

## Download data dependencies
```
cd data
wget -O edges_file.txt https://figshare.com/ndownloader/files/38264646
wget -O mapping_file.txt https://figshare.com/ndownloader/files/38264649
```

## Run flopro
```
mkdir flow_results
flow.py --min-sources 1 --min-targets 1 --edges-file data/edges_file.txt --mapping-file data/mapping_file.txt --sources-file data/hf_curve_shape_ensp_stringdb_filter.txt --targets-file data/gene_lists/mehle_targets_ensp.txt --outdir flow_results
```
