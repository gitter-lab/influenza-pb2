# FloPro
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7342881.svg)](https://doi.org/10.5281/zenodo.7342881)

FloPro is a **flo**w **pro**gram. This repository is a supplement to our manuscript:

[Alternative splicing liberates a cryptic cytoplasmic isoform of mitochondrial MECR that antagonizes influenza virus](https://doi.org/10.1371/journal.pbio.3001934).
Steven F Baker, Helene Meistermann, Manuel Tzouros, Aaron Baker, Sabrina Golling, Juliane Siebourg Polster, Mitchell P Ledwith, Anthony Gitter, Angelique Augustin, Hassan Javanbakht, Andrew Mehle.
_PLoS Biology_, 20:12, 2022.

## Contents
- [`code`](code): a Python library and scripts for network flow and enrichment analysis
- [`condor`](condor): scripts for running the analyses in a HTCondor cluster
- [`data`](data): source and target gene lists for the network flow analysis

## System dependencies
- Currently only `Linux` is officially supported and this package requires `graphviz` for some of its features, particularly the visualization of sets of enriched genes.

## Install python dependencies
We have tested `Python >= 3.8`, but other Python 3's may work as well.

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
cd ../../data
wget -O edges_file.txt https://figshare.com/ndownloader/files/38264646
wget -O mapping_file.txt https://figshare.com/ndownloader/files/38264649
```

## Run flopro
```
# back to repo root
cd ../
mkdir flow_results
python code/python/scripts/flow.py --min-sources 1 --min-targets 1 --edges-file data/edges_file.txt --mapping-file data/mapping_file.txt --sources-file data/hf_curve_shape_ensp_stringdb_filter.txt --targets-file data/gene_lists/mehle_targets_ensp.txt --outdir flow_results
```
