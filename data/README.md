# Data sources

`Top_50_PB2_ICC-MS_interactors.xlsx` contains the top 50 hits of PB2 interactors by affinity purification mass spectrometry.
The hits are ranked by confidence combined from two different summed datasets.

A subset of the top 50 hits that have the desired curve shape were mapped to HGNC identifiers and then ENSP names.
Each of these nodes was checked for membership in the STRINGdb network (experimental edges only) and included in the file `hf_curve_shape_ensp_stringdb_filter.txt` if it existed in the network.
This file is used as the `--sources-file` for `flow.py`.

The [`gene_lists`](gene_lists) subdirectory contains the host factor screening data used to create targets for the network flow analysis.

The `--edges-file` can be obtained at https://doi.org/10.6084/m9.figshare.21588270 and is a transformation of STRINGdb data into the ABC file format (nodeA, nodeB, weightC).

The `--mapping-file` can be obtained at https://doi.org/10.6084/m9.figshare.21588285 and is a snapshot of data from biomaRt, mapping HNGC to ENSP.

`Meistermann_2014_FlaviVirus_sources.txt` contains flavivirus sources from https://doi.org/10.1074/mcp.M113.028977. `Marceau_2016_FlaviVirus_targets.txt` contains targets from https://www.nature.com/articles/nature18631.
