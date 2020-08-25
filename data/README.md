# Data sources

`Top_50_PB2_ICC-MS_interactors.xlsx` contains the top 50 hits of PB2 interactors by affinity purification mass spectrometry.
The hits are ranked by confidence combined from two different summed datasets.

A subset of the top 50 hits that have the desired curve shape were mapped to HGNC identifiers and then ENSP names.
Each of these nodes was checked for membership in the STRINGdb network (experimental edges only) and included in the file `hf_curve_shape_ensp_stringdb_filter.txt` if it existed in the network.
This file is used as the `--sources-file` for `flow.py`.

The [`gene_lists`](gene_lists) subdirectory contains the host factor screening data used to create targets for the network flow analysis.
