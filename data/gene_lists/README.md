# Target gene lists

Host factor data were obtained from the `Host_Factors_Summary.xlsx` spreadsheet.
These host factors were the targets in the network flow analysis.
The sheet `Comparison of RNAi Screens, And` was used to create initial gene list files.
The first five columns are screen hits compiled from other researchers.
Mehle lab screens for influenza host factors are provided in the `Our Hits` sheet.

A sorted, concatenated list of the above screening results was then mapped to ENSP, filtered to only include the ENSP identifiers present in the STRINGdb network (experimental edges only), and named `mehle_targets_ensp.txt`.
That file was used as the `--targets-file` in `flow.py`.

`mehle_validation_sources.txt` is a short list of genes whose pathway activities are known.
It was used to verify that the flow method produces networks consistent with the well-studied pathways those genes are members of.
`mehle_validation_sources_ensp.txt` is the above file after mapping from HGNC to ENSP identifiers.
