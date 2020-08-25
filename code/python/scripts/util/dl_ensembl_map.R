#!/usr/bin/env Rscript
library(biomaRt)

main = function() {
  usage = "usage: dl_ensembl_map.R <hgnc_complete_set> <ensembl_map>
hgnc_complete_set is a mapping of HGNC identifiers to Ensembl gene identifiers provided by EBI at ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt
ensembl_map is the output file to write the hgnc_symbol, ensembl_gene_id, ensembl_transcript_id, ensembl_peptide_id associations to
"
  args = commandArgs(TRUE)
  if(length(args) != 2) {
    write(usage, stderr())
    quit(status = 1)
  }
  hgnc_complete_set_fp = args[1]
  ensembl_map_fp = args[2]

  # extract ensemble gene ids from hgnc_complete_set
  hgnc_complete_set_df = read.table(hgnc_complete_set_fp, sep='\t', header=TRUE, quote="", fill=TRUE)
  hgnc_ensembl_df = hgnc_complete_set_df[,c("symbol","ensembl_gene_id")]
  ensembl_ids_proto = hgnc_complete_set_df[,"ensembl_gene_id"]
  ensembl_gene_ids = Filter(function(x) { x != "" }, ensembl_ids_proto)

  # query biomaRt
  mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  results <- getBM(
    attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id"),
    filters = "ensembl_gene_id", values=ensembl_gene_ids,
    mart = mart
  )

  # append hgnc identifiers to results
  rv_df = merge(hgnc_ensembl_df, results)
  rv_df = rv_df[c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id", "symbol")] # reorder columns
  write.table(rv_df, file=ensembl_map_fp, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
}

main()
