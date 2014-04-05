# pretender: A RNA-Seq simulator for biological replicates

## Important note

This project is currently a **work in progress**. I'm updating it fairly
regularly and expect to have it working in the next few days. It's basically
not even alpha.

## Description

Pretender is a simulator of RNA-Seq reads. It is very idealized and attempts to
simulate *reads* with the distributional assumptions that many methods make
(i.e. negative binomial for biological replicates).

## Getting gene names

```R
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset
      = "hsapiens_gene_ensembl")
mart <- useMart(biomart = "ensembl", dataset
      = "mmusculus_gene_ensembl")
listAttributes(mart)
geneNames <- getBM(attributes
      = c("ensembl_transcript_id","external_gene_id"), mart
      = mart)
write.csv(geneNames, file = "geneNames.csv")
```

