

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

