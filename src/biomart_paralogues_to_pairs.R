library(plyr)
library(biomaRt)
library(data.table)

ens <- useMart("ENSEMBL_MART_ENSEMBL", dataset="gaculeatus_gene_ensembl")

# load paralogue info downloaded from Biomart online
para <- read.csv(gzfile("../data_annotation/ensemblGenes93.stickleback.paralogues.csv.gz"))
colnames(para) <- c("ensembl_gene_id", "paralogue_ensembl_gene_id", "target_query_pct_id", "query_target_pct_id", "paralogue_homology_type", "last_common_ancestor")

para <- para[complete.cases(para), ]
para <- para[para$paralogue_homology_type!="gene_split", ]

# average the 2 percent identity values (1 is target vs query, other is vice versa)
para$avg_pct_id <- (para$target_query_pct_id + para$query_target_pct_id)/2

# get Ensembl protein family IDs corresponding to target gene IDs in paralgoues dataset
famIDsOfTargetGenes <- getBM(attributes=c("ensembl_gene_id", "family"), filters="ensembl_gene_id", values=para$ensembl_gene_id, mart=ens)

# if gene belongs to more than one family, collapse family IDs into a single string per gene
famIDsOfTargetGenes.agg <- aggregate(family ~ ensembl_gene_id, famIDsOfTargetGenes, paste, sep=",")

# get Ensembl protein family IDs corresponding to target gene IDs in paralgoues dataset
famIDsOfQueryGenes <- getBM(attributes=c("ensembl_gene_id", "family"), filters="ensembl_gene_id", values=para$paralogue_ensembl_gene_id, mart=ens)

famIDsOfQueryGenes <- famIDsOfQueryGenes[famIDsOfQueryGenes$family!="",]
famIDsOfQueryGenes.agg <- aggregate(family ~ ensembl_gene_id, famIDsOfQueryGenes, function(x) paste(unique(x)))
colnames(famIDsOfQueryGenes.agg) <- c("paralogue_ensembl_gene_id", "paralogue_family")

# join
para <- join(para, famIDsOfTargetGenes.agg, by="ensembl_gene_id")
para <- join(para, famIDsOfQueryGenes.agg, by="paralogue_ensembl_gene_id")

# convert family and paralogue_family from lists to strings
para$family <- laply(para$family, toString)
para$paralogue_family <- laply(para$paralogue_family, toString)

# keep rows where target family is the same as query family
para_TargetQuerySameFam <- para[para$family==para$paralogue_family, ]

# group by family and select the row with the highest percent identity
dup_pairs <- setDT(para_TargetQuerySameFam)[, .SD[which.max(avg_pct_id)], by=family]
dup_pairs <- dup_pairs[!dup_pairs$family=="",]
write.csv(dup_pairs, "../data_annotation/dup.pairs.csv")