


# create table with ensembl_gene_id, promoter_methylation, and gene_methylation for each gene in an individual, then average across all individuals

directory_path <- commandArgs(trailingOnly = TRUE)[1]
files <- list.files(path=directory_path, pattern="*.bismark.cov.gz", full.names=TRUE, recursive=FALSE)

lapply(files, function(x) {
  meth_data <- read.table(gzfile(x), col.names=c("chr", "start", "end", "pct", "numCs", "numTs"))
  meth_data$coverage <- meth_data$numCs + meth_data$numTs
  print(summary(meth_data$coverage))
})




# goal: promoter-wide methylation of each gene and gene-wide methylation of each gene for an individual


# load the data
RS1 <- read.table(gzfile("../data_methylation/RS1.trimmed_cutadapt_bismark_bt2.bismark.cov.gz"), col.names = c("chr", "start", "end", "pct", "numCs", "numTs"))

# create coverage column and filter by coverage
RS1$coverage <- RS1$numCs + RS1$numTs
RS1.fil <- RS1[RS1$coverage >= 10, ]
RS1.fil <- RS1.fil[RS1.fil$coverage <= quantile(RS1$coverage, .999), ] # making sure that taking out the lower ones doesn't interfere w/ 99.9th percentile calculation here


library(genomation)
library(methylKit)
library(plyr)


RS1.gr <- as(RS1.fil, "GRanges") # GRanges object needed for findOverlap()
RS1.fil$region <- factor("intergenic_nonpromoter", levels=c("intergenic_nonpromoter", "TSSes", "introns", "exons", "promoters")) # default region; later, sites in other regions will have this column overwritten


## ASSIGN ENSEMBL TRANSCRIPT AND GENE IDs
# this step is somewhat redundant; the feature annotation step adds gene IDs for sites in gene features, but this step should get all the intergenic nonpromoter sites in addition to genic sites

# import transcript IDs and gene IDs table
t2gIDs <- read.csv(gzfile("../data_annotation/t2gIDs.csv.gz"))
colnames(t2gIDs) <- c("ensembl_transcript_id", "ensembl_gene_id")

bed <- read.table("../data_annotation/ensemblGenes.gasAcu1.bed")
bed <- bed[,1:4]
colnames(bed) <- c("seqname", "start", "end", "ensembl_transcript_id")
bed$ensembl_transcript_id <- sub('\\..', '', bed$ensembl_transcript_id)
bed <- join(bed, t2gIDs, by="ensembl_transcript_id")

bed.gr <- as(bed, "GRanges")
geneOverlaps <- as.data.frame(findOverlaps(RS1.gr, bed.gr))

RS1.fil[geneOverlaps$queryHits, "ensembl_gene_id"] <- as.character(bed.gr[geneOverlaps$subjectHits,]$ensembl_gene_id)
RS1.fil[geneOverlaps$queryHits, "ensembl_transcript_id"] <- as.character(bed.gr[geneOverlaps$subjectHits,]$ensembl_transcript_id)


## BED file gives functional regions (promoters, introns, exons, TSSes)

# mean(width(ranges(features$exons[features$exons$score==1,]))) # average length of first exon is 208.49 bp, so we'll include 200bp down.flank in promoters
features <- readTranscriptFeatures("../data_annotation/ensemblGenes.gasAcu1.bed", up.flank=1000, down.flank=200, unique.prom=FALSE)

# order matters; promoters overwrite exons, which overwrite TSSes, etc.
feats <- c("introns", "TSSes", "exons", "promoters")

# create columns to be populated
RS1.fil[ , c("feature_start", "feature_end", "feature_score")] <- NA

for (i in 1:length(feats)) {
  feature.gr <- as(features[[feats[i]]], "GRanges")
  overlaps <- as.data.frame(findOverlaps(RS1.gr, feature.gr))
  
  # ensembl_transcript_id, ensembl_gene_id
  feature.df <- as.data.frame(features[[feats[i]]])
  colnames(feature.df)[7] <- "ensembl_transcript_id"
  feature.df$ensembl_transcript_id <- sub('\\..', '', feature.df$ensembl_transcript_id) # remove version number of transcript
  feature.df <- join(feature.df, t2gIDs, by="ensembl_transcript_id")
  
  RS1.fil[overlaps$queryHits, "region"] <- feats[i]
  
  # feature_start, feature_end, ensembl_transcript_id, ensembl_gene_id
  RS1.fil[overlaps$queryHits, "feature_start"] <- start(ranges(feature.gr[overlaps$subjectHits, ]))
  RS1.fil[overlaps$queryHits, "feature_end"] <- end(ranges(feature.gr[overlaps$subjectHits, ]))
  RS1.fil[overlaps$queryHits, "ensembl_gene_id"] <- as.character(feature.df[overlaps$subjectHits, "ensembl_gene_id"])
  RS1.fil[overlaps$queryHits, "ensembl_transcript_id"] <- feature.df[overlaps$subjectHits, "ensembl_transcript_id"]
  
  # feature_score
  RS1.fil[overlaps$queryHits, "feature_score"] <- feature.gr[overlaps$subjectHits, ]$score
}



prom <- aggregate(pct ~ feature_start + ensembl_gene_id, RS1.fil[RS1.fil$region=="promoters", ], mean)
sum(duplicated(prom$ensembl_gene_id)) # 988
colnames(prom)[3] <- "promoter_methylation"
prom[,1] <- NULL

gene <- aggregate(pct ~ ensembl_gene_id, RS1.fil, mean)
colnames(gene)[2] <- "gene_methylation"

meth <- join(prom, gene, by="ensembl_gene_id")
write.csv(meth, "../data_annotation/RS1.meth.prom.gene.csv")