# PURPOSE:
# Annotate .bismark.cov data with information on:
# - Functional regions (promoters, introns, exons, intergenic), with no overlap
# - Ensembl gene and transcript IDs, and count of CpG sites in gene
# - Gene ontology term accession and names
# - Paralogues (Ensembl gene ID of paralogue, last common ancestor, paralogue homology type)

library(genomation)
library(methylKit)
library(plyr)

## first, import data in setup.r and filter by coverage

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

mean(width(ranges(features$exons[features$exons$score==1,]))) # average length of first exon is 208.49 bp, so we'll include first exons in promoters
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

featureInfo <- aggregate(pct ~ feature_start + region + chr, RS1.fil, mean)
colnames(featureInfo)[4] <- "feature_pct_meth"
featureInfo <- join(featureInfo, aggregate(pct ~ feature_start + region + chr, RS1.fil, length), by=c("feature_start", "region", "chr"))
colnames(featureInfo)[5] <- "feature_count_sites"
RS1.fil <- join(RS1.fil, featureInfo, by=c("feature_start", "region", "chr"))

## ADD GENE ONTOLOGY ANNOTATIONS

GO <- read.delim("../data_annotation/GO.txt.gz", header=TRUE) # import gene ontology data
GO <- GO[!GO$GO.term.name=="", ] # remove rows w/ no GO.term.name; note that there are some rows that have accessions but lack a term name
GO$Transcript.stable.ID <- NULL # remove Ensembl transcript IDs

GO$Gene.stable.ID <- as.character(GO$Gene.stable.ID)
GO$GO.term.accession <- as.character(GO$GO.term.accession)
GO$GO.term.name <- as.character(GO$GO.term.name)

GO.agg.term <- aggregate(GO.term.name ~ Gene.stable.ID, GO, paste, collapse=";") # put all GO term names in one row
# length(unique(GO$Gene.stable.ID)) is 15,561, and so is nrow(GO.agg.term)
GO.agg.acc <- aggregate(GO.term.accession ~ Gene.stable.ID, GO, paste, collapse=";")
GO.agg <- join(GO.agg.acc, GO.agg.term, by="Gene.stable.ID")
colnames(GO.agg) <- c("ensembl_gene_id", "GO_term_accession", "GO_term_name")

RS1 <- join(RS1, GO.agg, by="ensembl_gene_id") # join methylation and GO data by matching Ensembl gene ID
# sum(!is.na(RS1$GO_term_accession)) and sum(!is.na(RS1$GO_term_name)) both equal 339,566

## DUPLICATE GENE DATA

dup_pairs <- read.csv(gzfile("../data_annotation/dup.pairs.csv"))
dup_pairs[ , c("X", "paralogue_homology_type", "target_query_pct_id", "query_target_pct_id", "paralogue_family")] <- NULL
RS1.fil <- join(RS1.fil, dup_pairs, by="ensembl_gene_id")

## SAVE DATAFRAME W/ DATA ABOUT METHYLATION, REGIONS, DUPLICATES, AND GO TERMS

write.csv(RS1.fil, "../data_annotation/RS1_meth_prom_dup.csv")

inGenicRegion <- RS1.fil[!is.na(RS1.fil$promoters_start) | !is.na(RS1.fil$exons_start) | !is.na(RS1.fil$introns_start) | !is.na(RS1.fil$TSSes),] # n=1,016,529
#> nrow(inGenicRegion[is.na(inGenicRegion$ensembl_gene_id),])
#[1] 597,862

# > annotateWithFeature(RS1.gr, features$promoters)
# summary of target set annotation with feature annotation:
# Rows in target set: 3865721
# ----------------------------
# percentage of target elements overlapping with features:
# features$promoters              other 
#               2.55              97.45 
#
# percentage of feature elements overlapping with target:
# [1] 38.29

# > annotateWithGeneParts(RS1.gr, features)
# Summary of target set annotation with genic parts
# Rows in target set: 3865721
# -----------------------
# percentage of target features overlapping with annotation:
# promoter       exon     intron intergenic 
#     2.55      23.32      29.16      46.27 
#
# percentage of target features overlapping with annotation:
# (with promoter > exon > intron precedence):
# promoter       exon     intron intergenic 
#     2.55      23.05      28.13      46.27 
#
# percentage of annotation boundaries with feature overlap:
# promoter     exon   intron 
#    38.29    33.53    36.26 
#
# summary of distances to the nearest TSS:
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0    1896    5907   13852   16007  381038