### Are there differences in methylation across gene features?

library(reshape2)
library(ggplot2)

# methylation by gene feature
png("../figures/meth_by_feat.png", width=800, height=600)
ggplot(RS1.fil, aes(x=region, y=pct)) +
  geom_violin(aes(fill=region), scale="width", show.legend=FALSE) +
  scale_x_discrete(labels=c("intergenic_nonpromoter"="Intergenic", "introns"="Introns", "exons"="Exons", "promoters"="Promoters")) +
  labs(title="Distribution of Percent Methylation by Feature", x="Feature", y="Percent Methylation") +
  fte_theme()
dev.off()

# first vs non-first exon
exon.rows <- RS1.fil[RS1.fil$region=="exons", ]
exon.rows$feature_score_cat <- cut(exon.rows$feature_score, breaks=c(0,1,max(exon.rows$feature_score)), labels=c("First Exon", "Not a First Exon"))

png("../figures/first_vs_nonfirst_exons.png", width=800, height=600)
ggplot(exon.rows, aes(x=feature_pct_meth, fill=feature_score_cat)) +
  geom_density(alpha=.25) +
  labs(title="Distribution of Percent Methylation in First and Non-First Exons", x="Percent Methylation", y="Density") +
  fte_theme() +
  guides(fill=guide_legend(title="")) +
  theme(legend.position=c(0.5, 0.3))
dev.off()

# methylation by exon score
exon.rows.25maxscore <- exon.rows[exon.rows$feature_score <= 25, ]
exon.rows.25maxscore$feature_score <- factor(exon.rows.25maxscore$feature_score)
png("../figures/meth_by_exon_score.png", width=800, height=600)
ggplot(exon.rows.25maxscore, aes(x=feature_score, y=pct)) +
  geom_boxplot(outlier.size=.25) +
  labs(title="Methylation Levels by Exon Score", x="Exon Score", y="Percent Methylation") +
  fte_theme()
dev.off()

# gene methylation by number of exons in gene
geneMethByExonCount <- aggregate(pct ~ ensembl_gene_id, RS1.fil[RS1.fil$region=="exons", ], mean)
geneMethByExonCount <- join(geneMethByExonCount, aggregate(feature_score ~ ensembl_gene_id, RS1.fil[RS1.fil$region=="exons", ], max), by="ensembl_gene_id")
geneMethByExonCount <- geneMethByExonCount[geneMethByExonCount$feature_score <= 25, ]
geneMethByExonCount$feature_score <- factor(geneMethByExonCount$feature_score)
png("../figures/gene_meth_by_exon_count.png", width=800, height=600)
ggplot(geneMethByExonCount, aes(x=feature_score, y=pct)) +
  geom_boxplot(outlier.size=.25) +
  labs(title="Methylation Levels of Genes by Number of Exons", x="Number of Exons in Gene", y="Percent Methylation Across Gene") +
  fte_theme()
dev.off()


#### workin on it
# first vs non-first intron
intron.rows <- RS1.fil[RS1.fil$region=="introns", ]
intron.rows$feature_score_cat <- cut(intron.rows$feature_score, breaks=c(0,1,max(intron.rows$feature_score)), labels=c("First Intron", "Not a First Intron"))

png("../figures/first_vs_nonfirst_introns.png", width=800, height=600)
ggplot(intron.rows, aes(x=feature_pct_meth, fill=feature_score_cat)) +
  geom_density(alpha=.25) +
  labs(title="Distribution of Percent Methylation in First and Non-First Introns", x="Percent Methylation", y="Density") +
  fte_theme() +
  guides(fill=guide_legend(title="")) +
  theme(legend.position=c(0.5, 0.3))
dev.off()





##### old stuff
exon.rows <- RS1.fil[!is.na(RS1.fil$exons_start), ]
# keep those w/ # of sites in at least the 25th percentile (9 sites for exons)
exon.rows.fil <- exon.rows[exon.rows$count_bases_in_exon >= quantile(exon.rows$count_bases_in_exon, .25), ]

intron.rows <- RS1.fil[!is.na(RS1.fil$introns_start), ]
intron.rows.fil <- intron.rows[intron.rows$count_bases_in_intron >= quantile(intron.rows$count_bases_in_intron, .25), ] # 11 sites

## methylation first exons (score=1) vs other exons

summary(exon.rows.fil[exon.rows.fil$exons_score==1, ]$exons_pct_meth)
#   Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.0000   0.4275   1.4829  28.6779  79.6162 100.0000

summary(exon.rows.fil[exon.rows.fil$exons_score!=1, ]$exons_pct_meth)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   80.22   88.29   78.91   92.58  100.00

summary(exon.rows.fil[exon.rows.fil$exons_score==2, ]$exons_pct_meth)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   5.855  83.153  60.012  90.357  99.605

# looks like first exons have lower methylation (possibly more like promoters)
# by the second promoter, the methylation level is similar to the rest of the non-first exons

# plot distributions of exon methylation

# sites must be in an exon & must not be in a promoter, intron, or be a TSS
exonScoreAndMeth <- exon.rows.fil[!is.na(exon.rows.fil$exons_start) & is.na(exon.rows.fil$promoters_start) & is.na(exon.rows.fil$introns_start) & is.na(exon.rows.fil$TSSes), c("exons_score","exons_pct_meth")]

exonScoreAndMeth$exons_score_cat <- cut(exonScoreAndMeth$exons_score, breaks=c(0,1,max(exonScoreAndMeth$exons_score)), labels=c("First Exon", "Not a First Exon"))

png("../figures/first_vs_nonfirst_exons.png", width=800, height=600)
ggplot(exonScoreAndMeth, aes(x=exons_pct_meth, fill=exons_score_cat)) +
  geom_density(alpha=0.25) +
  labs(title="Methylation in First Exons and Non-First Exons", x="Mean Percent Methylation Across Exon", y="Density") +
  guides(fill=guide_legend(title="")) +
  fte_theme() +
  theme(legend.position=c(0.5, 0.3))
dev.off()

## methylation first introns (score=1) vs other intros

summary(intron.rows.fil[intron.rows.fil$introns_score==1, ]$introns_pct_meth)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.000   5.128  51.880  43.401  71.184  98.263

summary(intron.rows.fil[intron.rows.fil$introns_score!=1, ]$introns_pct_meth)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0.00   64.62   76.40   71.42   84.18   99.62

summary(intron.rows.fil[intron.rows.fil$introns_score==2, ]$introns_pct_meth)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0.00   57.18   70.85   64.68   79.35   99.47

# first introns have lower methylation - expecting very high peak near 0
# as with exons, 2nd introns are similar to all the other non-1st introns

# plot distributions of intron methylation

# sites must be in an intron & must not be in a promoter, exon, or be a TSS
exon.rows <- RS1.fil[RS1.fil$region=="exons",]
summary(exon.rows)

exon.rows$feature_score_cat <- cut(exon.rows$feature_score, breaks=c(0,1,max(exon.rows$feature_score)), labels=c("First Exon", "Not First Exon"))
ggplot(exon.rows, aes(x=pct, fill=feature_score_cat)) +
  geom_density(alpha=.25) +
  labs(title="Distribution of Percent Methylation in First and Non-First Exons", x="Percent Methylation in Exonic Sites", y="Density") +
  fte_theme() +
  guides(fill=guide_legend(title="")) +
  theme(legend.position=c(0.5, 0.3))

intron.rows <- RS1.fil[RS1.fil$region=="introns",]
summary(intron.rows$feature_count_sites)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    3.00    6.00   12.34   13.00  582.00
summary(intron.rows$feature_score)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   3.000   6.888   9.000 211.000
summary(intron.rows[intron.rows$feature_score==1,]$feature_pct_meth)
summary(intron.rows[intron.rows$feature_score!=1,]$feature_pct_meth)
intron.rows <- intron.rows[intron.rows$feature_count_sites >= 3,]
intron.rows$feature_score_cat <- cut(intron.rows$feature_score, breaks=c(0,1,max(intron.rows$feature_score)), labels=c("First Intron", "Not First Intron"))
ggplot(intron.rows, aes(x=feature_pct_meth, fill=feature_score_cat)) + geom_density()

intronScoreAndMeth <- intron.rows.fil[!is.na(intron.rows.fil$introns_start) & is.na(intron.rows.fil$promoters_start) & is.na(intron.rows.fil$exons_start) & is.na(intron.rows.fil$TSSes), c("introns_score","introns_pct_meth")]

intronScoreAndMeth$introns_score_cat <- cut(intronScoreAndMeth$introns_score, breaks=c(0,1,max(intronScoreAndMeth$introns_score)), labels=c("First Intron", "Not a First Intron"))

png("../figures/first_vs_nonfirst_introns.png", width=800, height=600)
ggplot(intronScoreAndMeth, aes(x=introns_pct_meth, fill=introns_score_cat)) +
  geom_density(alpha=0.25) +
  labs(title="Methylation in First Introns and Non-First Introns", x="Mean Percent Methylation Across Intron", y="Density") +
  guides(fill=guide_legend(title="")) +
  fte_theme() +
  theme(legend.position=c(0.3, 0.6))
dev.off()

## if a gene has highly methylated exons, does it have a highly methylated introns too?

intronsAndExonsMeth <- join(aggregate(introns_pct_meth ~ensembl_gene_id, RS1.fil, mean), aggregate(exons_pct_meth~ensembl_gene_id, RS1.fil, mean), by="ensembl_gene_id")

#     ensembl_gene_id introns_pct_meth exons_pct_meth
#1 ENSGACG00000000014         85.48885       85.36052
#2 ENSGACG00000000027         93.22687       87.42811
#3 ENSGACG00000000122         85.35713       82.54397
#4 ENSGACG00000000163         90.95238       93.14452
#5 ENSGACG00000000198         88.88889       86.81874
#6 ENSGACG00000000218         84.21053       95.20431

png("../figures/intron_exon_meth_scatter.png", width=800, height=600)
ggplot(intronsAndExonsMeth, aes(x=exons_pct_meth, y=introns_pct_meth)) +
  geom_point() +
  geom_smooth() +
  fte_theme() +
  labs(title="Exon and Intron Methylation", x="Mean Percent Methylation Across Exon", y="Mean Percent Methylation Across Intron")
dev.off()

cor.test(intronsAndExonsMeth$exons_pct_meth, intronsAndExonsMeth$introns_pct_meth)
#t = 35.731, df = 1282, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.6778616 0.7327636
#sample estimates:
#  cor 
#0.7063734

## if a gene has highly methylated exons, does it have a highly methylated introns too?

## -what is the distribution of exon number per gene in the genome? E.g. how many exons does the average gene have?

exonNumberPerGene <- aggregate(score ~ name, as.data.frame(features$exons), max)
summary(exonNumberPerGene$score)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  1.000   4.000   8.000   9.916  13.000 213.000

# -what if you separate the analysis by exon number in a gene and draw boxplots?
# So for genes with only 1 exon, what is the mean %methylation
# For genes with 2 exons, what is the mean %methylation for exon1, intron1, exon2?
# For genes with 3 exons, what is the mean %methylation for exon1, intron1, exon2, intron2, exon3?

RS1.fil.genebody

geneMethByExonScore <- join( aggregate(pct ~ ensembl_gene_id, RS1.fil, mean), aggregate(exons_score ~ ensembl_gene_id, RS1.fil, max), by="ensembl_gene_id" )
geneMethByExonScore.max25exons <- geneMethByExonScore[geneMethByExonScore$exons_score <= 25,]
geneMethByExonScore.max25exons$exons_score <- factor(geneMethByExonScore.max25exons$exons_score)

png("../figures/gene_meth_by_exon_count.png", width=800, height=600)
ggplot(geneMethByExonScore.max25exons, aes(x=exons_score, y=pct)) +
  geom_boxplot() +
  labs(title="Methylation Levels of Genes by Number of Exons", x="Number of Exons in Gene", y="Mean Percent Methylation Across Gene") +
  fte_theme()
dev.off()

# were all those genes protein-coding? noncoding genes tend to only have 1 exon, which would influence methylation level of the first exon

geneIDsTypes <- read.csv(gzfile("../data_annotation/geneIDsTypes.csv.gz"), col.names=c("ensembl_gene_id", "gene_type"))
geneMethByExonScore <- join(geneMethByExonScore, geneIDsTypes, by="ensembl_gene_id")

aggregate(ensembl_gene_id ~ exons_score, geneMethByExonScore, length)
