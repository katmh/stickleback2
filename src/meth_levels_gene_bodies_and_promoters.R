# PURPOSE:
# Compare methylation levels across gene bodies, promoters, and other variations (e.g. gene body w/ or w/o introns, whole gene)
# in all genes (with a minimum count of sites), singletons, duplicates (younger vs older), etc.

RS1_ann <- read.csv(gzfile("../data_annotation/RS1_features_dup.csv.gz"))


## EXPLORING HOW COVERED GENES AND GENE FEATURES ARE

summary(RS1_ann$count_sites_in_gene)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 1.0    26.0    46.0    64.5    76.0  1181.0 1464600

png("../figures/sites_per_gene_id.png", width=800, height=600)
ggplot(RS1_ann, aes(x=count_sites_in_gene)) +
  geom_histogram(bins=50, color="blue", fill="white") +
  labs(title="Distribution of Number of Sites per Gene ID in RS1", x="Number of Sites per Gene ID", y="Count") +
  fte_theme()
dev.off()

summary(RS1_ann$count_sites_in_promoter)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 1.0    22.0    38.0    44.2    60.0   195.0 1522909
# Note that this is using the 2kb definition of promoters

png("../figures/sites_per_promoter.png", width=800, height=600)
ggplot(RS1_ann, aes(x=count_sites_in_promoter)) +
  geom_histogram(bins=50, color="blue", fill="white") +
  labs(title="Distribution of Number of Sites per Promoter in RS1", x="Number of Sites per Promoter", y="Count") +
  fte_theme()
dev.off()


## ANALYZE METHYLATION LEVELS IN WHOLE GENES W/ >= 26 SITES
RS1WithGeneIDs <- RS1[!is.na(RS1$ensembl_gene_id), ]

rows.genes.fil <- RS1WithGeneIDs[RS1WithGeneIDs$count_sites_in_gene >= quantile(RS1WithGeneIDs$count_sites_in_gene, .25), ] # 26 sites

prom.rows <- RS1[!is.na(RS1$count_sites_in_promoter), ]

rows.proms.fil <- prom.rows[prom.rows$count_sites_in_promoter >= quantile(prom.rows$count_sites_in_promoter, .25), ]

rows.proms.fil <- RS1WithGeneIDs[RS1WithGeneIDs$count_sites_in_promoter >= quantile(RS1WithGeneIDs$count_sites_in_promoter, .25), ] # 26 sites




genes26Sites <- RS1WithGeneIDs[RS1WithGeneIDs$count_sites_in_gene >= 26, ]
methGenes26Sites <- aggregate(pct ~ ensembl_gene_id, genes26Sites, mean)

png("../figures/whole_gene_meth_all_genes_26sites_RS1.png", width=800, height=600)
ggplot(methGenes26Sites, aes(x=pct)) +
  geom_histogram(bins=50, color="blue", fill="white") +
  labs(title="Distribution of Whole Gene Methylation Level in RS1\n(Genes with >=26 Sites)", x="Mean Percent Methylation Across Whole Gene", y="Count") +
  fte_theme()
dev.off()

## ANALYZE METHYLATION LEVELS IN PROMOTERS W/ >= 22 SITES
prom22Sites <- RS1WithGeneIDs[RS1WithGeneIDs$count_sites_in_promoter >= 22, ]

methPromoters22Sites <- aggregate(pct ~ ensembl_gene_id + promoter_start, prom22Sites, mean)

png("../figures/promoter_meth_22sites_RS1.png", width=800, height=600)
ggplot(methPromotersOfGenes26Sites, aes(x=pct)) +
  geom_histogram(bins=50, color="blue", fill="white") +
  labs(title="Distribution of Promoter Methylation Level in RS1\n(Promoters with >=22 Sites)", x="Mean Percent Methylation Across Promoter", y="Count") +
  fte_theme()
dev.off()




## CORRELATION BTWN PROMOTER METHYLATION AND GENE BODY (INTRONS+EXONS, NO PROMOTER)

RS1WithGeneIDs <- RS1[!is.na(RS1$ensembl_gene_id), ]

genebody.rows <- RS1[!is.na(RS1$ensembl_gene_id) & is.na(RS1$promoter_start), ]

rows.genes.fil <- RS1WithGeneIDs[RS1WithGeneIDs$count_sites_in_gene >= quantile(RS1WithGeneIDs$count_sites_in_gene, .25), ] # 26 sites

prom.rows <- RS1[!is.na(RS1$promoter_start), ]

rows.proms.fil <- prom.rows[prom.rows$count_sites_in_promoter >= quantile(prom.rows$count_sites_in_promoter, .25), ] # 22 sites

promMethPerGene <- aggregate(pct ~ ensembl_gene_id + promoter_start, rows.proms.fil, mean)
colnames(promMethPerGene)[3] <- "prom_pct_meth"

wholeGeneMethPerGene <- aggregate(pct ~ ensembl_gene_id, rows.genes.fil, mean)
colnames(wholeGeneMethPerGene)[2] <- "whole_gene_meth"

promAndWholeGeneMeth <- join(promMethPerGene[,c("ensembl_gene_id", "prom_pct_meth")], wholeGeneMethPerGene, by="ensembl_gene_id")

plot(x=promAndWholeGeneMeth$whole_gene_meth, y=promAndWholeGeneMeth$prom_pct_meth)

ggplot(promAndWholeGeneMeth, aes(x=whole_gene_meth, y=prom_pct_meth)) + geom_point(alpha=0.5) + geom_smooth()
