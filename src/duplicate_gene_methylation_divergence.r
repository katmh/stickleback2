# p = promoters
# goal: get table with these columns: ensembl_gene_id, paralogue_ensembl_gene_id, pct_id, prom_meth, paralogue_prom_meth, avg_prom_meth, diff_prom_meth

RS1.fil.p <- RS1.fil[RS1.fil$region=="promoters", ] # 11.77% of RS1.fil

# each promoter is identified by its start position and gene ID; add metadata: pct_meth and count_sites
pMethByGene <- aggregate(feature_pct_meth ~ feature_start + ensembl_gene_id, RS1.fil.p, unique)
pMethByGene <- join(pMethByGene, aggregate(feature_count_sites ~ feature_start + ensembl_gene_id, RS1.fil.p, unique), by=c("feature_start", "ensembl_gene_id"))

# some gene IDs have multiple transcripts -> multiple TSSes -> multiple promoters
# choose the promoter w/ the highest amount of sites that we have methylation data for
library(dplyr)
pMethByGene <- pMethByGene %>% group_by(ensembl_gene_id) %>% top_n(1, feature_count_sites)

# join w/ methylation data
colnames(pMethByGene)[2:4] <- c("paralogue_ensembl_gene_id", "paralogue_prom_meth", "paralogue_prom_count_sites")
RS1.fil.p <- join(RS1.fil.p, pMethByGene[,2:4], by="paralogue_ensembl_gene_id")

write.csv(pMethByGene, "../data_annotation/promoter_meth_by_gene.csv")

# avg_prom_meth, diff_prom_meth

RS1.fil.p <- RS1.fil.p[complete.cases(RS1.fil.p), ]
head(rowMeans(RS1.fil.p[, c("feature_pct_meth","paralogue_prom_meth")]))

prom <- as.data.table(RS1.fil.p)[ , .(paralogue_ensembl_gene_id=unique(paralogue_ensembl_gene_id), pct_id=unique(avg_pct_id), prom_meth=unique(feature_pct_meth), paralogue_prom_meth=unique(paralogue_prom_meth), last_common_ancestor=unique(last_common_ancestor)), by=ensembl_gene_id]
prom <- as.data.frame(prom)
prom$avg_prom_meth <- rowSums(prom[c("prom_meth", "paralogue_prom_meth")])
prom$diff_prom_meth <- abs(prom$prom_meth - prom$paralogue_prom_meth)

# LCA order
prom$last_common_ancestor <- factor(prom$last_common_ancestor, levels=c("Gasterosteus aculeatus", "Eupercaria", "Percomorphaceae", "Acanthomorphata", "Clupeocephala", "Neopterygii", "Euteleostomi", "Vertebrata", "Chordata", "Bilateria"))

# categorize LCA into stickebacks, other fish, and older
prom[prom$last_common_ancestor=="Gasterosteus aculeatus", "LCA_cat"] <- "Sticklebacks"
prom[prom$last_common_ancestor=="Eupercaria" | prom$last_common_ancestor=="Percomorphaceae" | prom$last_common_ancestor=="Acanthomorphata" | prom$last_common_ancestor=="Clupeocephala" | prom$last_common_ancestor=="Neopterygii", "LCA_cat"] <- "Other Fish"
prom[prom$last_common_ancestor %in% c("Euteleostomi", "Vertebrata", "Chordata", "Bilateria"), "LCA_cat"] <- "Older Last Common Ancestor"




##### PROMOTER DIVERGENCE PLOTS

png("../figures/promoter_meth_div_by_LCA.png", width=800, height=600)
ggplot(prom, aes(x=diff_prom_meth, col=LCA_cat)) + geom_density(size=.75) +
  fte_theme() +
  theme(legend.position = c(0.8, 0.5), legend.title=element_text(size=14)) +
  labs(title="Promoter Methylation Divergence by Last Common Ancestor Groups", x="Difference in Promoter Methylation Between Duplicate Partners", y="Density", col="Last Common Ancestor Group")
dev.off()

png("../figures/promoter_meth_avg_by_LCA.png", width=800, height=600)
ggplot(prom, aes(x=avg_prom_meth, col=LCA_cat)) + geom_density(size=.75) +
  fte_theme() +
  theme(legend.position = c(0.8, 0.5), legend.title=element_text(size=14)) +
  labs(title="Average Promoter Methylation by Last Common Ancestor Groups", x="Average Promoter Methylation Between Duplicate Partners", y="Density", col="Last Common Ancestor Group")
dev.off()



# plots
png("../figures/promoter_meth_div_dup_genes.png", width=800, height=600)
ggplot(prom, aes(x=pct_id, y=diff_prom_meth)) +
  geom_point() + 
  geom_smooth(method="lm", se=F) +
  labs(title="Percent Identity vs. Promoter Methylation Divergence in Duplicate Gene Pairs", x="Percent Identity Between Duplicate Partners", y="Difference in Promoter Methylation of Duplicate Partners") +
  fte_theme()
dev.off()

png("../figures/promoter_meth_dup_genes.png", width=800, height=600)
ggplot(prom, aes(x=pct_id, y=avg_prom_meth)) +
  geom_point() + 
  geom_smooth(method="lm", se=F) +
  labs(title="Percent Identity vs. Promoter Methylation in Duplicate Gene Pairs", x="Percent Identity Between Duplicate Partners", y="Average Promoter Methylation of Duplicate Partners") +
  fte_theme()
dev.off()

png("../figures/LCA_vs_pct_ID.png", width=800, height=600)
ggplot(prom, aes(x=last_common_ancestor, y=pct_id)) +
  geom_boxplot() +
  labs(title="Last Common Ancestor vs. Percent Identity of Duplicate Gene Pairs", x="Last Common Ancestor Between Duplicate Partners", y="Percent Identity Between Duplicate Partners") +
  fte_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



### GENE BODY METHYLATION

# format the data
RS1.fil.gb <- RS1.fil[RS1.fil$region==c("introns", "exons"), ]
nrow(RS1.fil.gb)/nrow(RS1.fil) # 22.85%
methByGeneID <- aggregate(pct~ensembl_gene_id, RS1.fil.gb, mean)
colnames(methByGeneID)[2] <- "gb_meth"

dgGbMeth <- as.data.table(RS1.fil.gb)[ , .(paralogue_ensembl_gene_id=unique(paralogue_ensembl_gene_id), pct_id=unique(avg_pct_id), last_common_ancestor=unique(last_common_ancestor)), by=ensembl_gene_id]
dgGbMeth <- join(dgGbMeth, methByGeneID, by="ensembl_gene_id")

colnames(methByGeneID) <- c("paralogue_ensembl_gene_id", "paralogue_gb_meth")
dgGbMeth <- join(dgGbMeth, methByGeneID, by="paralogue_ensembl_gene_id")

dgGbMeth <- dgGbMeth[complete.cases(dgGbMeth), ] # 2605 pairs
dgGbMeth$avg_gb_meth <- rowMeans(as.data.frame(dgGbMeth)[c("gb_meth", "paralogue_gb_meth")])
dgGbMeth$diff_gb_meth <- abs(dgGbMeth$gb_meth - dgGbMeth$paralogue_gb_meth)

# LCA order
dgGbMeth$last_common_ancestor <- factor(dgGbMeth$last_common_ancestor, levels=c("Gasterosteus aculeatus", "Eupercaria", "Percomorphaceae", "Acanthomorphata", "Clupeocephala", "Neopterygii", "Euteleostomi", "Vertebrata", "Chordata", "Bilateria"))

# categorize LCA into stickebacks, other fish, and older
dgGbMeth[dgGbMeth$last_common_ancestor=="Gasterosteus aculeatus", "LCA_cat"] <- "Sticklebacks"
dgGbMeth[dgGbMeth$last_common_ancestor=="Eupercaria" | dgGbMeth$last_common_ancestor=="Percomorphaceae" | dgGbMeth$last_common_ancestor=="Acanthomorphata" | dgGbMeth$last_common_ancestor=="Clupeocephala" | dgGbMeth$last_common_ancestor=="Neopterygii", "LCA_cat"] <- "Other Fish"
dgGbMeth[dgGbMeth$last_common_ancestor %in% c("Euteleostomi", "Vertebrata", "Chordata", "Bilateria"), "LCA_cat"] <- "Older Last Common Ancestor"

# gene body divergence between duplicate partners, colored by LCA group
png("../figures/gb_meth_div_by_LCA.png", width=800, height=600)
ggplot(dgGbMeth, aes(x=diff_gb_meth, col=LCA_cat)) +
  geom_density(size=.75) +
  fte_theme() +
  theme(legend.position = c(0.8, 0.5), legend.title=element_text(size=14)) +
  labs(title="Gene Body Methylation Divergence by Last Common Ancestor Groups", x="Difference in Gene Body Methylation Between Duplicate Partners", y="Density", col="Last Common Ancestor Group")
dev.off()

# average gene body methylation between duplicate partners, colored by LCA group
png("../figures/avg_gb_meth_by_LCA.png", width=800, height=600)
ggplot(dgGbMeth, aes(x=avg_gb_meth, col=LCA_cat)) +
  geom_density(size=.75) +
  fte_theme() +
  theme(legend.position = c(0.8, 0.5), legend.title=element_text(size=14)) +
  labs(title="Average Gene Body Methylation by Last Common Ancestor Groups", x="Average Gene Body Methylation Between Duplicate Partners", y="Density", col="Last Common Ancestor Group")
dev.off()

#### age (percent identity) vs gene body methylation
png("../figures/gb_meth_dup_genes.png", width=800, height=600)
ggplot(dgGbMeth, aes(x=pct_id, y=avg_gb_meth)) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  labs(title="Gene Body Methylation vs. Percent Identity in Duplicate Genes", x="Percent Identity Between Duplicate Partners", y="Average Gene Body Methylation of Duplicate Partners") +
  fte_theme()
dev.off()

png("../figures/gb_meth_div_dup_genes.png", width=800, height=600)
ggplot(dgGbMeth, aes(x=pct_id, y=diff_gb_meth)) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  labs(title="Gene Body Methylation Divergence vs. Percent Identity in Duplicate Genes", x="Percent Identity Between Duplicate Partners", y="Difference in Gene Body Methylation of Duplicate Partners") +
  fte_theme()
dev.off()

png("../figures/gb_meth_div_dup_genes.png", width=800, height=600)
ggplot(dgGbMeth, aes(x=pct_id, y=diff_gb_meth)) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  labs(title="Gene Body Methylation Divergence vs. Percent Identity in Duplicate Genes", x="Percent Identity Between Duplicate Partners", y="Difference in Gene Body Methylation of Duplicate Partners") +
  fte_theme()
dev.off()

methByGeneID <- join(methByGeneID, RS1.fil.gb[c("ensembl_gene_id", "paralogue_ensembl_gene_id", "last_common_ancestor", "avg_pct_id")], by="ensembl_gene_id")
