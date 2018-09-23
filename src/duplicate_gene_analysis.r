RS1.hasdups <- RS1.fil[!is.na(RS1.fil$paralogue_ensembl_gene_id), ] # 10.8% of RS1


# promoters
# ensembl_gene_id, prom_meth, paralogue_ensembl_gene_id, paralogue_prom_meth, pct_id, avg_prom_meth, diff_prom_meth
RS1.fil.prom <- RS1.fil[RS1.fil$region=="promoters", ]
promByGene <- aggregate(pct ~ ensembl_gene_id, RS1.fil.prom, mean)
colnames(promByGene) <- c("paralogue_ensembl_gene_id", "paralogue_prom_meth")
RS1.fil <- join(RS1.fil, promByGene, by="paralogue_ensembl_gene_id")
colnames(promByGene) <- c("ensembl_gene_id", "prom_meth")
RS1.fil <- join(RS1.fil, promByGene, by="ensembl_gene_id")

library(data.table)
setDT(RS1.fil)
promMethDG <- RS1.fil[ , .(prom_meth=unique(prom_meth), paralogue_ensembl_gene_id=unique(paralogue_ensembl_gene_id), paralogue_prom_meth=unique(paralogue_prom_meth), pct_id=unique(avg_pct_id)), by=ensembl_gene_id]
promMethDG <- promMethDG[complete.cases(promMethDG), ]
promMethDG$avg_prom_meth <- rowMeans(as.data.frame(promMethDG)[c("prom_meth", "paralogue_prom_meth")])
promMethDG$diff_prom_meth <- abs(promMethDG$prom_meth - promMethDG$paralogue_prom_meth)

ggplot(promMethDG, aes(x=pct_id, y=diff_prom_meth)) + geom_point() + geom_smooth(method="lm", se=F)



# methylation levels of promoters, gene bodies, and promoters + gene bodies for all genes
wholeGeneMeth <- aggregate(pct ~ ensembl_gene_id, RS1.fil, mean) # n = 17,344
colnames(wholeGeneMeth) <- c("paralogue_ensembl_gene_id", "paralogue_whole_gene_meth")
RS1.fil <- join(RS1.fil, wholeGeneMeth, by="paralogue_ensembl_gene_id")
colnames(wholeGeneMeth) <- c("ensembl_gene_id", "whole_gene_meth")
RS1.fil <- join(RS1.fil, wholeGeneMeth, by="ensembl_gene_id")

library(data.table)
setDT(RS1.fil)
wgMethDGP <- RS1.fil[ , .(whole_gene_meth = unique(whole_gene_meth), paralogue_ensembl_gene_id = unique(paralogue_ensembl_gene_id), paralogue_whole_gene_meth = unique(paralogue_whole_gene_meth), pct_id = unique(avg_pct_id)), by = ensembl_gene_id]
# ensembl_gene_id, whole_gene_meth, paralogue_ensembl_gene_id, paralogue_whole_gene_meth, pct_id
wgMethDGP <- wgMethDGP[complete.cases(wgMethDGP), ]
wgMethDGP$avg_whole_gene_meth <- rowMeans(as.data.frame(wgMethDGP)[c('whole_gene_meth', 'paralogue_whole_gene_meth')])
wgMethDGP$whole_gene_meth_diff <- abs(wgMethDGP$whole_gene_meth - wgMethDGP$paralogue_whole_gene_meth)

ggplot(wgMethDGP, aes(x=pct_id, y=avg_whole_gene_meth)) +
  geom_point(alpha=0.75,size=.75) +
  geom_smooth(method="lm", se=F) +
  labs(title="Percent Identity vs. Gene-Wide Methylation (")
cor(wgMethDGP$pct_id, wgMethDGP$avg_whole_gene_meth)

promoterMethByGene <- aggregate(pct ~ ensembl_gene_id, RS1.fil[RS1.fil$region=="promoters", ], mean) # n = 5,806
colnames(promoterMethByGene) <- c("paralogue_ensembl_gene_id", "paralogue_promoter_meth")
RS1.fil <- join(RS1.fil, promoterMethByGene, by="paralogue_ensembl_gene_id")

geneBodyMeth <- aggregate(pct ~ ensembl_gene_id, RS1.fil[RS1.fil$region %in% c("introns", "exons"), ], mean) # n = 17,030
colnames(geneBodyMeth) <- c("paralogue_ensembl_gene_id", "paralogue_gene_body_meth")
RS1.fil <- join(RS1.fil, geneBodyMeth, by="paralogue_ensembl_gene_id")

# categorize LCA into stickebacks, other fish, and older
RS1.hasdups[RS1.hasdups$last_common_ancestor=="Gasterosteus aculeatus", "LCA_cat"] <- "Sticklebacks"
RS1.hasdups[RS1.hasdups$last_common_ancestor=="Eupercaria" | RS1.hasdups$last_common_ancestor=="Percomorphaceae" | RS1.hasdups$last_common_ancestor=="Acanthomorphata" | RS1.hasdups$last_common_ancestor=="Clupeocephala" | RS1.hasdups$last_common_ancestor=="Neopterygii", "LCA_cat"] <- "Other Fish"
RS1.hasdups[RS1.hasdups$last_common_ancestor %in% c("Euteleostomi", "Vertebrata", "Chordata", "Bilateria"), "LCA_cat"] <- "Older Last Common Ancestor"

# why do these 2 look so diff? sus
ggplot(RS1.hasdups, aes(x=LCA_cat, y=pct)) + geom_boxplot(outlier.size=.25, varwidth = T)

png("../figures/prom_meth_DG_by_LCA.png", width=800, height=600)
ggplot(RS1.hasdups[RS1.hasdups$region=="promoters", ], aes(x=LCA_cat, y=feature_pct_meth)) +
  geom_boxplot(outlier.size=.25, varwidth=T) +
  labs(title="Promoter Methylation by Age (Last Common Ancestor) of Duplicate Genes", x="Last Common Ancestor", y="Average Percent Methylation Across Promoter") +
  fte_theme()
dev.off()

ggplot(RS1.hasdups, aes(x=avg_pct_id, y=pct)) + geom_point()
cor(RS1.hasdups$avg_pct_id, RS1.hasdups$pct)



ggplot(RS1.fil[RS1.fil$region=="promoters" & !is.na(RS1.fil$paralogue_ensembl_gene_id), ], aes(x=avg_pct_id, y=feature_pct_meth)) + geom_point() + geom_smooth()

RS1.fil$last_common_ancestor <- factor(RS1.fil$last_common_ancestor, levels=c("Gasterosteus aculeatus", "Eupercaria", "Percomorphaceae", "Acanthomorphata", "Clupeocephala", "Neopterygii", "Euteleostomi", "Vertebrata", "Chordata", "Bilateria"))



RS1.fil$LCA_cat <- cut(RS1.fil$last_common_ancestor, )

ggplot(RS1.hasdups[RS1.hasdups$region=="promoters", ], aes(x=feature_pct_meth)) + geom_histogram()
ggplot(RS1.fil[RS1.fil$region=="promoters", ], aes(x=feature_pct_meth)) + geom_histogram()
