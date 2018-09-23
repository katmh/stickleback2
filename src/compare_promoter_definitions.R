# purpose: different studies have used different definitions of promoters
#
# here I'll compare:
# 1) promoters defined as 200bp upstraeam of TSS
# 2) promoters defined as 1kb upstream and 1kb downstream of TSS
#
# by looking at:
# 1) distribution (histogram) of promoter-wide methylation (dist_prom_meth_2kb/200bp_RSx)
# 2) bar plot of promoter-wide methylation by LCA of duplicate pairs (prom_meth_by_LCA_2kb/200bp_RSx)
#
# this file is only a code template for making those 2 plots. you need to change the promoter definition in dataset_annotation.R and run that file

###  PLOT #1
png("../figures/dens_prom_meth_2kb_RS1.png", width=800, height=600)
ggplot(RS1, aes(x=promoter_pct_meth)) +
  geom_density(color="blue", fill="white") +
  fte_theme() +
  labs(title="Distribution of Mean Promoter-Wide Percent Methylation in Sample RS1\n(Promoter = 1kb upstream + 1kb downstream of TSS)", x="Mean Percent Methylation Across Promoter", y="Density")
dev.off()

### PLOT #2
prom_by_LCA$last_common_ancestor <- factor(prom_by_LCA$last_common_ancestor, levels=c("Gasterosteus aculeatus", "Eupercaria", "Percomorphaceae", "Acanthomorphata", "Clupeocephala", "Neopterygii"))

prom_by_LCA <- aggregate(promoter_pct_meth ~ last_common_ancestor, RS1, mean)
n_prom_by_LCA <- aggregate(promoter_pct_meth ~ last_common_ancestor, RS1, length)
colnames(n_prom_by_LCA) <- c("last_common_ancestor", "count_proms")
prom_by_LCA <- join(prom_by_LCA, n_prom_by_LCA)

df <- aggregate(start~ensembl_gene_id,RS1,length)

png("../figures/prom_meth_by_LCA_2kb_RS1.png", width=800, height=600)
ggplot(prom_by_LCA, aes(x=last_common_ancestor, y=promoter_pct_meth)) +
  geom_bar(stat="identity", color="blue", fill="white") +
  coord_flip() +
  geom_text(aes(label=paste("n = ", count_proms, sep="")), hjust = -0.15, size=6.5) +
  fte_theme() +
  labs(title="Mean Promoter-Wide Methylation in Duplicate Genes\n(Promoter = 1kb upstream + 1kb downstream of TSS)", x="Mean Percent Methylation Across Promoter", y="Last Common Ancestor") +
  scale_y_continuous(limits=c(0,60))
dev.off()
