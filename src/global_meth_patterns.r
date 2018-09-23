ggplot(RS1.fil, aes(x=pct)) + geom_histogram(bins=20, color="blue", fill="white")

# methylation by chromosome

RS1.fil$chr <- factor(RS1.fil$chr, levels=c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI", "chrXVII", "chrXVIII", "chrXIX", "chrXX", "chrXXI", "chrUn"))

png("../figures/meth_by_chr.png", width=800, height=600)
ggplot(RS1.fil[RS1.fil$chr!="chrM" & !is.na(RS1.fil$chr),], aes(x=chr, y=pct)) +
  geom_boxplot() +
  labs(title="Distribution of Percent Methylation by Chromosome", x="Chromosome Number", y="Percent Methylation") +
  fte_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()