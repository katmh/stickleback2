df <- data.frame(
  organism=c("A. thaliana (flowering plant)", "M. musculus (mouse)", "E. coli (bacterium)", "D. melanogaster (fruit fly)", "D. discoideum (slime mold)", "H. sapiens (human)", "G. aculeatus (stickleback)"),
  percent_Cs_in_CpGs_methylated=c(24,7.6,2.3,.0034,.006,70,60.3)
)

library(ggplot2)

png("../figures/cpg_meth_diff_species.png", width=800, height=600)
ggplot(data=df, aes(x=organism, y=percent_Cs_in_CpGs_methylated, fill=organism)) +
  geom_bar(stat="identity", color="blue") +
  geom_text(aes(label=percent_Cs_in_CpGs_methylated), hjust = -0.15, size=6.5) +
  coord_flip() +
  labs(title="CpG Methylation Across Species", y="% of Cs in CpG Contexts that are Methylated", x="Species") +
  fte_theme()
dev.off()
