para <- read.csv(gzfile("../data_annotation/stickleback_paralogues_Ensembl93.csv.gz"))
para <- para[complete.cases(para), ]

# average difference between the 2 percent identity metrics

para$diff <- abs(para$Paralogue..id..target.Stickleback.gene.identical.to.query.gene - para$Paralogue..id..query.gene.identical.to.target.Stickleback.gene)

png("../figures/diff_btwn_percent_ID_paralogues.png", width=800, height=600)
ggplot(para, aes(x=diff)) + geom_histogram(bins=50, color="blue", fill="white") + labs(title="Distribution of Absolute Differences Between Percent Identity Metrics", x="|(% ID of Target to Query) - (% ID of Query to Target)|", y="Count") + fte_theme()
dev.off()

# the distributions of the 2 percent identity metrics are identical: summary
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00   35.73   88.80   70.31   98.32  100.00

png("../figures/dist_pct_ID.png", width=800, height=600)
ggplot(para, aes(x=Paralogue..id..target.Stickleback.gene.identical.to.query.gene)) + geom_density(color="blue", fill="white") + fte_theme() + labs(title="Distribution of Percent Identity of Target to Query Genes in Sticklebacks", x="Percent Identity of Gene (Target) to its Paralogue (Query)", y="Density")
dev.off()