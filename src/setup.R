setwd("~/stickleback-methylation/src")

# only needed the first time, before dataset annotation
RS1 <- read.table(gzfile("../data_methylation/RS1.trimmed_cutadapt_bismark_bt2.bismark.cov.gz"), col.names = c("chr", "start", "end", "pct", "numCs", "numTs"))
RS1$coverage <- RS1$numCs + RS1$numTs
RS1.fil <- RS1[RS1$coverage >= 10, ]
RS1.fil <- RS1.fil[RS1.fil$coverage <= quantile(RS1$coverage, .999), ] # making sure that taking out the lower ones doesn't interfere w/ 99.9th percentile calculation here

# after performing dataset annotation, just load the file that already has annotations
RS1.fil <- read.csv(gzfile("../data_annotation/RS1_features_dup.csv.gz"))

RS1_obj <- methRead("../data_methylation/RS1.trimmed_cutadapt_bismark_bt2.bismark.cov", sample.id="", assembly="", pipeline="bismarkCoverage", header=FALSE)

library(ggplot2)
library(RColorBrewer)
fte_theme <- function() {
  
  # Generate the colors for the chart procedurally with RColorBrewer
  
  palette <- brewer.pal("Greys", n=9)
  color.background = palette[2]
  color.grid.major = palette[3]
  color.axis.text = palette[6]
  color.axis.title = palette[7]
  color.title = palette[9]
  
  # Begin construction of chart
  
  theme_bw(base_size=9) +
    
    # Set the entire chart region to a light gray color
    
    theme(panel.background=element_rect(fill=color.background, color=color.background)) +
    theme(plot.background=element_rect(fill=color.background, color=color.background)) +
    theme(panel.border=element_rect(color=color.background)) +
    
    # Format the grid
    
    theme(panel.grid.major=element_line(color=color.grid.major,size=.25)) +
    theme(panel.grid.minor=element_blank()) +
    theme(axis.ticks=element_blank()) +
    
    # Format the legend
    theme(legend.background = element_rect(fill=color.background)) +
    theme(legend.text = element_text(size=15,color=color.axis.title)) +
    
    # Set title and axis labels, and format these and tick marks
    
    theme(plot.title=element_text(color=color.title, size=20, vjust=1.25, face="bold")) +
    theme(axis.text.x=element_text(size=15,color=color.axis.text)) +
    theme(axis.text.y=element_text(size=15,color=color.axis.text)) +
    theme(axis.title.x=element_text(size=17,color=color.axis.title, vjust=0, face="bold")) +
    theme(axis.title.y=element_text(size=17,color=color.axis.title, vjust=1.25, face="bold")) +
    
    # Plot margins
    
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
}
