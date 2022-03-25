require(Biostrings)
require(stringr)
require(tidyverse)

#read in dataset with genotypes (2 Alleles per individual even if homozygous)
test_genotypes <- readDNAStringSet(filepath = 'data/ArGa_MHC_clone_sample_genotypes_20210114.fas') 
                                             #"data/ArGa_MHC_test_sample_genotypes_truncate_20201110.fas"

allele_genotypes <- as.vector(test_genotypes)
allele_genotypes <- as.factor(allele_genotypes)

alleles <- levels(allele_genotypes)
allele_number <- length(alleles)

test <- allele_genotypes %in% alleles

# substitute with a MiSeq allele check where each DQB Allele is compared with each unique allele found in this data analysis
sample_genotype <- read.csv(file = "data/ArGa_MHC_DQB_genotypes_20201111.csv", sep = ";")
head(sample_genotype)

allele_count_total <- c(as.vector(sample_genotype$Allele1), as.vector(sample_genotype$Allele2))

unique_alleles <- unique(allele_count_total)
# sort allele names by number of first occurence in MiSeq
unique_alleles <- str_sort(unique_alleles, numeric = T)

which(allele_count_total %in% allele_count_total[1])

allele_freq <- 1:length(unique_alleles)
allele_names <- 1:length(unique_alleles)
for (i in 1:length(unique_alleles)) {
  allele_freq[i] <- length(which(allele_count_total %in% unique_alleles[i]))
  allele_names[i] <- unique_alleles[i]
}
names(allele_freq) <- allele_names
allele_freq
sum(allele_freq)

allele_freq <- (allele_freq/40)*100
sum(allele_freq)

allele_freq.df <- data.frame(allele_names, allele_freq, pointer = 1:length(allele_freq))

allele_freq_plot <- ggplot(data = allele_freq.df, aes(x = pointer, y = allele_freq)) +
  geom_col(color = "black", fill = "turquoise") +
  ylab("Frequency [%]") +
  # xlab("Allele") +
  scale_x_continuous(breaks = allele_freq.df$pointer,
                   labels = allele_freq.df$allele_names) +
  ylim(c(0,20)) +
  theme_minimal() +
  theme(panel.grid = element_line(color = "darkgrey"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black", margin = margin(10,10,20,10)),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x.bottom = element_text(angle = 45, vjust = 1.3, hjust = 1.05),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(.15,"cm"),
        plot.background = element_rect(color = "white", fill = "white"),
        legend.position = "none",
        plot.margin = unit(c(0.5,0.5,0.5,0,5), "cm")
  )  
allele_freq_plot
ggsave("graphics/allele_frequency_cloning_data.png", allele_freq_plot, dpi = 400)



