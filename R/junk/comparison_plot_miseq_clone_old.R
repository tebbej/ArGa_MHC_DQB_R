## NGS AND CLONING DATA COMPARISON (HEATMAP) ##
## run clone_analysis.R before running this script, as objects in this script
## are dependant on the objects created there!

miseq_data <- readDNAStringSet("data/ArGa-DQB-NGS_Artemis_20210301.fas")
miseq_data <- miseq_data[-length(miseq_data)]
clone_genotypes <- genotype_info <- readDNAStringSet("data/ArGa-DQB_clone-alleles_20210430.fas")


miseq_in_genotypes <- as.matrix(
  sapply(seq_along(miseq_data), 
         function(x) vcountPattern(as.vector(miseq_data[x]), clone_genotypes)))
colnames(miseq_in_genotypes) <- names(miseq_data)
rownames(miseq_in_genotypes) <- names(clone_genotypes)




miseq_in_genotypes <- miseq_in_genotypes %>% # convert matrix to long data frame for ggplot heatmap format
  as.data.frame() %>%
  rownames_to_column("clone_alleles") %>%
  pivot_longer(-c(clone_alleles), names_to = "ngs_alleles", values_to = "counts") %>%
  mutate(., clone_alleles = factor(clone_alleles, levels = str_sort(unique(clone_alleles), numeric = T)),
         ngs_alleles = factor(ngs_alleles, levels = rev(str_sort(unique(ngs_alleles), numeric = T)))) %>%
  arrange(., desc(ngs_alleles)) %>% #sort both alleles for tidy plotting
  arrange(., clone_alleles) %>% 
  # create column with information whether one or both of the listed allele is a putative artefact or not
  mutate(., artefact_ngs = ifelse((ngs_alleles %in% paste0("ArGa-DQB*", c(15,20,21))) == T, 1, 0)) %>%
  mutate(., artefact_clone = ifelse((clone_alleles %in% paste0("ArGa-DQB*", 20:30)) == T, 1,0))

# color scheme for the axis tick labels
x_col <- c(rep("black", 19), rep("red", 11))
# 1,2,7 index result position from sorting: 21-21 or 20 or 15+1 = 1,2, and 7
y_col <- ifelse((1:21 %in% c(1,2,7)) == T, "red", "black")


ngs_clone_comparison <- ggplot(miseq_in_genotypes, aes(x = clone_alleles,y = ngs_alleles, fill = counts)) + 
  geom_tile(aes(width=0.9, height=0.9)) + 
  scale_fill_viridis_c(option = "cividis") +
  xlab("Putative clone alleles") +
  ylab("Putative MiSeq alleles") +
  theme(panel.background = element_blank(),
        axis.text.x.bottom = element_text(angle = 45, 
                                          vjust = 1, 
                                          hjust = 1, 
                                          size = 11,
                                          color = x_col),
        axis.text.y = element_text(size = 11,
                                   color = y_col),
        plot.margin = unit(c(0.5, 0, 0, 0), "cm"),
        legend.position = "none"
  )
ngs_clone_comparison
# ggsave(filename = "graphics/ngs_clone_comparison.png",
#        plot = ngs_clone_comparison, dpi = 400)
