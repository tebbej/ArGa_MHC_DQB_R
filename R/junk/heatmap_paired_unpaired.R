# create separate heatmaps for mum-pup pairs and single individuals

# create vector with info about mum pup pairings
paired <- c(rep(0,24), rep(1,7), 0, 0, rep(1,4), 0, rep(1,14), 0,1,1,0) %>%
  factor(., levels = c(0,1))

# bind to df and arrange for paired
allele_summary <- cbind(allele_summary, paired = rep(paired,30)) %>%
  arrange(., desc(paired)) %>%
  arrange(., alleles) %>%
  mutate(., sample_id = fct_reorder(sample_id, as.numeric(paired)))

allele_summary_paired <- allele_summary[allele_summary$paired == 1,] %>%
  arrange(., sample_id)

allele_summary_unpaired <- allele_summary[allele_summary$paired == 0,] %>%
  arrange(., sample_id) %>%
  mutate(., sample_id = fct_reorder(sample_id, as.numeric(paired)))

# plot heatmaps
clone_heatmap_unpaired <- ggplot(allele_summary_unpaired, aes(x = alleles,y = sample_id, fill = log(counts+1))) + 
  geom_tile() + 
  coord_fixed(ratio = 0.6) +
  scale_fill_viridis_c() +
  xlab("Alleles") +
  ylab("Sample ID") +
  labs(fill = "Log \nclone \nnumber") +
  theme(
    axis.text.x.bottom = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = unit(c(0.5, 0, 0, 0), "cm"),
    legend.position = "none"
  )
clone_heatmap_unpaired

clone_heatmap_paired <- ggplot(allele_summary_paired, aes(x = alleles,y = sample_id, fill = log(counts+1))) + 
  geom_tile() + 
  coord_fixed(ratio = 0.6) +
  scale_fill_viridis_c() +
  xlab("Alleles") +
  ylab("Sample ID") +
  labs(fill = "Log \nclone \nnumber") +
  theme(
    axis.text.x.bottom = element_text(angle = 45, 
                                      vjust = .98, 
                                      hjust = .94),
    axis.title.y = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )
clone_heatmap_paired


clone_heatmap <- egg::ggarrange(clone_heatmap_unpaired,
                                clone_heatmap_paired,
                                ncol=1)

clone_heatmap <- ggpubr::as_ggplot(clone_heatmap)




































# plot heatmap for ssb
# figures[[2]] <- 
  ggplot(
  data = allele_summaryX,
  aes(x = alleles,
      y = sample_id, 
      fill = log(counts+1))) + 
  geom_tile() + 
    # facet_grid(vars(family)) +
  coord_fixed(ratio = 0.6) +
  scale_fill_viridis_c() +
  xlab("Alleles") +
  # ylab("SSB\n") +
  labs(fill = "Log \nclone \nnumber") +
  theme_minimal() +
  theme(
    axis.text.x.bottom = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_line(color = "white",
                               size = 1),
    # axis.title.y = element_blank(),
    axis.ticks.y = element_line(),
    # axis.text.y  = element_text(color = color_ssb),
    plot.margin  = unit(c(0.5, 0, 0, 0), "cm"),
    legend.position = "right"
  )
names(figures)[2] <- "clone_heatmap_ssb"
