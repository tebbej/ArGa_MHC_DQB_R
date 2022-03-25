gtypes_ae <- read.table(file = "data/ngs_genotypes_AE.txt", 
                        sep = "\t", 
                        skip = 1,
                        col.names = c("id", "allele", "colony", "paired"))
str(gtypes_ae)

keeep <- c("W8584pup",	"ArGa-DQB*11_var2")

ii <- match(gtypes_ae$allele, shared_alleles$ngs_alleles)
gtypes_ae <- gtypes_ae[-which(is.na(ii) == T),] %>%
  mutate(., id = as.factor(id))

gtypes_ngs <- data.frame(
  id = as.factor(gtypes_ae$id),
  allele = shared_alleles$clone_alleles[na.omit(ii)],
  colony = as.factor(gtypes_ae$colony),
  paired = gtypes_ae$paired
) %>%
  mutate(., id = factor(id, levels = rev(levels(id)))) %>%
  arrange(., desc(id))

gtypes_ngs <- gtypes_ngs %>%
  mutate(., allele = factor(gtypes_ngs$allele, levels = levels(allele)[1:19])) %>%
  mutate(., allele_var = sapply(
    as.character(gtypes_ngs$allele), 
    function(x) {
      stringr::str_split(x, "\\*")[[1]][2] %>%                                  # "\\" double escape regex * in order for R to recognize it as char
        paste0(., collapse = "\\*") %>%
        as.numeric()
    }))

# include only paired individuals
gtypes_ngs <- gtypes_ngs[which(gtypes_ngs$paired == 1),]

# plot heatmap for ssb
ngs_heat_ssb <- ggplot(
  data = gtypes_ngs[which(gtypes_ngs$colony == "SSB"),],
  aes(x = allele, y = id)) + 
  geom_tile(fill = "#FDE725") + 
  coord_fixed(ratio = 0.6) +
  xlab("Alleles") +
  ylab("Sample ID") +
  scale_x_discrete(drop = F) +
  # labs(fill = "Log \nclone \nnumber") +
  # theme_minimal() +
  theme(
    panel.background   = element_rect(fill = "#440154"),
    panel.grid = element_line(color = "#440154"),
    axis.text.x.bottom = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_line(),
    plot.margin = unit(c(0.5, 0.5, 0, 0), "cm"),
    legend.position = "none"
  )
ngs_heat_ssb

# plot heatmap for fwb
ngs_heat_fwb <- ggplot(
  data = gtypes_ngs[which(gtypes_ngs$colony == "FWB"),],
  aes(x = allele, y = id)) + 
  geom_tile(fill = "#FDE725") + 
  coord_fixed(ratio = 0.6) +
  xlab("Alleles") +
  ylab("Sample ID") +
  scale_x_discrete(drop = F) +
  # labs(fill = "Log \nclone \nnumber") +
  # theme_minimal() +
  theme(
    panel.background   = element_rect(fill = "#440154"),
    panel.grid = element_line(color = "#440154"),
    axis.text.x.bottom = element_text(angle = 45, 
                                      vjust = 1, 
                                      hjust = 1, 
                                      size = 11),
    plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"),
    legend.position = "none"
  )
ngs_heat_fwb

ngs_heat <- egg::ggarrange(ngs_heat_ssb,
                               ngs_heat_fwb,
                               ncol=1)

ngs_heat <- ggpubr::as_ggplot(ngs_heat)



