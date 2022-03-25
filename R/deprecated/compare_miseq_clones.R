## NGS AND CLONING DATA COMPARISON (HEATMAP?) ##
miseq_data <- readDNAStringSet("data/ArGa-DQB-NGS_Artemis_20210301.fas")
miseq_data <- miseq_data[-length(miseq_data)]
miseq_info <- read.table(file = "data/ngs/Acacia_data_artemis.txt", header = T)[1:21,]
clone_genotypes <- readDNAStringSet("data/ArGa-DQB_clone-alleles_20210430.fas") #order by occurence frequency in cloning seqs
# in case allele ordering is desired to be after actual genotypes:
# clone_genotypes <- readDNAStringSet("data/allele_order_from_genotypes.fas")


#assign frequency values from plot to correct alleles (unsorted ngs alleles)
freq_info <- data.frame(ngs_alleles = names(miseq_data),
                        freq = miseq_info$freq[1:21]/100)
                   # freq =   c(.15,.08,.10,.07,.11,.08,.07,.05,.06,.04,.02,
                   #            .03,.02,.03,.01,.02,.02,.02,.00,.01,.01))

# along sequences in miseq_data, compare the sequences in miseq_data with 
# sequences in clone_genotypes and create a matrix with 1's indicating correspondence
# of sequences and 0's non-correspondence
miseq_in_genotypes <- as.array(sapply(seq_along(miseq_data), 
         function(x) vcountPattern(as.vector(miseq_data[x]), clone_genotypes)))
# adjust col and row names according to their corresponding data.frame of origin
colnames(miseq_in_genotypes) <- names(miseq_data)
rownames(miseq_in_genotypes) <- names(clone_genotypes)

miseq_in_genotypes <- miseq_in_genotypes %>% # convert matrix to long data frame for ggplot heatmap format
  as.data.frame() %>%
  rownames_to_column("clone_alleles") %>%
  pivot_longer(-c(clone_alleles), names_to = "ngs_alleles", values_to = "counts") %>%
  mutate(., clone_alleles = factor(clone_alleles, levels = str_sort(unique(clone_alleles), numeric = T)),
         ngs_alleles = factor(ngs_alleles, levels = rev(str_sort(unique(ngs_alleles), numeric = T)))) %>%
  arrange(., desc(ngs_alleles)) %>% #sort both alleles for tidy plotting
  arrange(., clone_alleles)  %>%
  # create column with information whether one or both of the listed allele is a putative artefact or not
  mutate(., artefact_ngs = ifelse((ngs_alleles %in% paste0("ArGa-DQB*", c(15,20,21))) == T, 1, 0)) %>%
  mutate(., artefact_clone = ifelse((clone_alleles %in% paste0("ArGa-DQB*", 20:30)) == T, 1,0))

shared_alleles <- miseq_in_genotypes[-which(miseq_in_genotypes$counts != 1),]

shared_alleles <- miseq_in_genotypes[-which(miseq_in_genotypes$counts != 1),1:2] %>%
  mutate(sequence = as.vector(
    miseq_data[na.omit(
      match(shared_alleles$ngs_alleles, names(miseq_data)))]))

# double check whether sequences now coincide
iter_bin <- vector(length=dim(shared_alleles)[1])
for (i in seq_along(shared_alleles[[3]])) {
  iter_bin[i] <- str_detect(as.character(clone_genotypes[i]), shared_alleles[[3]][[i]])
}
if (any(iter_bin != T)) {
  warning("Allele matches share erroneous sequences")
} else {
  message("Allele sequences coincide")
}

# update frequency of occurences of MiSeq alleles and overwrite names in clone allele manner
# (miseq alleles, are named according to their matching clone allele sequence)
miseq_alleles <- shared_alleles %>% 
  mutate(
    freq = miseq_info$freq[na.omit(match(shared_alleles$ngs_alleles, miseq_info$allele))]/100,
    variant_count = miseq_info$count[na.omit(match(shared_alleles$ngs_alleles, miseq_info$allele))]
    ) %>%
  mutate(ngs_alleles = clone_alleles) %>%
  select(!(clone_alleles))

miseq_alleles <- as.data.frame(miseq_alleles)

miseq_alleles_variant_counter <- unlist(sapply(miseq_alleles$variant_count, function(x) 1:x))

miseq_alleles_expand <- as.data.frame(lapply(miseq_alleles, rep, miseq_alleles$variant_count)) %>%
  mutate(., variant_counter = miseq_alleles_variant_counter) 

miseq_alleles_expand <- miseq_alleles_expand %>%
  mutate(mutate(.,variant_no  = sapply(
    miseq_alleles_expand$ngs_alleles, 
    function(x) {
      stringr::str_split(x, "\\*")[[1]][2] %>% 
        paste0(., collapse = "\\*") %>% 
        as.numeric()
    })))


figures[[7]] <- ggplot(miseq_alleles_expand, 
       aes(x = variant_no, 
           group = desc(variant_counter), 
           fill = variant_counter)) + 
  geom_bar(aes(y = stat(count) / sum(count))) + 
  scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                     limits = c(0,0.2), #swapped x & y values for reversing!
                     expand = c(0,0)) +
  scale_fill_viridis_c(option = "cividis",
                       begin = 0,
                       end = 1) +
  ylab("Frequency\n'MiSeq'\n") +
  labs(fill = "Allele\ncounts") + 
  scale_x_continuous(breaks = seq_along(unique(miseq_alleles$ngs_alleles)),
                     labels = str_sort(unique(miseq_alleles$ngs_alleles), numeric = T),
                     expand = c(0, 0.3)) +  
  theme_minimal() +
  theme(panel.grid = element_line(color = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black", 
                                  margin = margin(10,10,20,10)),
        axis.ticks = element_line(color = "black", 
                                  size = 0.2),
        axis.line.x = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x.bottom = element_text(angle = 45, 
                                          vjust = 1, 
                                          hjust = 1, 
                                          size = 11),
        axis.line.y = element_line(color = "black"),
        axis.title.y = element_text(size = 15.5),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(.15,"cm"),
        plot.background = element_rect(color = "white", 
                                       fill = "white"),
        legend.position = "none",
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")
  )
figures[[7]]

ggsave(filename = "graphics/allele_frequency_ngs_artemis_preliminary_scale.png",
       figures[[7]],
       dpi = 400)

freq_compare <- data.frame(clone_freq = alleles$frequency[1:19]/100,
                           ngs_freq =miseq_alleles$freq)

ngs_clone_correlation <- ggplot(freq_compare, aes(x=ngs_freq, y=clone_freq)) +
  geom_point(shape = 15,
             size = 3) +
  geom_smooth(method = "lm", 
              se = T,
              colour = "orange") + 
  scale_y_continuous(name   = "Allele frequency\ncloning data\n",
                     limits = c(0,0.2),
                     expand = c(0,0)) +
  scale_x_continuous(name   = "\nAllele frequency\nMiSeq data",
                     limits = c(0,0.2),
                     expand = c(0,0)) +
  annotate(geom = "text", x = 0.181, y = 0.01,
           label = "italic(R²)==0.84", parse = T) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.margin = unit(c(0.5,0.5,0.5,0,5), "cm"),
    axis.title = element_text(color = "black", 
                              margin = margin(10,10,20,10))
  )
ngs_clone_correlation
ggsave("graphics/freq_correlation.png",
       ngs_clone_correlation,
       dpi = 400)

corr <- lm(ngs_freq~clone_freq, data = freq_compare)
summary(corr)
