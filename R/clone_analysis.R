# rm(list = ls())
require(Biostrings) # easily work with genetic string sets
require(tidyverse) # package collection for easy and pretty data science with R
require(gridExtra) # ggplot grid manipulations
require(egg) # ggplot grid and plot alignment functions
require(ggpubr) # ggplot grid and plot alignment functions
require(ape)
require(poppr)

## BUILDING ANALYSIS DATA SETS ##
# Data set to analyse: Fasta file that contains multiple clones for one individual
genotype_info <- readDNAStringSet("data/Clones_MHC_ArGa_exon_20210319.fas")

# Data set containing alleles identified by Illumina NGS analysis
ngs_alleles <- readDNAStringSet("data/ArGa-DQB-NGS_Artemis_20210301.fas") %>%
  as.vector()

# Data set containing meta data about the sampled Antarctic fur seals
metadata_df <- read.table(file   = "data/sample_list.txt", 
                          header = T) %>%
  mutate(
    real_id  = factor(real_id, 
                     levels = str_sort(real_id, 
                                       numeric = T)),
    colony   = as.factor(colony),
    maturity = as.factor(maturity),
    family   = as.factor(family)
    ) %>%
  arrange(real_id) %>%
  arrange(colony)

# Subset of metadata_df containing only paired individuals (mother-pup pairs)
metadata_df_pairs <- metadata_df[
  which(
    duplicated(metadata_df$family) | 
    duplicated(metadata_df$family, fromLast = T) 
    == T)
  ,]


# identify index position and name of uniques seqs in the data set
get_unique_seqs <- function(raw_clones, allowed.unique.seq = 1, max.mismatch = 0){
  
  if (is.vector(raw_clones) != T){
    raw_clones <- as.vector(raw_clones)
  } 
  # sort all instances that only occur once
  unique_indicator <- vector(length = length(raw_clones))
  for (i in seq_along(raw_clones)) {
    # unique_indicator[i] <- sum(match(raw_clones, raw_clones[i]), na.rm = T)
    unique_indicator[i] <- sum(vcountPattern(as.character(raw_clones[i]), raw_clones,
                                             max.mismatch = max.mismatch, min.mismatch = 0)
                               , na.rm = T)    
  }
  
  index_unique_seq <- which(unique_indicator == allowed.unique.seq | unique_indicator <= allowed.unique.seq)
  names_unique_seq <- names(raw_clones[index_unique_seq])
  
  return(list(index_unique_seq = index_unique_seq, 
              names_unique_seq = names_unique_seq))
}#end get_unique_seqs

get_allele_info <- function(data, allele_name = "ArGa-DQB*", rm.unique = T){
  if (rm.unique == T) {
    unique_identifier <- get_unique_seqs(data, 
                                         allowed.unique.seq = 1, 
                                         max.mismatch = 0)
    
    # delete unique sequences from the data
    if (!purrr::is_empty(unique_identifier[[1]])) {
      genotype_fas <- data[-unique_identifier$index_unique_seq]
    }  
  } else {
    genotype_fas <- data
  }
  
  
  # create data set with allele sequences
  allele_seq <- unique(genotype_fas)
  # console output for no. of alleles that can be identified
  message(paste0("Data set contains ", length(unique(genotype_fas)), " unique sequences/alleles"))
  
  # count occurences of identified alleles in the data
  allele_count <- as.vector(
    sapply(seq_along(allele_seq), function(i) 
      sum(match(genotype_fas, allele_seq[[i]]), na.rm = T))
  )
  
  # create data.frame where alleles will be named after its decreasing frequency in the data
  # sorted by the prior allele count
  alleles <- data.frame(seq = as.vector(allele_seq), 
                        counts = allele_count,
                        row.names = NULL) %>%
    arrange(., desc(counts)) %>%
    # find shared sequence between cloning data and ngs data
    mutate(., shared_ngs = apply(
      # go through each sequence in ngs_alleles
      sapply(ngs_alleles, 
             # apply vcountPattern to find where an ngs_alleles sequence is found in alleles$seq
             vcountPattern, 
             subject = seq),  
      # result is a matrix with rows x col = length of alleles$seq x ngs_alleles
      # apply sums each of those rows (share the same indexing as alleles$seq) to identify where
      # sequences are shared between cloning and ngs data sets
      1, sum)) %>% 
    
    mutate(frequency = (counts/sum(counts))*100) %>%
    `rownames<-`(., sapply(seq_along(allele_count), function(i) 
      paste0(allele_name,i))) %>%
    rownames_to_column("name")
  
  out = list(alleles = alleles, 
             genotype_fas = genotype_fas)
  
  return(out)
} # end get_allele_info

out <- get_allele_info(genotype_info)

alleles <- out[[1]]
genotype_fas <- out[[2]]

# create data.frame where each cloning sequence has the info on its corresponding allele number
# allele numbering follows the exact same naming convention as allele names (1 most frequent, 30 least frequent)

clone_allele_df <- as.data.frame(genotype_fas) %>%
  transmute(sequence = x) %>%
  rownames_to_column(var = "clone_var")

allele_index_in_df <- as.vector(
  sapply(
    clone_allele_df$sequence, 
    function(x) match(x, alleles$seq)))

clone_allele_df <- {clone_allele_df %>%
  transmute(.,
            id         = sapply(
              clone_allele_df$clone_var, 
              function(x) {
                stringr::str_split(x, "-")[[1]][1] %>% 
                  paste0(., collapse = "-") %>% 
                  as.factor()
              }),
            clone_var  = clone_var,
            allele     = alleles$name[allele_index_in_df],
            variant_no = allele_index_in_df,
            variant_count  = alleles$counts[allele_index_in_df],
            sequence   = sequence
  ) %>%
  mutate(., allele = factor(
    allele, 
    levels = str_sort(
      unique(allele), 
      numeric = T))
  ) %>%
  arrange(
    ., allele
  ) %>% 
  mutate(
    variant_counter = as.vector(
      unlist(
        sapply(alleles$counts, 
               function(x) seq(1:x))
      ))
  ) %>%
  relocate(., sequence, .after = last_col())
}


#update clone_allele_df with metadata information
# create index vector where sample ids correspond to the correct names in the metadata data.frame
index <- match(as.character(clone_allele_df$id), metadata_df$clone_id)

# rearrange columns in clone_allele_df based on 'index'
clone_allele_df <- clone_allele_df %>%
  mutate(
    id       = metadata_df$real_id[index],
    colony   = metadata_df$colony[index],
    maturity = metadata_df$maturity[index],
    family   = metadata_df$family[index]
  ) %>%
  relocate(., sequence, .after = last_col())

clone_allele_df <- clone_allele_df %>% arrange(., variant_no)

# create suitable data frame for a heatmap that contains allele names, sample ids and the respective number an allele occurs in
# sample id
allele_summary <- matrix(nrow = length(unique(clone_allele_df$id)),
                         ncol = length(unique(clone_allele_df$allele))) %>%
  `rownames<-`(., as.character(unique(clone_allele_df$id))) %>%
  `colnames<-`(., as.character(
    str_sort(
      levels(
        clone_allele_df$allele), 
      numeric = T)))

# fill matrix with info on which and how many alleles are found in the clones for each individual fur seal
for (i in seq_along(unique(clone_allele_df$id))) {
  alleles_in_id <- summary(clone_allele_df$allele[clone_allele_df$id == unique(clone_allele_df$id)[i]])
  allele_summary[i, ] <- alleles_in_id[str_sort(names(alleles_in_id), numeric = T)]
}

# convert to data.frame and create a "tidy" version, ggplot and tidyverse can handle easily
allele_summary <- allele_summary %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column("alleles") %>%
  pivot_longer(-c(alleles), 
               names_to = "sample_id", 
               values_to = "counts") %>%
  mutate(., alleles = factor(
    alleles, 
    levels = str_sort(
      unique(alleles), 
      numeric = T)),
    sample_id = as.factor(sample_id)) %>%
  arrange(., sample_id) %>%
  arrange(., alleles)

index2 <- match(as.character(allele_summary$sample_id), metadata_df$real_id)
allele_summary <- allele_summary %>%
  mutate(colony = metadata_df$colony[index2],
         maturity = metadata_df$maturity[index2],
         family = metadata_df$family[index2],
         sample_id = factor(sample_id, 
                            levels = rev(
                              levels(sample_id)))) %>%
  arrange(desc(sample_id)) %>%
  arrange(alleles)

pair_match_index <- match(allele_summary$sample_id, metadata_df_pairs$real_id)
pair_match_index <- which(is.na(pair_match_index) == T)
allele_summaryX <- allele_summary[-pair_match_index,]

allele_summaryX <- allele_summaryX[allele_summaryX$alleles %in% levels(allele_summaryX$alleles)[1:19],] %>%
  # create color palette of 5 colors for mum and corresponding pup
  # 20 families, so each color iterates 4 times in a plot for a pair)
  mutate(color_code = rep(c(rep("#f58231", 2), 
                            rep("#f032e6", 2), 
                            rep("#000000", 2),
                            rep("#bfef45", 2), 
                            rep("#4363d8", 2)), 
                          # length of allele_summaryX after trimming divided by 5 = 152 divided by 2 = 76
                          76))

## Create data.frame with genotype information
clone_genotype_df <- allele_summary[
  allele_summary$alleles %in% levels(allele_summary$alleles)[1:19],]

clone_genotype_df <- clone_genotype_df[which(clone_genotype_df$counts != 0),]

f1 <- function(x){
  length(
    na.omit(
      match(clone_genotype_df$alleles, x)
    ))
}

clone_genotype_df <- clone_genotype_df %>%
  mutate(., 
         variant_no = clone_allele_df$variant_no[
           match(clone_genotype_df$alleles, 
                 clone_allele_df$allele)],
         freq = clone_allele_df$allele_frequency[
           match(clone_genotype_df$alleles, 
                 clone_allele_df$allele)],
  ) %>%
  arrange(., alleles) %>%
  mutate(., 
         variant_counts = unlist(
           sapply(clone_genotype_df$alleles, f1))) %>%
  mutate(.,
         variant_counter = unlist(
           sapply(
             sapply(
               unique(
                 clone_genotype_df$alleles), f1), 
             function(x) seq(1:x)))) %>%
  arrange(., desc(variant_counter))


#### PLOTS AND FIGURES RAW ####
# store figures
figures <- vector(mode = "list")

# frequency plot with clone sequences
figures[[1]] <- ggplot(clone_allele_df[1:771,], # indexing only variants that occur more than once
                       aes(x = variant_no, 
                           group = dplyr::desc(variant_counter), 
                           fill = variant_counter)) + 
  geom_bar(aes(y = stat(count) / sum(count))) + 
  scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                     limits = c(0,.2),
                     expand = c(0,0)) +
  scale_fill_viridis_c(option = "viridis",
                       begin = 0,
                       end = 1) +
  ylab("Frequency\n") +
  labs(fill = "No. of\nclones") + 
  scale_x_continuous(name = "MHC DQB II allele",
                     breaks = seq_along(unique(clone_allele_df$allele)),
                     labels = str_sort(unique(clone_allele_df$allele), numeric = T),
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
        axis.title.x = element_text(color = "black", size = 15.5),
        axis.text.x.bottom = element_text(angle = 45, 
                                          vjust = 1, 
                                          hjust = 1, 
                                          size = 13),
        axis.line.y = element_line(color = "black"),
        axis.title.y = element_text(size = 15.5),
        axis.ticks.y = element_blank(),
        axis.text.y.left = element_text(size = 13),
        axis.ticks.length = unit(.15,"cm"),
        plot.background = element_rect(color = "white", 
                                       fill = "white"),
        legend.position = c(0.88,.8),
        legend.background = element_rect(fill = "white",
                                         color = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,0,5), "cm")
  )
names(figures)[1] <- "allele_freq_hist"

ggsave(file = "graphics/allele_frequency_cloning_data_no_artefacts.png", 
       figures$allele_freq_hist, dpi = 400)

figures[[1]]

# plot heatmap for all
# extract color pallettes
color_ssb <- with(allele_summaryX, color_code[which(colony == "SSB")])
color_fwb <- with(allele_summaryX, color_code[which(colony == "FWB")])

sep_heat <- function(x, color = "blue"){ 
  df <- allele_summaryX[which(allele_summaryX$family == x),]
  gg <- ggplot( data = df,
                aes(x = alleles,
                    y = sample_id, 
                    fill = log(counts+1))) + 
    geom_tile() + 
    coord_fixed(ratio = 0.6) +
    scale_fill_viridis_c(limits = c(0, 3.5),
                         begin = 0.05, 
                         breaks = 0:3,
                         labels = c(0, 2.7, 7.4, 20.1)) +
    # xlab("Alleles") +
    # ylab("SSB\n") +
    labs(fill = "Log \nclone \nnumber") +
    # theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.x.bottom = element_blank(),
      axis.text.y = element_text(size = 10,
                                      hjust = 0),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.y.left = element_line(color = color,
                                 size = 1),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      # axis.text.y  = element_text(color = color_ssb),
      plot.margin  = unit(c(-.1, 0, 0, 0), "cm"),
      legend.position = "none"
    )
  
  return(gg)
  
}

vec <- unique(allele_summaryX$family)
vec_ssb <- unique(allele_summaryX$family[allele_summaryX$colony == "SSB"])
vec_fwb <- unique(allele_summaryX$family[allele_summaryX$colony == "FWB"])

plots_ssb <- lapply(vec_ssb, sep_heat, color = "white")
plots_fwb <- lapply(vec_fwb, sep_heat, color = "white")

plot_list <- c(plots_ssb, plots_fwb)

# (heat_all <- ggarrange(plotlist = plot_list, ncol = 1, align = "v"))
# # ggsave("graphics/test_heat.png", heat_all, dpi = 400, width = 6.5, height = 6.5, units = "in")
# 
# egg::ggarrange(plots = plot_list, ncol = 1)

# plot heatmap for last individual with x axis
l <- length(vec)
plot_list[[l]] <- ggplot(
  data = allele_summaryX[which(allele_summaryX$family == vec[l]),],
  aes(x = alleles,
      y = sample_id, 
      fill = log(counts+1))) + 
  geom_tile() + 
  coord_fixed(ratio = 0.6) +
  scale_fill_viridis_c(name = "No. of\nclones",
                       limits = c(0,3.5),
                       begin = 0.05, 
                       breaks = 0:3,
                       labels = c(0, 2.7, 7.4, 20.1)) +
  xlab("Alleles") +
  ylab("FWB\n") +
  # labs(fill = "Log \nclone \nnumber") +
  # theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_line(color = "#000000"),
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(size = 10,
                               hjust = 0), 
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y.left = element_line(color = "white", 
                                    size = 1),
    axis.text.x.bottom = element_text(angle = 45,
                                      vjust = 1,
                                      hjust = 1,
                                      size = 10),
    axis.title.x = element_text(size = 14),
    plot.margin = unit(c(-1, 0, 0, 0), "mm"),
    legend.position = "right",
    legend.background = element_rect(fill = "white",
                                     color = "black")
  )

# increase headspace for first list element
plot_list[[1]] <- plot_list[[1]] + theme(plot.margin = unit(c(10,0,0,0), "mm"))

figures[[4]] <- egg::ggarrange(plots = plot_list, ncol = 1)

names(figures)[4] <- "clone_heatmap"

ggsave(filename = "graphics/clone_heat_map_pairs_only_-artefacts.svg",
       figures[[4]],
       width = 254.6,
       height = 180,
       units = "mm",
       dpi = 400)

figures[[4]]


# ## PHYLOGENETIC TREE ##
# # read in tree that was created in MEGAX (.newick file format)
# 
# seal_dqb <- read.tree(file = "data/phylogeny/seal_dqb.nwk")
# 
# node_name <- read.table(file = "data/phylogeny/node_name.txt", sep = "\t", header = T) %>%
#   as.data.frame() %>%
#   mutate(., grp = as.factor(grp))
# 
# seal_dqb$tip.label <- node_name$id
# 
# species.exp <- c(
#   expression(italic("Arctocephalus forsteri")), 
#   expression(italic("Arctocephalus gazella")), 
#   expression(italic("Halicheorus grypus")), 
#   expression(italic("Homo sapiens")), 
#   expression(italic("Mirounga angustirostris")), 
#   expression(italic("Mirounga leonina")), 
#   expression(italic("Monachus schauinslandi")), 
#   expression(italic("Neophoca cinerea")), 
#   expression(italic("Odobenus rosmarus")), 
#   expression(italic("Phocarctos hookeri")), 
#   expression(italic("Zalophus californianus")), 
#   expression(italic("Zalophus wollebaeki"))             
#                  )
# 
# legend.palette <- RColorBrewer::brewer.pal(12,name = "Paired")
# legend.palette[11] <- "#000000"
#   
# 
# p <- ggtree(seal_dqb, branch.length='none', layout='circular')
# p <- p %<+% node_name + geom_tippoint(aes(color = grp), size = 3) + 
#   geom_tiplab(hjust = -.08,
#               size = 5) +
#   xlim(0,28) +
#   scale_color_manual("Species", values = legend.palette,
#                      labels = species.exp) + 
#   theme(text = element_text(size = 18), 
#         legend.position = "right",
#         legend.background = element_rect(size = .5, color = "black"),
#         legend.justification = "top",
#         legend.title = element_blank()
#   ) + 
#   guides(color = guide_legend(ncol = 1, 
#                               label.hjust = 0,
#                               title.position = "bottom"))
# 
# figures[[5]] <- p
# names(figures)[5] <- "seal_dqb_phylogeny"
# rm(p)
# 
# figures[[5]]
# 
# ggsave(filename = "graphics/seal_dqb_phylogeny.png", figures[[5]],
#        width = 478.5,
#        height = 338.3,
#        units = "mm",
#        dpi = 400)
# # größe stimmt noch nicht


#####
# store data frames
mhc_analysis <- list(allele_summary = allele_summary,
                     allele_summaryX = allele_summaryX,
                     alleles = alleles,
                     clone_allele_df = clone_allele_df,
                     clone_genotype_df = clone_genotype_df,
                     genotype_fas = genotype_fas,
                     metadata_df = metadata_df,
                     metadata_df_pairs = metadata_df_pairs)


## delete everything but figures container and data frame container
# rm(list = setdiff(ls(),c("figures", "mhc_analysis")))
