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
               unique(clone_genotype_df$alleles), f1), 
             function(x) seq(1:x)))) %>%
  arrange(desc(variant_counts))

# allele reordering after occuring frequencies in the actual genotypes and not just in the cloning sequences

order_name <- as.character(unique(clone_genotype_df$alleles))

i <- match(order_name, alleles$name)
alleles_new <- alleles[i,] %>%
  mutate(., counts = clone_genotype_df$variant_counts[match(order_name, clone_genotype_df$alleles)]) %>%
  mutate(frequency = (counts/sum(counts))*100) %>%
  mutate(new_name = paste0("ArGa-DQB*", rep(1:length(order_name))))

# use order of matches to easily rename all alleles in clone_genotype_df$alleles
name_pos <- match(clone_genotype_df$alleles, order_name)

clone_genotype_df <- clone_genotype_df %>%
  mutate(., alleles = alleles_new$new_name[name_pos])

clone_genotype_df <- clone_genotype_df %>%
  mutate(.,variant_no       = sapply(
    clone_genotype_df$alleles, 
    function(x) {
      stringr::str_split(x, "\\*")[[1]][2] %>% 
        paste0(., collapse = "\\*") %>% 
        as.numeric()
    }))



figures[[6]] <- ggplot(clone_genotype_df, 
       aes(x = variant_no, 
           group = desc(variant_counter), 
           fill = variant_counter)) + 
  geom_bar(aes(y = stat(count) / sum(count))) + 
  scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                     limits = c(0,.2),
                     expand = c(0,0)) +
  scale_fill_viridis_c(option = "viridis",
                       begin = 0,
                       end = 1) +
  ylab("Cloning data") +
  labs(fill = "Allele\ncounts") + 
  scale_x_continuous(breaks = seq_along(unique(clone_genotype_df$alleles)),
                     labels = str_sort(unique(clone_genotype_df$alleles), numeric = T),
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
        # axis.line.x = element_line(color = "black"),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        # axis.text.x.bottom = element_text(angle = 45, 
        #                                   vjust = 1, 
        #                                   hjust = 1, 
        #                                   size = 11),
        axis.line.y = element_line(color = "black"),
        axis.title.y = element_text(size = 15.5),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(.15,"cm"),
        plot.background = element_rect(color = "white", 
                                       fill = "white"),
        legend.position = "none",
        plot.margin = unit(c(0.5,0.5,0.5,0), "cm")
  )
figures[[6]]
