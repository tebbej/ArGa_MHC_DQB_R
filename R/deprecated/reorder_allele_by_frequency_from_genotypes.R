# allele reordering after occuring frequencies in the actual genotypes and not just in the cloning sequences

order_name <- as.character(unique(clone_genotype_df$alleles))

i <- match(order_name, alleles$name)
alleles_new <- alleles[i,] %>%
  mutate(., counts = clone_genotype_df$variant_counts[match(order_name, clone_genotype_df$alleles)]) %>%
  mutate(frequency = (counts/sum(counts))*100) %>%
  mutate(new_name = paste0("ArGa-DQB*", rep(1:length(order_name))))

alleles_fas <- DNAStringSet(alleles_new$seq) %>%
  `names<-`(alleles_new$new_name)

writeXStringSet(alleles_fas,
                filepath = "data/allele_order_from_genotypes.fas")

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
