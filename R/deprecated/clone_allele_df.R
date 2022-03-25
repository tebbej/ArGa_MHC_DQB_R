clone_allele_df <- as.data.frame(genotype_fas) %>%
  transmute(sequence = x) %>%
  rownames_to_column(var = "clone_var")

allele_index_in_df <- as.vector(sapply(clone_allele_df$sequence, function(x) match(x, alleles$seq)))

clone_allele_df <- clone_allele_df %>%
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

#update clone_allele_df with new info
index <- match(as.character(clone_allele_df$id), metadata_df$clone_id)
clone_allele_df <- clone_allele_df %>%
  mutate(
    id       = metadata_df$real_id[index],
    colony   = metadata_df$colony[index],
    maturity = metadata_df$maturity[index],
    family   = metadata_df$family[index]
  ) %>%
  relocate(., sequence, .after = last_col())