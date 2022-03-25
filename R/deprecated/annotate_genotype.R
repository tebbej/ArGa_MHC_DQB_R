annotate_genotype <- function(seq_data){
  # annotate sample genotypes from a matrix containing information about each haplotype present in a cloning sequence
  
  # maybe create a metadata file that contains individual IDs to sort easier?
  
  sample_names <- lapply(names(seq_data), function(x) {
    stringr::str_split(x, "_")[[1]][1] %>% 
      paste0(., collapse = "_")
  }) %>% unlist() %>%
    unique()
  
  # create ID list for alle unique sample (with mutliple clones)
  id_for_genotypes <- vector("list", length(sample_names))
  names(id_for_genotypes) <- sample_names
  
  # pre-allocate new list where genotypes will be stored eventually
  genotypes <- id_for_genotypes
  for (i in seq_along(id_for_genotypes)) {
    # extract unique alleles from all clones subset from one individual
    alleles <- str_match(names(seq_data), names(id_for_genotypes[i])) %>%
    {which(. != is.na(.))} %>%
      genotype_fas[.] %>%
      as.vector() %>%
      unique() 
    
    # name alleles by:
    # finding the alleles in the unique clone sequences
    # and get their allele names and index position to re-assign the names in the genotypes.list
    allele_name_index <- lapply(as.vector(miseq_in_genotypes$dqb_alleles), function(x){vcountPattern(as.character(x), alleles, 
                                                                                                     algorithm = "auto")}) %>%
      lapply(., function(x){which(x != 0)}) %>%
      .[which(. != is_empty(.))]

    # actually re-assign the correct names to the correct allele sequences for an individual
    allele_info <- c(as.vector(unlist(allele_name_index)))
    names(allele_info) <- names(allele_name_index[which(allele_name_index != is_empty(allele_name_index))])
    names(alleles)[allele_info] <- names(allele_info)  
    
    # save the allele info to the genotype.list
    genotypes[[i]] <- alleles
    # repeat for each individual 
  }
  
  return(genotypes)
  
} #end function




