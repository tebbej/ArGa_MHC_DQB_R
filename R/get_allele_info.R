# get_allele_info
# function to remove unique clone from multiple single sequences within a bigger .fas-file of multiple sequences
# dependant on function get_unique_seqs.R.
# requires installed packages Biostrings, tidyverse and purrr.

# input:
# data        a DNAStringset containing multiple sequences
# allele_name a character to name alleles in the data return
# rm.unique   a logical indicating whether unique sequences in data shall be removed


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
