get_clone_alleles <- function(seq_input_fas, genotype_fas){
  # get_clone_alleles compares a list of unsorted genetic sequences with a list of known genetic sequences. The function allows previously unknown alleles to be
  # sorted into the list of known alleles after passing control barriers to check uniqueness of the new sequence(s) in the data set. Output is a boolean matrix with 1
  # indicating row x column associationg the presence of a known sequence on that matrix position or 0 when there is no sequence correspondence.
  
  
  # convert input data into conventions used for this function
  miseq_data <- seq_input_fas
  clone_genotypes <- genotype_fas
  # remove unneeded tmp data
  rm(seq_input_fas, genotype_fas)

  # for-loop preambula and dependant variables
  # pre-allocate a correct data matrix that stores information from the for-loop
  miseq_in_genotypes <- matrix(nrow = length(clone_genotypes), ncol = length(miseq_data))
  # make data identifiable by their names in the matrices' row and column
  colnames(miseq_in_genotypes) <- names(miseq_data)
  rownames(miseq_in_genotypes) <- names(clone_genotypes)
  # set the starting point for an increasing counter that runs in the loop outside of the loop
  iteration_count <- c(1)
  
  # for-loop to check which sample in genotype_fas/clone_genotypes contains an allele listed in seq_input_fas/miseq_data and which one
  # runs columnwise through each sample listed in clone_genotypes
  for (seqs in as.vector(miseq_data)) {
    # fills data matrix sequentially with information 1 = 'allele found' and 0 = 'allele not found'
    miseq_in_genotypes[ ,iteration_count] <- vcountPattern(seqs, clone_genotypes)
    # increase iteration counter to not overwrite any data
    iteration_count <- iteration_count + 1
  } #end for
  
  if (any(apply(miseq_in_genotypes, 1, sum) == 0)){ #if 1
  new_seqs <- unique(clone_genotypes[which(apply(miseq_in_genotypes, 1, sum) == 0)])
  message("Unknown allele(s) detected")
  
    if (!any(new_seqs %in% miseq_data)) {
      names(new_seqs) <- as.vector(lapply(seq_along(as.vector(new_seqs)), function(i) paste0( "ArGa-DQB*", c(i+length(miseq_data)))))
      update_alleles <- c(miseq_data, new_seqs)
      message(paste(names(update_alleles)[(1+length(miseq_data)):length(update_alleles)], "has not been identified before\n"))
    } # end if
  
  return(get_clone_alleles(seq_input_fas = update_alleles, genotype_fas = clone_genotypes))
  
  } else {
    return(list(miseq_matrix = miseq_in_genotypes,
                dqb_alleles = miseq_data))
  }
}# end function 
