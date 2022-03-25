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
