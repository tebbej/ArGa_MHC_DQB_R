rm(list = ls()); library(magrittr); library(tidyverse); library(ape)
source("R/functions.R")

Hamming.dist <- function(seq, ref, method = c("rel", "abs"), type = c("DNA", "Protein")) {
  method <- match.arg(method)
  type <- match.arg(type)
  # discard gaps and binding N or X if protein
  if (type == "DNA") {
    gaps_seq <- which(seq %in% c("-", "N"))
    gaps_ref <- which(ref %in% c("-", "N"))
    gaps <- unique(c(gaps_seq, gaps_ref))
  } else {
    gaps_seq <- which(seq %in% c("-", "X"))
    gaps_ref <- which(ref %in% c("-", "X"))
    gaps <- unique(c(gaps_seq, gaps_ref))
  }

  if (length(gaps) > 0) {
    seqx <- seq[-gaps]
    refx <- ref[-gaps]
  }  else {
    seqx <- seq
    refx <- ref
  }


  # estimate diff
  diff <- 0
  for (i in 1:length(seqx))  diff <- diff + ifelse(seqx[i] == refx[i], 0, 1)
  # correct for sequence length
  if (method == "rel") {
    diff <-
      ifelse(length(diff) > 0,diff/length(seqx), NA)
  }
  return(diff)
}# end Hamming.dist

seqs <- ape::read.dna("data/phylogeny/Pinnipedia-DQB-DNA.fasta", "fasta") %>%
  as.character(.) %>%
  as.matrix(.)

# 2. Sort by identity
tag <- rownames(seqs)
hd <- lapply(1:length(tag), function(x) {
  out <- data.frame(
    ref = tag[x],
    focal = tag[-x],
    hamming.dist = apply(seqs[-x,], 1,
                         Hamming.dist,
                         type = "DNA",
                         ref = seqs[x,],
                         method = "abs") %>%
      unlist())
}) %>%
  do.call("rbind",.)

df <- lapply(1:nrow(hd), function(x) {
  s <- c(hd$ref[x], hd$focal[x]) %>% sort()
  data.frame(A = s[1], B = s[2], hd = hd$hamming.dist[x])
}) %>%
  do.call("rbind",.) %>%
  unique.data.frame()

df2 <- filter(df, hd == 0)
