require(Biostrings)
require(stringr)
require(tidyverse)

seq_input_fas <- readDNAStringSet("data/ArGa-DQB-NGS-Meinolf_20200929.fas")
seq_input_fas <- readDNAStringSet("data/ArGa-DQB-NGS_Artemis_20210301.fas")
genotype_info <- readDNAStringSet("data/Clones_MHC_ArGa_exon_20210319.fas")
genotype_info <- readDNAStringSet("data/ArGa-DQB_clone-alleles_20210430.fas")


source(file="R/get_unique_seqs.R")
unique_identifier <- get_unique_seqs(genotype_info, allowed.unique.seq = 1, max.mismatch = 0)

if (!purrr::is_empty(unique_identifier[[1]])) {
  genotype_fas <- genotype_info[-unique_identifier$index_unique_seq]
}
message(paste0("Data set contains ", length(unique(genotype_fas)), " unique sequences/alleles"))


source(file = "R/get_clone_alleles.R")
miseq_in_genotypes  <- get_clone_alleles(seq_input_fas = seq_input_fas,
                                 genotype_fas = genotype_fas)

# Results for NGS allele comparison
arga_dqb_info <- apply(miseq_in_genotypes[[1]], 2, sum)
arga_dqb_info

# genotype information
source(file = "R/annotate_genotype.R")
genotypes <- annotate_genotype(genotype_fas)
genotype_list <- as.matrix(names(unlist(genotypes)))

# next step: paternity assignment and allele frequency plots