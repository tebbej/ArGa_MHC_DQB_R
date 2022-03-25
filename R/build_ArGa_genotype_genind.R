#### CREATE AND CONVERT A LIST OF INDIVIDUAL GENOTYPES OF ANTARCTIC FUR SEALS INTO A GENIND-OBJ.

library(adegenet)
library(genepop)
library(hierfstat)
library(tidyverse)
library(poppr)

# exclude putative artefacts from main data.frame
clone_allele_df <- clone_allele_df[1:771,] %>%
  mutate(., variant_no = str_pad(variant_no, 2, pad = "0"))

# create a list of genotypes
called_clones <- vector(mode = "list", length = 1)
called_clones[[1]] <- seq_along(unlist(attributes(clone_allele_df$id)[1]))
names(called_clones[[1]]) <- unlist(attributes(clone_allele_df$id)[1])


id <- as.character(unique(clone_allele_df$id))
called_clones <- lapply(id, function(x) 
  as.character(
    unique(
      clone_allele_df$variant_no[which(!is.na(match(clone_allele_df$id,x)))]
      )
    )
  ) %>%
  `names<-`(., id)

# filter out individuals that do not fit the presumed ploidy of the genotyped locus
# by deleting the least likely allele as we presume diploidy
ploidy_mismatches <- which(lapply(called_clones, length) > 2)
called_clones[ploidy_mismatches] <- lapply(called_clones[ploidy_mismatches], function(x) x[1:2])
# called_clones <- called_clones[-ploidy_mismatches]


called_clones <- lapply(called_clones, function(x){
  c(x[1], tail(x,1))
})



# build a data frame like
#       locusA locusB locusC
#       genotype1     11   <NA>     22
#       genotype2     11     34     22
#       genotype3     12     55     21
#       genotype4     32     15     22
# that can be coerced into a "genind"
clone_df <- lapply(called_clones, function(x)
  paste0(x, collapse = "/")) %>%
  as.data.frame(.) %>%
  t(.)



# build data frame with additional info for strata in genind class object
n <- rownames(clone_df)
ind_n <- match(n,clone_genotype_df$sample_id)
strata_df <- data.frame(
  id   = n, 
  pops = clone_genotype_df$colony[ind_n],
  mtry = clone_genotype_df$maturity[ind_n],
  fmly = clone_genotype_df$family[ind_n])


# coerce into genind
clone_gen <- df2genind(clone_df, 
                          ploidy = 2, 
                          sep = "/",
                          pop = strata_df$pops,
                          strata = strata_df)

save.df <- genind2df(clone_gen, sep = "/")
# write.table(save.df, file = "data/clone_gen.txt", sep = "\t")

################################################################################
# try analyses

# pairwise fixation index
fst <- genet.dist(clone_gen, method = "Fst")



# private alleles per site (per locus)
private_alleles(clone_gen) %>% apply(MARGIN = 1, FUN = sum)

# allelic richness per site (per locus)
allelic.richness(genind2hierfstat(clone_gen))$Ar %>%
  apply(MARGIN = 2, FUN = mean) %>% 
  round(digits = 3)

# all kinds of basic stats
clone_gen_stats <- basic.stats(clone_gen, diploid = TRUE)

# mean observed heterozygosity per site
apply(clone_gen_stats$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE)

# mean expected heterozygosity per site
apply(clone_gen_stats$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE)

# inbreeding coefficient F_IS
apply(clone_gen_stats$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE)

# pairwise F_st
clone_fst <- genet.dist(clone_gen, method = "WC")







