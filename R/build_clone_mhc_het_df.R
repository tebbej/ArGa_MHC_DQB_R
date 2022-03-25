# exclude putative artefacts from main data.frame
clone_allele_df <- clone_allele_df[1:771,] 

# create a list of genotypes
called_clones <- vector(mode = "list", length = 1)
called_clones[[1]] <- seq_along(unlist(attributes(clone_allele_df$id)[1]))
names(called_clones[[1]]) <- unlist(attributes(clone_allele_df$id)[1])


id <- as.character(unique(clone_allele_df$id))
called_clones <- lapply(id, function(x) as.character(unique(
  clone_allele_df$variant_no[which(
    !is.na(match(clone_allele_df$id,x)))]))) %>%
  `names<-`(., id)

# filter out individuals that do not fit the presumed ploidy of the genotyped locus
# by deleting the least likely allele as we presume diploidy
ploidy_mismatches <- which(lapply(called_clones, length) > 2)
called_clones[ploidy_mismatches] <- lapply(called_clones[ploidy_mismatches], function(x) x[1:2])
# called_clones <- called_clones[-ploidy_mismatches]


called_clones <- lapply(called_clones, function(x){
  c(x[1], tail(x,1))
})

clones_het <- sapply(called_clones, function(x) ifelse(x[1] == x[2], 0, 1))
dat <- called_clones %>% as.data.frame() %>%
  t() %>% as.data.frame()
clones_het <- cbind(dat, clones_het) %>%
  `colnames<-`(c("a1", "a2", "het"))

write.table(clones_het, file = "data/clone_mhc_het.txt", sep = "\t")
