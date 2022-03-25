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

t <- t(as.data.frame(called_clones)) %>% as.data.frame()
it1 <- match(t$V1, clone_allele_df$variant_no)
it2 <- match(t$V2, clone_allele_df$variant_no)

t <- t %>%
  mutate(., n1 = clone_allele_df$allele[it1],
         a1 = clone_allele_df$sequence[it1], 
         n2 = clone_allele_df$allele[it2],
         a2 = clone_allele_df$sequence[it2],
         )

t <- t %>%
  unite("p", c(a1,a2) ,sep = "")
t <- t %>% 
  transmute(., p = p)

n <- rownames(t)

i <- match(rownames(msats_gp), rownames(t))

n <- n[i]
t <- t[i,] %>% as.matrix %>%
  `rownames<-`(n)

n <- rownames(t)

t <- DNAStringSet(t, use.names = T) %>%
  `names<-`(n)

writeXStringSet(t, "data/ArGa_p-distance.fas")



