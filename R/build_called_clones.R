called_clones <- vector(mode = "list", length = 1)
called_clones[[1]] <- seq_along(unlist(attributes(clone_allele_df$id)[1]))
names(called_clones[[1]]) <- unlist(attributes(clone_allele_df$id)[1])


id <- as.character(unique(clone_allele_df$id))
called_clones <- lapply(id, function(x) as.character(unique(
  clone_allele_df$allele[which(
    !is.na(match(clone_allele_df$id,x)))]))) %>%
  `names<-`(., id)

# called_clones <- list(clone_exon = called_clones)
# called_clones
# 
# genotypes_df <- lapply(called_clones, function(x){
#   a <- 6-length(x)
#   b <- c(x, rep(NA, a))
#   return(b)
# }) %>% 
#   as.data.frame() %>%
#   t() %>%
#   as.data.frame

gtypes_miseq <- read.table("data/Genotypes-Artemis_ngs2.txt", 
                           header = T, 
                           sep    = "\t",
                           na.strings = "")
gtypes_miseq <- gtypes_miseq[,-2]
test <- as.list(gtypes_miseq[,2:5])

called_alleles_artemis <- lapply(seq_along(gtypes_miseq$individual.ID), function(x){
  a <- na.omit(unlist(gtypes_miseq[x,2:5]))
  b <- na.omit(match(a, as.character(shared_alleles$clone_alleles)))
  c <- as.character(shared_alleles$ngs_alleles[b])
  return(c)
}) %>% 
  `names<-`(., gtypes_miseq$individual.ID)

called_alleles_artemis <- list(miseq_artemis = called_alleles_artemis)
