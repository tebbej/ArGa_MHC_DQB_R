df.sats <- genind2df(msats_gen, sep = "/")
df.sats <- df.sats[match(metadata_df$real_id, rownames(df.sats)),]

df.mom <- df.sats[which(metadata_df$maturity =="M"),] 
df.mom <- df2genind(df.mom[,-1], sep = "/", pop = df.mom[,1])

df.pup <- df.sats[which(metadata_df$maturity =="P"),]
df.pup <- df2genind(df.pup[,-1], sep = "/", pop = df.pup[,1])

##### df.mom diversity estimates ######
# private alleles per site (per locus)
private_alleles(df.mom) %>% apply(MARGIN = 1, FUN = sum)

# allelic richness per site (per locus)
allelic.richness(genind2hierfstat(df.mom))$Ar %>%
  apply(MARGIN = 2, FUN = mean) %>% 
  round(digits = 3)

# all kinds of basic stats
df.mom_stats <- basic.stats(df.mom, diploid = TRUE)

# mean observed heterozygosity per site
apply(df.mom_stats$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE)

# mean expected heterozygosity per site
apply(df.mom_stats$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE)

# inbreeding coefficient F_IS
apply(df.mom_stats$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE)

# pairwise F_st
genet.dist(df.mom, method = "WC")




####### df.pup diversity estimates ########
# private alleles per site (per locus)
private_alleles(df.pup) %>% apply(MARGIN = 1, FUN = sum)

# allelic richness per site (per locus)
allelic.richness(genind2hierfstat(df.pup))$Ar %>%
  apply(MARGIN = 2, FUN = mean) %>% 
  round(digits = 3)

# all kinds of basic stats
df.pup_stats <- basic.stats(df.pup, diploid = TRUE)

# mean observed heterozygosity per site
apply(df.pup_stats$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE)

# mean expected heterozygosity per site
apply(df.pup_stats$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE)

# inbreeding coefficient F_IS
apply(df.pup_stats$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE)

# pairwise F_st
genet.dist(df.pup, method = "WC")






df.clone <- genind2df(clone_gen, sep = "/")
df.clone <- df.clone[match(metadata_df$real_id, rownames(df.clone)),]

df.mom <- df.clone[which(metadata_df$maturity =="M"),] 
df.mom <- df2genind(as.matrix(df.mom[,-1]), sep = "/", pop = df.mom[,1])

df.pup <- df.clone[which(metadata_df$maturity =="P"),]
df.pup <- df2genind(as.matrix(df.pup[,-1]), sep = "/", pop = df.pup[,1])

##### df.mom diversity estimates ######
# private alleles per site (per locus)
private_alleles(df.mom) %>% apply(MARGIN = 1, FUN = sum)

# allelic richness per site (per locus)
allelic.richness(genind2hierfstat(df.mom))$Ar %>%
  apply(MARGIN = 2, FUN = mean) %>% 
  round(digits = 3)

# all kinds of basic stats
df.mom_stats <- basic.stats(df.mom, diploid = TRUE)

# mean observed heterozygosity per site
apply(df.mom_stats$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE)

# mean expected heterozygosity per site
apply(df.mom_stats$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE)

# inbreeding coefficient F_IS
apply(df.mom_stats$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE)

# pairwise F_st
genet.dist(df.mom, method = "WC")




####### df.pup diversity estimates ########
# private alleles per site (per locus)
private_alleles(df.pup) %>% apply(MARGIN = 1, FUN = sum)

# allelic richness per site (per locus)
allelic.richness(genind2hierfstat(df.pup))$Ar %>%
  apply(MARGIN = 2, FUN = mean) %>% 
  round(digits = 3)

# all kinds of basic stats
df.pup_stats <- basic.stats(df.pup, diploid = TRUE)

# mean observed heterozygosity per site
apply(df.pup_stats$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE)

# mean expected heterozygosity per site
apply(df.pup_stats$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE)

# inbreeding coefficient F_IS
apply(df.pup_stats$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE)

# pairwise F_st
genet.dist(df.pup, method = "WC")
