library(inbreedR)
library(tidyverse)
library(poppr)
library(adegenet)
library(vegan)
library(phyloseq)
source("R/perm_fst.R")

# load msats genotype df in a format easily convertible into a genind
msats_gp <- read.table(file = "data/msats_genind.txt", sep = "\t") 
sites_msats <- msats_gp[,1]
msats_gen <- msats_gp[,-1]
msats_gen <- df2genind(msats_gen, ploidy = 2, sep = "/", NA.char = NA, pop = sites_msats)

set.seed(111)
perm.fst(msats_gen, nperm = 99) # 9999 for publication

# repeat for mhc clone genotypes
clone_gp <- read.table(file = "data/clone_genind.txt", sep = "\t")
sites_clones <- clone_gp[,1]
clone_gen <- clone_gp[,-1] %>% as.data.frame() %>% 
  `rownames<-`(rownames(clone_gp)) %>%
  `colnames<-`("dqbII")
clone_gen <- df2genind(clone_gen, ploidy = 2, sep = "/", pop = sites_clones)

set.seed(111)
perm.fst(clone_gen, nperm = 99) # 9999 for publication

# calculate distance measures
# absolute genetic distance {poppr}: number of allelic difference
msats_gen.abs <- poppr::diss.dist(msats_gen) %>% as.matrix()

# euclidian distance: individual genetic distance
msats_gen.euc <- dist(msats_gen) %>% as.matrix()
clone_gen.euc <- dist(clone_gen) %>% as.matrix()

# individual genetic distance: number of loci for which individuals differ
msats_gen.loc <- ape::dist.gene(msats_gp) %>% as.matrix()
clone_gen.loc <- ape::dist.gene(clone_gp) %>% as.matrix()

# load p-distances of the MHC genotypes
clone_gen.p <- read.csv("data/p-distance_pairwise.csv", sep = ";", row.names = 1) %>%
  as.matrix()

clone_allele_df <- clone_allele_df[1:771,]
# create UniFrac distances for MHC data
allele_mat <- matrix(nrow = length(unique(clone_allele_df$id)),
                     ncol = length(unique(clone_allele_df$allele))) %>%
  `rownames<-`(., as.character(unique(clone_allele_df$id))) %>%
  `colnames<-`(., as.character(
    str_sort(
      levels(
        clone_allele_df$allele)[1:19], 
      numeric = T)))

# fill matrix with info on which and how many alleles are found in the clones for each individual fur seal
for (i in seq_along(unique(clone_allele_df$id))) {
  alleles_in_id <- summary(
    clone_allele_df$allele[clone_allele_df$id == unique(clone_allele_df$id)[i]])[1:19]
  
  allele_mat[i, ] <- alleles_in_id[str_sort(names(alleles_in_id), numeric = T)]
}

allele_mat <- ifelse(allele_mat != 0, 1, 0) %>% 
  t() 
allele_mat <- allele_mat[1:19,]

phyloseq_tree <- ape::read.tree("data/unifrac_tree_p.nwk")
plot(phyloseq_tree)

arga_phylseq <- otu_table(allele_mat, taxa_are_rows = T)
arga_phylseq <- merge_phyloseq(arga_phylseq, phyloseq_tree)

clone_gen.ufrac <- UniFrac(arga_phylseq, weighted = F) %>% as.matrix()

# correlate pairwise p-distances with pairwise number of allelic difference

set.seed(111)
mantel.corr <- vegan::mantel(clone_gen.ufrac,msats_gen.abs)
mantel.corr

# load msats genotype df to inbreedR and calculate smlh and identity 
# disequilibrium g2

df <- read.table("data/msats/msats_genotypes_inbreedR.txt", sep = "\t") %>%
  convert_raw()

check_data(df)

# standardized multi loc het
sMLH_res <- sMLH(df)
hist(sMLH_res)

# g2
g2 <- g2_microsats(df, nperm = 1000, nboot = 100 )
plot(g2)

# convert mhc data into categorical hom/het values: T/F

clones_het <- read.table(file = "data/clone_mhc_het.txt", sep = "\t")

n <- names(sMLH_res)
n_c <- rownames(clones_het)
n_in <- match(n, n_c)
clones_het <- clones_het[n_in,]

# proportion of het individuals
# het = 1, hom = 0, thus, summing het gives the correct number of het samples
hets <- sum(clones_het$het)
n_ind <- dim(clones_het)[1]
(perc <- (hets/n_ind)*100)

corr_het <- cbind(sMLH_res, clones_het$het) %>%
  `colnames<-`(c("smlh","mhc_het")) %>%
  as.data.frame()

# glm with binomially distributed data
het_glm <- glm(cbind(corr_het$mhc_het, 1-corr_het$mhc_het) ~ corr_het$smlh, family = "binomial")
summary(het_glm)
anova(het_glm, test = "Chisq") # test glm; chi square due to binomial data
chi_glm <- qchisq(1-0.1551, df = 54, lower.tail = T)

# test with hypothetical 

corr_plot <- ggplot(data = corr_het, 
                    aes(y = smlh, x = mhc_het)) + 
  geom_jitter(height = 0.02, width = 0.07) +
  geom_smooth(method = "glm", color = "orange", alpha = 0.4) +
  ylab("sMLH") +
  scale_x_continuous(name = "MHC-DQBII\nheterozygosity",
                     breaks = c(0,1),
                     labels = c("hom", "het")) +
  theme_light(base_size = 15) +
  theme(
   panel.grid = element_blank() 
  )
corr_plot
