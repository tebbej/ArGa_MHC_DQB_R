rm(list = ls())
library(adegenet)
library(tidyverse)
library(hierfstat)
library(reshape2)
library(RColorBrewer)
library(scales)

# create column labels for the markers
msat.names <- vector()
for (i in 1:41){
  a.var <- paste0("sat", i, "_a")
  b.var <- paste0("sat", i, "_b")
  vars <- c(a.var, b.var)
  msat.names <- c(msat.names, vars)
}

# read in genotypes
msats <- read.table(file = "data/msats/msats_genotypes.txt", sep = "\t", row.names = NULL)[-1,] %>%
  `colnames<-`(c("id", "beach", msat.names)) %>%
  remove_rownames() %>%
  column_to_rownames(., var = "id")

# site <- msats$beach %>%
#   sapply(., function(x) ifelse(x == "1 ", "SSB", "FWB")) %>%
#   as.factor()

# remove population/colony info in df
msats <- msats[, -1]



# configure msat genotypes so each genotype consists of 3 characters

for (i in 1:dim(msats)[1]) {
  for (j in 1:dim(msats)[2]) {

    # msats[i,j] <- gsub(" ", "", msats[i,j])
    # # if there's no " " there needs to be another separator for
    # # merge.col to work, as that is substantial for the function to work
    # # correctly! Thus, keep " " for now

    # nchar == 3 is important and needs to be adjusted if whitespace is removed
    # before logical statement checks for which character lengths the df needs
    # to be adjusted !
    # as well as != "NA " whitespace needs to be adjusted
    if (nchar(msats[i,j], keepNA = F) == 3 & msats[i,j] != "NA " & !is.na(msats[i,j])) {
            msats[i,j] <- paste0("0", msats[i,j])
    } else {
      msats[i,j] <- msats[i,j]
    }

  }#end j
}# end i

# test whether there are still entries not NA and >3
h <- unlist(msats) %>%
  gsub(" ", "", .)# remove whitespace 

if(is_empty(which(nchar(h, keepNA = T) == 3 & h != "NA"))){ 
  message("Zero entries nchar > 3 that != NA")
}

# extract original rownames from df
rnames <- rownames(msats) %>%
  gsub(" ", "", .) # remove whitespace

rownames(msats) <- rnames

# function to merge columns in the df for each row
# columns must be united to a single column first
# original number of columns must be even, i.e. for each satellite marker two columns
merge.col <- function(m){
  s <- str_split(m, pattern = " ") %>%
    unlist()
  i = 1
  j = 2
  a = 1
  new_col <- vector(mode = "character", length = 41) #(length(s)/2)
  new_col[a] <- paste0(s[i], "/", s[j])
  while (!isTRUE(a == 41)) { # (length(s)/2)
    a = 1 + a
    i = 2 + i
    j = 2 + j
    new_col[a] <- paste0(s[i], "/", s[j])
  }
  return(new_col)
}

# merga all column info into a single column
msats <- unite(msats, "merge", sep = "")


# apply merging function and coerce to dataframe with dim(individuals x markers)
msats <- lapply(msats$merge, merge.col) %>% 
  as.data.frame() %>%
  `colnames<-`(rnames) %>%
  `rownames<-`(paste0("sat_", 1:41)) %>%
  t() %>%
  as.data.frame()

# Data set containing meta data about MHC clone samples of Antarctic fur seals
metadata_df <- {read.table(file   = "data/sample_list.txt", 
                           header = T) %>%
    mutate(
      real_id  = factor(real_id, 
                        levels = str_sort(real_id, 
                                          numeric = T)),
      colony   = as.factor(colony),
      maturity = as.factor(maturity),
      family   = as.factor(family)
    ) %>%
    arrange(real_id) %>%
    arrange(colony)}

# Data set containing msat genotype data of Antarctic fur seals
# pcrID:  ID used in amplifiction protocols
# pnasID: IDs that matches individuals and formatting ind msats df
# colony: Information for population membership
id_df <- read.table(file = "data/msats/msats_sample_list.txt", sep = "\t") %>%
  `colnames<-`(c("pcrID", "pnasID", "colony")) %>%
  .[-1,] # first row contains column names of the .txt file that for soem reason cannot be imported diretly while file reading

# update metadata_df with pnas IDs 
pos <- match(metadata_df$clone_id, id_df$pcrID) %>% na.omit()
# double check
metadata_df$clone_id == id_df$pcrID[pos]
any(is.na(pos))

metadata_df <- metadata_df %>%
  mutate(., pnasID = id_df$pnasID[pos])

# trim msats df to only contain individuals that have both data
# mhc cloning genotypes and msat genotypes

# check whether all individuals can be found in the msat genotypes
pos_df <- match(metadata_df$pnasID, rownames(msats))
any(is.na(pos_df)) # must be FALSE, otherwise there are non overlaps in the data

# build msats df matching names in metadata_df
msats <- msats[pos_df,] %>%
  `rownames<-`(metadata_df$real_id)

# create vector with matching colony membership for overwritten rownames
# to correctly assign colony membership in genind-object
site <- metadata_df$colony

# build genind of msats df
msats_gen <- df2genind(msats, sep = "/", 
                          ploidy = 2,
                          pop = site, 
                       NA.char = NA)

save.df <- genind2df(msats_gen, sep = "/")
write.table(save.df, file = "data/msats_genind.txt", sep = "\t")

# Calculate the percentage of complete genotypes per loci 
locmiss_msats <- propTyped(msats_gen, by = "loc")
locmiss_msats[which(locmiss_msats < 0.80)] # print loci with < 80% complete genotypes
## named numeric(0)

# Barplot
barplot(locmiss_msats, ylim = c(0,1), ylab = "Complete genotypes (proportion)", xlab = "Locus", las = 2, cex.names = 0.7)

# Check genotypes are unique
mlg(msats_gen)

msats_gen
# number of alleles per each locus
table(msats_gen$loc.fac)

#  mean allelic richness per beach across all loci
allelic.richness(genind2hierfstat(msats_gen))$Ar %>%
  apply(MARGIN = 2, FUN = mean) %>% 
  round(digits = 3)

# basic stats
basic_msats <- basic.stats(msats_gen, diploid = TRUE)

# mean observed heterozygosity per beach
Ho_msats <- apply(basic_msats$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
Ho_msats

# mean expected heterozygosity per beach
He_msats = apply(basic_msats$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
He_msats


# Create a data.frame of site names, Ho and He and then convert to long format
Het_msats_df = data.frame(Site = names(Ho_msats), Ho = Ho_msats, He = He_msats) %>%
  melt(id.vars = "Site")

# Custom theme for ggplot2
custom_theme = theme(
  axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, face = "bold"),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 12),
  axis.title.x = element_blank(),
  axis.line.y = element_line(size = 0.5),
  legend.title = element_blank(),
  legend.text = element_text(size = 12),
  panel.grid = element_blank(),
  panel.background = element_blank(),
  plot.title = element_text(hjust = 0.5, size = 15, face="bold")
)

# Italic label
hetlab.o = expression(italic("H")[o])
hetlab.e = expression(italic("H")[e])

# msats heterozygosity barplot
ggplot(data = Het_msats_df, aes(x = Site, y = value, fill = variable))+
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), colour = "black")+
  scale_fill_manual(values = c("royalblue", "#bdbdbd"), labels = c(hetlab.o, hetlab.e))+
  ylab("Heterozygosity")+
  custom_theme

# inbreeding coefficient Fis
apply(basic_msats$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)

# Compute pairwise Fsts
msats_fst = genet.dist(msats_gen, method = "WC84")
msats_fst

# regular Fst
fst <- genet.dist(msats_gen, method = "Fst")
fst

# PCA
# Replace missing data with the mean allele frequencies
x = tab(msats_gen, NA.method = "mean")

# Perform PCA
pca1 = dudi.pca(x, scannf = FALSE, scale = FALSE, nf = 3)

# Analyse how much percent of genetic variance is explained by each axis
percent = pca1$eig/sum(pca1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,5),
        names.arg = round(percent, 1))

# Create a data.frame containing individual coordinates
ind_coords = as.data.frame(pca1$li)

# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2","Axis3")

# Add a column containing individuals
ind_coords$Ind = indNames(msats_gen)

# Add a column with the site IDs
ind_coords$Site = msats_gen$pop

# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords, FUN = mean)

# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))

# Define colour palette
cols = brewer.pal(nPop(msats_gen), "Set1")

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Custom theme for ggplot2
ggtheme = theme(axis.text.y = element_text(colour="black", size=12),
                axis.text.x = element_text(colour="black", size=12),
                axis.title = element_text(colour="black", size=12),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                panel.background = element_blank(),
                plot.title = element_text(hjust=0.5, size=15) 
)

# Scatter plot axis 1 vs. 2
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Site, fill = Site), size = 4, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  # custom theme
  ggtheme

