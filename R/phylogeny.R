library(ape)
library(ggtree)
library(treeio)
library(magrittr)
library(tidyverse)

# nwk <- system.file("extdata", "data/phylogeny/tree.nwk", package="treeio")

tree <- ape::read.tree("data/phylogeny/tree.nwk")
tree <- root(tree, "AM259941.1_Homo_Sapiens", resolve.root = T)
# tree <- root(tree, node = 64, edgelabel = F, resolve.root = T)

is.rooted(tree)
# tree <- drop.tip(tree, "AM259941.1_Homo_Sapiens")

# nwk <- read.newick("data/phylogeny/tree.nwk")



# create tip label and group information
tip_name <- tree$tip.label
tip_species <- sapply(
  tip_name, 
  function(xx){
    z <- strsplit(xx, "_")[[1]][2]
    strsplit(z, "-")[[1]][1]
  }
) %>% unname()

tip_acc <- sapply(
  tip_name, 
  function(xx)
    strsplit(xx, "_")[[1]][1]
) %>% unname()

tip_df <- paste0(tip_acc, " ", tip_species)
tip_df <- data.frame(label = tip_df, 
                       accession = tip_acc, 
                       species = tip_species) %>% 
  `rownames<-`(NULL)
tip_df[63,] <- c("AM259941.1 HoSa", "AM259941.1", "HoSa")
tip_df %<>% mutate(species = as.factor(species))

# Seemingly: for `%<+%` to work and join information of tip label with
# metadata, the first column of the assigned metadata frame and the tree$tip.label
# need to be same. Thus, we assign the new labels first and only afterwards join
# the metadata with the tree `Formal class treedata`
tree$tip.label <- tip_df$label

# PHYLOGENY
# create initial tree
{
  p <- ggtree(tree, 
              layout = "circular",
              ladderize = F) 
  # provide metadata: `%<+%`
  p <- p %<+% tip_df + 
    # geom_aline() +
    geom_tiplab2(size = 3, 
                 align = T, 
                 hjust = -0.05) +
    geom_rootedge(rootedge = 0.05) +
    geom_tippoint(aes(color = species)) +
    # scale_color_viridis_c() +
    geom_treescale(x = 0.03, y = -1) +
    theme(
      legend.position = "bottom"
    )
  p
}
