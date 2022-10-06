library(ape)
library(ggtree)
library(treeio)
library(magrittr)
library(tidyverse)

# nwk <- system.file("extdata", "data/phylogeny/tree.nwk", package="treeio")

tree <- ape::read.tree("data/phylogeny/tree.nwk")
# tree <- root(tree, "AM259941.1_Homo_Sapiens", resolve.root = T)
tree <- root(tree, node = 130, edgelabel = F, resolve.root = T)

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

tip_df %<>% mutate(species = as.factor(species))

# Seemingly: for `%<+%` to work and join information of tip label with
# metadata, the first column of the assigned metadata frame and the tree$tip.label
# need to be same. Thus, we assign the new labels first and only afterwards join
# the metadata with the tree `Formal class treedata`
tree$tip.label <- tip_df$label

# dataframe needs update, too! Drop unused names and correct entries
tip_drop <- tip_df$label[63:67]

tree <- drop.tip(tree, tip_drop)
tree$tip.label[63] <- "Canis lupus familiaris"
tip_df <- tip_df[1:63,]
tip_df[63,] <- c("Canis lupus familiaris", "AF016904-AF016909", "CaFa")

# PHYLOGENY

species_names <- c("Arctocephalus forsteri", 
                   "Arctocephalus gazella", 
                   "Canis lupus familiaris", 
                   "Halichoerus grypus", 
                   "Mirounga angustirostris",
                   "Mirounga leonina", 
                   "Monachus schauinslandi", 
                   "Neophoca cinerea", 
                   "Odobenus rosmarus",
                   "Phocarctos hookeri",
                   "Zalophus californianus", 
                   "Zalophus wollebaeki")

species_palette <- RColorBrewer::brewer.pal(12, "Set3")
species_palette[c(2,3, 12)] <- c("black", "#a76437","red")

# create initial tree
{
  p <- ggtree(tree, 
              layout = "circular",
              ladderize = F) 
  # provide metadata: `%<+%`
  p <- p %<+% tip_df + 
    # geom_aline() +
    geom_tiplab2(size = 5, 
                 align = F, 
                 hjust = -0.05) +
    geom_rootedge(rootedge = 0.05) +
    scale_color_manual(values = species_palette,
                       labels = species_names) +    
    geom_tippoint(aes(color = species),
                  size = 3) +
    # scale_color_viridis_c() +
    geom_treescale(x = 0,
                   linesize = 1.3,
                   offset = 1) +
    xlim_tree(.35) +
    labs(color = "Species") + 
    theme(
      legend.position = "right"
    )
  p + guides(color = guide_legend(title.position = "top"))
} 

addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize),
                                ncol = 3,
                                label.hjust = 0),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_blank(),
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

p <- addSmallLegend(p, pointSize = 3, textSize = 14)
p


ggsave("graphics/phyl_new.png", dpi = 400, width = 50.24, height = 32.88, units = "cm")
