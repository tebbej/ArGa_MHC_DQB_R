arrangeGrob(clone_heatmap_paired, clone_heatmap_unpaired, ncol = 1))

library(gtable)
g2 <- ggplotGrob(clone_heatmap_paired)
g3 <- ggplotGrob(clone_heatmap_unpaired)
g <- rbind(g2, g3, size = "first")
g$widths <- unit.pmax(g2$widths, g3$widths)
grid.newpage()
grid.draw(g)
# no es funktionieren

p1 <- clone_heatmap_paired
p2 <- clone_heatmap_unpaired

aligned <- align_plots(p1, p2, align = "v")
ggdraw(aligned[[1]])
ggdraw(aligned[[2]])

plot_grid(p2, p1, ncol = 1, align = "hv", axis = "lbtr")


set.seed(123)
nr1 = 4; nr2 = 8; nr3 = 6; nr = nr1 + nr2 + nr3
nc1 = 6; nc2 = 8; nc3 = 10; nc = nc1 + nc2 + nc3
mat = cbind(rbind(matrix(rnorm(nr1*nc1, mean = 1,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc1, mean = 0,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc1, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc2, mean = 0,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc2, mean = 1,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc2, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc3, mean = 0.5, sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc3, mean = 0.5, sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc3, mean = 1,   sd = 0.5), nr = nr3))
)
mat = mat[sample(nr, nr), sample(nc, nc)] # random shuffle rows and columns
rownames(mat) = paste0("row", seq_len(nr))
colnames(mat) = paste0("column", seq_len(nc))

allele_summary <- log(allele_summary+1)

ComplexHeatmap::Heatmap(allele_summary)

ComplexHeatmap::Heatmap(allele_summary, name = "hi", row_split = rep(c("A", "B"), 31), row_km = 2)

# subset sequences and clones IDs of clone_allele_df because in there, sequences are ordered according to allele
# frequences and makes haplotype maps look prettier. also, numbers in maps and alleles coincide that way
x <- clone_allele_df %>% select(clone_var, sequence) %>%
  column_to_rownames(., var = "clone_var")  %>% 
  Biostrings::as.matrix() %>% 
  Biostrings::DNAStringSet() %>% 
  ape::as.DNAbin()
require(pegas)
h <- haplotype(x)
net <- haploNet(h)
plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 0.8)

ind.hap<-with(
  stack(setNames(attr(h, "index"), rownames(h))), 
  table(hap=ind, pop=rownames(x)[values])
)


data(woodmouse)
x <- woodmouse[sample(15, size = 110, replace = TRUE), ]
h <- haplotype(x)
net <- haploNet(h)
plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 0.8)

dat <- readDNAStringSet("data/Clones_MHC_ArGa_all_20210318.fas")
name_dat <- names(unique(genotype_fas))
new_fas <- dat[which(!is.na(match(names(dat), name_dat)))]

test <- clone_allele_df[clone_allele_df$id == "AGF11008",]

called_clones <- vector(mode = "list", length = 1)
called_clones[[1]] <- seq_along(unlist(attributes(clone_allele_df$id)[1]))
names(called_clones[[1]]) <- unlist(attributes(clone_allele_df$id)[1])


id <- as.character(unique(clone_allele_df$id))
called_clones <- lapply(id, function(x) as.character(unique(
  clone_allele_df$allele[which(
    !is.na(match(clone_allele_df$id,x)))]))) %>%
  `names<-`(., id)

called_clones <- list(clone_exon = called_clones)
called_clones




ngs_freq <- ggplot(data = miseq_alleles, 
                   aes(x    = ngs_alleles, 
                       y    = freq, 
                       fill = freq)) +
  geom_col(fill = "#2de6e7") +
  geom_col(aes(x    = ngs_alleles, 
               y    = -freq), 
           fill = "red") +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                     limits = c(-.20,.20),
                     expand = c(0,0)) +
  
  ylab("Frequency\n") +
  labs(fill = "Allele\ncount") + 
  theme_minimal() +
  theme(panel.grid         = element_line(color = "white"),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line          = element_line(color = "black"),
        axis.text          = element_text(color = "black"),
        axis.title         = element_text(color = "black", 
                                          margin = margin(10,10,20,10)),
        axis.ticks         = element_line(color = "black", 
                                          size = 0.2),
        axis.line.x        = element_line(color = "black"),
        axis.ticks.x       = element_blank(),
        axis.title.x       = element_blank(),
        axis.text.x.bottom = element_text(angle = 45, 
                                          vjust = 1, 
                                          hjust = 1, 
                                          size = 11),
        axis.line.y        = element_line(color = "black"),
        axis.title.y       = element_text(size = 15.5),
        axis.ticks.y       = element_blank(),
        axis.ticks.length  = unit(.15,"cm"),
        plot.background    = element_rect(color = "white", 
                                          fill = "white"),
        legend.position    = "right",
        plot.margin        = unit(c(0.5,0.5,0.5,0,5), "cm")
  )
ngs_freq
# ggsave("graphics/allele_frequency_ngs_artemis_preliminary_scale.png", ngs_freq, dpi = 400)






nwk <- system.file("extdata", "data/phylogeny/seal_dqb.nwk", package="treeio")
tree <- read.tree(nwk)
p <- ggtree(tree)


clone_genotype_df <- allele_summary[allele_summary$alleles %in% levels(allele_summary$alleles)[1:19],]
clone_genotype_df <- clone_genotype_df[which(clone_genotype_df$counts != 0),]
clone_genotype_df <- mutate(clone_genotype_df, variant_no = clone_allele_df$variant_no[match(clone_genotype_df$alleles, clone_allele_df$allele)],
                            freq = clone_allele_df$allele_frequency[match(clone_genotype_df$alleles, clone_allele_df$allele)])

length(which(clone_genotype_df$alleles %in% "ArGa-DQB*1") == T)

ggplot(data = clone_genotype_df, aes(x = colony, y = freq)) +
  geom_boxplot()

# fst calculation and permutations
# parameters in r package diveRsity?
# better AMOVA in Arlequin

# test_index <- vector(
#   mode   = "integer", 
#   length = length(genotype_fas)
#   )
# for (seqs in alleles$seq[1]) {
#   placer <- which(match(genotype_fas, seqs) == T)
#   test_index[placer] <- placer
# }


test_allele_df <- as.data.frame(genotype_fas) %>%
  transmute(sequence = x) %>%
  rownames_to_column(var = "clone_var")


test <- sapply(test_allele_df$sequence, function(x) match(x, alleles$seq))










#### plot both, miseq and clone alleles and their freq
test_figure <- egg::ggarrange(figures[[6]],
                              figures[[7]],
                              ncol=1)

test_figure <- ggpubr::ggarrange(figures[[6]],
                                 figures[[7]],
                                 ncol=1)

ggpubr::as_ggplot(test_figure)



ttttest <- clone_genotype_df %>%
  transmute(., c_allele      = alleles,
            c_var_counter = variant_counter,
            c_var_no  = variant_no)

trrry <- miseq_alleles_expand
trrry[101:102,] <- c(NA)

ttttest <- ttttest %>%
  mutate(., mi_allele = trrry$ngs_alleles,
         mi_var_counter = trrry$variant_counter,
         mi_var_no = trrry$variant_no)

ggplot(ttttest, 
       aes(x = c_var_no, 
           group = desc(c_var_counter), 
           fill = c_var_counter)) + 
  geom_bar(aes(y = stat(count) / sum(count))) + 
  geom_bar(aes(x = mi_var_no, y = -(stat(count) / sum(count)))) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                     limits = c(-0.2, 0.2), #swapped x & y values for reversing!
                     expand = c(0,0)) +
  scale_fill_viridis_c(option = "cividis",
                       begin = 0,
                       end = 1) +
  ylab("MiSeq data") +
  labs(fill = "Allele\ncounts") + 
  scale_x_continuous(breaks = seq_along(unique(miseq_alleles$ngs_alleles)),
                     labels = str_sort(unique(miseq_alleles$ngs_alleles), numeric = T),
                     # expand = c(0, 0.3)) +  
                     theme_minimal() +
                       theme(panel.grid = element_line(color = "white"),
                             panel.grid.minor = element_blank(),
                             panel.grid.major.x = element_blank(),
                             axis.line = element_line(color = "black"),
                             axis.text = element_text(color = "black"),
                             axis.title = element_text(color = "black", 
                                                       margin = margin(10,10,20,10)),
                             axis.ticks = element_line(color = "black", 
                                                       size = 0.2),
                             axis.line.x = element_line(color = "black"),
                             axis.ticks.x = element_blank(),
                             axis.title.x = element_blank(),
                             axis.text.x.bottom = element_text(angle = 45, 
                                                               vjust = 1, 
                                                               hjust = 1, 
                                                               size = 11),
                             axis.line.y = element_line(color = "black"),
                             axis.title.y = element_text(size = 15.5),
                             axis.ticks.y = element_blank(),
                             axis.ticks.length = unit(.15,"cm"),
                             plot.background = element_rect(color = "white", 
                                                            fill = "white"),
                             legend.position = "none",
                             plot.margin = unit(c(0,0.5,0.5,0.5), "cm")
                       )
                     figures[[7]]
                     
                     freq_compare <- data.frame(clone_freq = alleles$frequency[1:19]/100,
                                                ngs_freq =miseq_alleles$freq)
                     
                     ngs_clone_correlation <- ggplot(freq_compare, aes(x=ngs_freq, y=clone_freq)) +
                       geom_point() +
                       geom_smooth(method = "lm", se = F) + 
                       scale_y_continuous(limits = c(0,0.2),
                                          expand = c(0,0)) +
                       scale_x_continuous(limits = c(0,0.2),
                                          expand = c(0,0))
                     ngs_clone_correlation
                     
                     
                     
                     
                     ### plot heatmaps by pairs function:
                     sep_heat <- function(x,y){ 
                       df <- allele_summaryX[which(allele_summaryX$family == x),]
                       gg <- ggplot( data = df,
                                     aes(x = alleles,
                                         y = sample_id, 
                                         fill = log(counts+1))) + 
                         geom_tile() + 
                         coord_fixed(ratio = 0.6) +
                         scale_fill_viridis_c() +
                         # xlab("Alleles") +
                         # ylab("SSB\n") +
                         labs(fill = "Log \nclone \nnumber") +
                         theme_minimal() +
                         theme(
                           panel.grid = element_blank(),
                           axis.text.x.bottom = element_blank(),
                           axis.title.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.line.y = element_blank(),
                           # axis.line.y = element_line(color = "blue",
                           #                            size = 2),
                           axis.title.y = element_blank(),
                           axis.ticks.y = element_line(),
                           # axis.text.y  = element_text(color = color_ssb),
                           plot.margin  = unit(c(-.1, 0, 0, 0), "cm"),
                           legend.position = "none"
                         )
                       
                       return(gg)
                       
                     }
                     
                     vec <- unique(allele_summaryX$family)
                     plot_list <- lapply(vec, sep_heat)
                     test_heat <- ggarrange(plotlist = plot_list, ncol = 1)
                     ggsave(filename = "graphics/test_heat.png", plot = test_heat, dpi = 400)
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     # convert clone_genind genind2df into inbreedR object: #
                     # no functions work there for only one locus though
                     func <- function(s) {
                       s = as.character(s)
                       a <- str_split(s, pattern = "/")[[1]][1]
                       b <- str_split(s, pattern = "/")[[1]][2]
                       d <- c(a,b)
                       return(d)
                     }
                     
                     t <- matrix(nrow = length(clone_genind), ncol = 2) %>%
                       `row.names<-`(names(clone_genind))
                     
                     t <- sapply(clone_genind, func) %>%
                       t() %>%
                       as.data.frame() %>%
                       `rownames<-`(rownames(clone_genind)) %>%
                       convert_raw()
                     
                     