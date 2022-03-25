# Mon Feb 21 16:57:39 2022 ------------------------------

# Allele detection curve: Update including cloning data
# ------------------------------------------------------------------------------

## Load libraries 
# ------------------------------------------------------------------------------
rm(list = ls())
library(ape)
library(dplyr)
library(ggplot2)
library(magrittr)
library(patchwork)
library(EnvStats)
library(ggpubr)
library(gridExtra)
# ------------------------------------------------------------------------------

## Define some functions
# ------------------------------------------------------------------------------
## calculate pairwise difference to primer sequences.
## Optional, account for variable alignment length
## 
Hamming.dist <- function(seq, ref, method = c("rel", "abs")) {
  method <- match.arg(method)
  # discard gaps and binding N
  gaps_seq <- which(seq %in% c("-", "N"))
  gaps_ref <- which(ref %in% c("-", "N"))
  gaps <- unique(c(gaps_seq, gaps_ref))
  
  seqx <- seq[-gaps]
  refx <- ref[-gaps]
  
  # estimate diff
  diff <- 0
  for (i in 1:length(seqx))  diff <- diff + ifelse(seqx[i] == refx[i], 0, 1)
  # correct for sequence length
  if (method == "rel") {
    diff <- 
      ifelse(length(diff) > 0,diff/length(seqx), NA)
  } 
  return(diff)
}# end Hamming.dist

## Pick alleles based on hamming value threshold
simulate_hoelzel <- function(data, n = 1:length(data), bs = 999, 
                             hamming = hamming_values, mismatch = 1) {
  
  hamming <- subset(hamming, x <= mismatch)
  x <- rep(n, each = bs)
  y <- lapply(x, function(temp) {
    # sample genotypes
    get <- data[sample(x = 1:length(data),
                       size = temp,
                       replace = T)] %>%
      unlist() %>%
      unique() 
    # keep alleles with < mismatch differences
    keep <- get[get %in% rownames(hamming)] %>%
      length()
  }) 
  
  df <- data.frame(x = x,y = unlist(y))
  df$x <- as.factor(df$x)
  return(df)
}# end simulate_hoelzel

#' @description Summarizes data 
#' @param data a data frame
#' @param measurevar character giving column name of data to summarise
#' @param groupvars character giving column names of grouping variables
#' @param na.rm boolean
#' @param conf.interval confidence interval (default 0.95)
#' @param .drop boolean
#'
#' @source
#' Taken from the R cookbook (cookbook-r.com/Manipulating_data/Summarizing_data/)
#'
summary_stats <- function(data = NULL, measurevar = NULL, groupvars = NULL, na.rm = TRUE, conf.interval = 0.95, .drop = TRUE) {
  
  length2 <- function(x, na.rm = FALSE) {
    if (na.rm) {
      sum(!is.na(x))
    } else {
      length(x)
    }
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop = .drop,
                       .fun = function(xx, col) {
                         c(N = length2(xx[[col]], na.rm = na.rm),
                           mean = mean(xx[[col]], na.rm = na.rm),
                           sd = sd(xx[[col]], na.rm = na.rm)
                         )
                       },
                       measurevar
  )
  
  # Rename the "mean" column
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N - 1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
## SORT ALIGNMENT OF SEQUENCE VARIANTS
# 1. Load data
seqs <- ape::read.dna("C:/Users/Jonas/sciebo/MHC/Alignments/MiSeq_Clones_Hoelzel_20211027.fas",
                      format = "fasta", as.character = T) %>% 
  apply(.,2, toupper) %>% 
  cbind(.,"-")  # dummy variable

# 2. Sort by identity
tag <- rownames(seqs)
hd <- lapply(1:length(tag), function(x) {
  out <- data.frame(
    ref = tag[x],
    focal = tag[-x],
    hamming.dist = apply(seqs[-x,], 1,
                         Hamming.dist,
                         ref = seqs[x,],
                         method = "abs") %>%
      unlist())
}) %>% 
  do.call("rbind",.)


## Putative alleles Cloning sequences (full exon, 267bp)
Clones <- ape::read.dna("C:/Users/Jonas/sciebo/MHC/Clones/ArGa_DQB-Hoelzel-primer-clones_20211027.fas",
                         format = "fasta") %>%
  as.character() %>%
  apply(.,2, toupper) %>% ## append a dummy colum
  cbind(., "-")
# ------------------------------------------------------------------------------

## Extract and remove primer from the alignment
ClonesPrimer <- Clones[1,]

## remove primer from matrix
Clones <- Clones[-1,]

Clones_hd <- data.frame(x = apply(Clones, 1,
                                  Hamming.dist,
                                  ref = ClonesPrimer,
                                  method = "abs") %>%
                          unlist())

# ------------------------------------------------------------------------------
Clones_glm_df <- data.frame(mismatches = Clones_hd$x,
                           binom = 0,
                           a_counts = c(145, 72, 62, 60, 57, 55, 54, 52, 48, 46, 18, 17, 17, 14, 13, 13, 12, 9, 7))
## Set alleles characterised in Hoelzel et al to 1
Clones_glm_df$binom[c(6, 17)] <- 1 

set.seed(98)
(boxplot.B <- ggplot(Clones_glm_df, aes(x = as.factor(binom), y = mismatches, fill = as.factor(binom))) +
    geom_boxplot(alpha = 0.9,
                 fatten = 3, 
                 outlier.shape = NA ) +
    geom_jitter(aes(size = a_counts), 
                shape = 21,
                alpha = 0.9,
                width = 0.4,
                height = 0.05,
                color = "black", 
                fill = "grey") +
    scale_size(range = c(3,7)) +
    theme_classic(base_size = 26,
                  base_line_size = 1,
                  base_rect_size = 1) +
    scale_x_discrete(name = "Allele detected in both studies",
                     labels = c("No", "Yes")) +
    ylab("Mismatches at primer binding site") +
    labs(tag = "A") +
    scale_fill_manual(values = c("#FDE725FF", "#481567FF")) + 
    theme(axis.ticks = element_line(color = "black"), 
          axis.line = element_line(color = "black"),
          axis.text = element_text(color = "black"),
          legend.position = "none") 
    # stat_n_text()
)

# ggsave(boxplot.B, filename = "graphics/HammingDistanceFigure-NEW.tiff",
#        device = "tiff", width = 9, height = 6, units = "in", dpi = 400)
## ------------------------------------------------------------------------------


## Allele detection curves on simulated primer-mismatches
## -----------------------------------------------------------------------------
## Load called alleles
load("data/called_clones-20211027.RData")
clone_genotypes <- called_clones

## Simulate datasaets
clone_simul <- lapply(0:max(Clones_hd), function(x) {
  simulate_hoelzel(data = clone_genotypes[["clone_exon"]],
                   bs = 99,
                   hamming = Clones_hd,
                   mismatch = x)
})

for (i in 1:length(clone_simul)) {
  clone_simul[[i]]$mismatches <- as.character(i - 1)
}

clone_simul <- do.call("rbind", clone_simul) 
clone_simul$x <- as.numeric(as.character(clone_simul$x))
clone_summary <- summary_stats(clone_simul,
                               measurevar = "y",
                               groupvars = c("x","mismatches"),
                               conf.interval = 0.99)
## add number of Hoelzel et al., 1999
clone_summary[nrow(clone_summary) + 1, ] <- c(13, 99, 0, 4, 0, 0, 0) 
clone_summary$x <- as.numeric(as.character(clone_summary$x))

hoelzel.exp <- c(expression("Hoelzel " *italic("et al.")))
(plot.B <- ggplot(clone_summary, aes(x,y)) +
    geom_linerange(ymin = clone_summary$y - clone_summary$sd, ymax = clone_summary$y + clone_summary$sd,
                   col = "grey0", alpha = 0.4) +
    geom_point(aes(shape = mismatches),  size = 4, fill = "black") + 
        xlab("Sample size") +
    ylab("Number of alleles detected") +
    scale_x_continuous(breaks = seq(0,60,5)) +
    scale_y_continuous(breaks = seq(0,20,5),
                       limits = c(0,22)) +
    labs(tag = "B") +
         # + labs(caption = "Mean and CI estimates based on 999 simulation runs") +
    scale_shape_manual(labels = c("0 bp", "1 bp", "2 bp", "3 bp", "4 bp", "5 bp", hoelzel.exp),
                       breaks = c(0, 1, 2, 3, 4, 5, 99),
                       values = c(15, 0, 17, 2, 16, 1, 8)) +  
    theme_classic(base_size = 26,
                  base_line_size = 1,
                  base_rect_size = 1) +
    theme(axis.ticks = element_line(color = "black"), 
          axis.line = element_line(color = "black"),
          axis.text = element_text(color = "black"),
          #axis.text.y = element_blank(),
          #axis.ticks.y = element_blank(),
          #axis.title.y = element_blank(),
          #axis.title = element_text(size = 12),
          # legend.box = element_rect(color = "black"),
          legend.title = element_blank(),
          legend.background = element_rect(linetype = 1, 
                                           color = "black"),# element_rect(fill = NA),  
          legend.position = c(.0,1.0),
          # legend.spacing.y = unit(0, "pt"),
          legend.box.margin = margin(-5,0,0,8, "pt"),
          legend.justification = c("left", "top")) +
  guides(shape = guide_legend(ncol = 3,
                              label.hjust = 0)))

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

addSmallLegend(plot.B, pointSize = 3, textSize = 14)

# ggsave(plot.B, filename = "graphics/AlleDetectionFigureNEW.tiff",
#        device = "tiff", width = 8, height = 7.14, units = "in", dpi = 400)

(panel <- ggpubr::ggarrange(boxplot.B, plot.B, nrow = 2, ncol = 1, align = "v"))
ggsave("graphics/sim_panel.png", panel, dpi = 400, width = 10, height = 16, units = "in")
