# load msats genotype df in a format easily convertible into a genind
msats_gp <- read.table(file = "data/msats_genind.txt", sep = "\t") 
sites_msats <- msats_gp[,1]
msats_gen <- msats_gp[,-1]

# repeat for mhc clone genotypes
clone_gp <- read.table(file = "data/clone_genind.txt", sep = "\t")
clone_gp <- clone_gp[match(rownames(msats_gen), rownames(clone_gp)),]
sites_clones <- clone_gp[,1]
clone_gen <- clone_gp[,-1] %>% as.data.frame() %>% 
  `rownames<-`(rownames(clone_gp)) %>%
  `colnames<-`("dqbII")

arga_gen <- cbind(msats_gen, clone_gen) %>%
  df2genind(., sep = "/", NA.char = NA, 
            pop = sites_msats)

arga_gen <- genind2df(arga_gen)

f <- function(x){
  if (is.na(x)) {
    x <- "000000"
  }
  return(x)
}

for (n in 1:ncol(arga_gen)) {
  arga_gen[,n] <- sapply(arga_gen[,n], f)
}

genepop::Fst("data/msats/arga_genepop.txt", outputFile = "data/msats/fst.txt.DIV")
genepop::test_HW("data/msats/arga_genepop.txt", outputFile = "data/msats/HW.txt.DIV")
genepop::basic_info("data/msats/arga_genepop.txt", outputFile = "data/msats/basic_stats.txt.INF")
