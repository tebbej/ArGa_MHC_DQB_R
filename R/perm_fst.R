require(hierfstat)
require(ade4)

perm.fst <- function(data, nperm = 100, resamples = 50){
  x <- data
  df <- genind2df(x, sep = "/")
  resamples <- resamples
  
  # # permute
  # x <- df[sample(1:nrow(df), resamples), ]
  # x <- df2genind(as.matrix(x[,-1]), ncode = 2, pop = x[,1])
  
  # calculate pairwise fst
  mat.obs <- genet.dist(x, method = "WC") %>% as.matrix()
  mat.obs <- mat.obs[1,2]
  cat("Fst: ", mat.obs)
  
  # calculate permute fst
  nperm <- nperm
  mat.perm <- sapply(1:nperm, function(i) {
    # permute
    x <- df[sample(1:nrow(df), resamples), ]
    x <- df2genind(as.matrix(x[,-1]), sep = "/", pop = x[,1])
    
    # calculate pairwise fst
    mat.fst <- genet.dist(x, method = "WC") %>% as.matrix()
    mat.fst <- mat.fst[1,2]
  })
  
  test.rand <- as.randtest(
    sim = na.omit(sapply(1:nperm, function(i) mat.perm[i])), 
    obs = mat.obs, 
    alter="greater"
    )
  
  return(test.rand)
  
} # end perm.fst


