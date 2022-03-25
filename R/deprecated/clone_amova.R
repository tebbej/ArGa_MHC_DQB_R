library(poppr)
# occurences of alleles that only occur in one population
private_alleles(clone_gen)

# diversity table per population
clone_div <- poppr(clone_gen)
write.table(clone_div, file ="data/clone_diversity_estimates.txt", 
            sep ="\t", row.names = F)

names(strata(clone_gen))

# amova
clone_amova <- poppr.amova(clone_gen, 
                    ~pops/fmly,
                    within = F,
                    sep = "/")

clone_amova


