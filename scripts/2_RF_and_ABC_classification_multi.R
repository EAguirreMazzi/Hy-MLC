# IF you have more than one pair of species, when the species are closely related it make sense to
# train the RF and ABC can be trained with a pooled simulated dataset. 
#(this would also increase the number of training observations)
# here is an example on how to achieve this concatenatig hi outputs:
prefix <- c("pair1","pair2","pair3","pair4")
path_to_hi <- paste0(outpath, prefix, ".hybrid_index.txt")