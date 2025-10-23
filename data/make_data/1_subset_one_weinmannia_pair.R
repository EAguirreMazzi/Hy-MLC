
library(vcfR)
vcfR <- read.vcfR("C:/Users/aedua/Documents/Weinmannia/PROJECTS/2_hybridization/central_andes/v5/v1/WCA470.SNPfiltr.usnps.vcf.gz")
test20 <- read.table("C:/Users/aedua/Documents/Weinmannia/PROJECTS/2_hybridization/central_andes/v6/gghybrid/hindex_min80/test20", header=T)
vcfR@gt <- vcfR@gt[,colnames(vcfR@gt) %in%  c("FORMAT", test20$INDLABEL)]
write.vcf(vcfR, file="data/example_data.vcf.gz")
write.table(data.frame(ind.names=test20$INDLABEL,pop=test20$POPID),
    file="data/example_data.popmap.txt",
    row.names= F, col.names = T, quote = F, sep = "\t")


