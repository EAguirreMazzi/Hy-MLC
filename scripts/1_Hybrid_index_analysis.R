rm(list = ls())
source("scripts/Hy-MLC_functions.R")

#modify this block and the rest should run if files format and info is compatible
path_to_vcf <- "data/example_data.vcf.gz"
path_to_popmap <- "data/example_data.popmap.txt"
outpath <- "output/"
prefix <- "example_data"
source0 <- "crassifolia"
source1 <-  "yungasensis"

#### 1. load vcf/gl data and filter ####
#ADJUST THIS BLOCK if neccesary 
#(different pipelines can be used as long as the input for prepare.data is genlight object)
#load  VCF file (as genlight object)
gl<-gl.read.vcf(path_to_vcf)
gl<-gl.compliance.check(gl)
#filter loci that are at least in 50% of samples (ADJUST THIS THRESHOLD ACCORDING TO DATASET)
gl<-gl.filter.callrate(gl,threshold = 0.5, recalc = TRUE, mono.rm = TRUE)
# popmap indicating species and putative hybrids
popmap<-read.table(path_to_popmap, header = T)
# asign pop field in the gl object, use popmap info here
popmap$ind.names[match(gl$ind.names, popmap$ind.names)] == gl$ind.names # just a check
gl$pop <- factor(popmap$pop[match(gl$ind.names, popmap$ind.names)])


##### 2. Simulate hybrids using sequences from source0 and source1 ####
#for example W. crassifolia will be source0 (s0) and W. yungasensis will be s1
gl_sim <- simulate.hybrids.gl(gl=gl, s0 = source0, s1 = source1, Nhybrids = 15)
#Customizable crosses can be configured editing the simulate.hybrids.gl function in "scripts/Hy-MLC_functions.R"


#### 3. Prepare data for hybrid index analisis ####
obs_prepdata <- prepare.data(gl = gl,
        source0 = source0 , source1 = source1,
        simulated = FALSE)
sim_prepdata<- prepare.data(gl = gl_sim, simulated = TRUE)

##### 4.Run hybrid index analysis for observed and simulated datasets ####
#depending on the number of loci and samples, this might be the biggest botleneck point on the analisis 
# we set include.Source to TRUE to get hybrid indices for the parental reference individuals
                  
sim_hindlabel=  esth(data.prep.object = sim_prepdata$data.prep,
                  read.data.precols = sim_prepdata$precols,
                  include.Source = TRUE,
                  nitt=3000,burnin=2000)

obs_hindlabel=  esth(data.prep.object = obs_prepdata$data.prep,
                  read.data.precols = obs_prepdata$precols,
                  include.Source = TRUE,
                  nitt=3000,burnin=2000)

# Check that hybrid index is  succesfully estimated for all individuals
# Sometimes mcmc fails in some individuals, in this case re-running the analysis is usually enough, 
# alternatively if not one can try a smaller init.var2 value than the one included in the result objects

#### 5. Merged simulated and observed hi results as dataframes and save these results 
# for Machine learning classification
obs_hi <- data.frame(obs_hindlabel$hi)
sim_hi <- data.frame(sim_hindlabel$hi)
obs_hi$dataset <- "observed"
sim_hi$dataset <- "simulated"
all_hi <- rbind.data.frame(obs_hi, sim_hi)

write.table(all_hi, file = paste0(outpath, prefix, ".hybrid_index.txt"), quote = F, sep = "\t", row.names = F)
