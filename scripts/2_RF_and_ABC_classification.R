rm(list = ls())

source("scripts/Hy-MLC_functions.R")

# modify this block and the rest should run if files format and info is compatible
# if paths have not change in relation to the 1_Hybrid_index_analysis.R script no need to define the path_to_hi 
outpath <- "output/"
prefix <- "example_data"
source0 <- "crassifolia"
source1 <-  "yungasensis"
path_to_hi <- paste0(outpath, prefix, ".hybrid_index.txt")

#Load simulated and observed Hibrid index(hi) data produced by the 1_Hybrid_index_analysis.R script
hi_data <- read.hi.data(path_to_hi = path_to_hi)
#RUN RANDOM FOREST prediction with all variables and variable importance analysis
rf_result <- rf.varImp(hi_data = hi_data) 
##RUN APROXIMATE BAYESIAN COMPUTATION
#fist calculate the category with the least observations for defining the number of observattions for the cross validation
mincat<-min(aggregate(rf_result$model_index, list(rf_result$model_index), FUN=length)$x)
nval<- round(mincat*0.2) # lets take  20% of the minimum category for the cross validation,
#this percentage can be adjusted if there is more observations 10% is enough  

#NOW build the abc model and predict hybrid generation
abc_result <- build.predict.abc(rf_result = rf_result,nval = n_val, tols = 0.01)
#for defining tols one must adjust according to data here I use high value 0.01 because i have few observations.
# ideally start with tols=0.005 then adjust depending on number of observations (pay attention to warnings)
# a good rule of thumb is use the maximum value at wich the posteriors show no more change

#Finally lets harvest these results. 
#This will gather observed hi and RF+ABC classifications into single hymlc_result and export summarized classification and error rates of abc and rf 
hymlc_result <- harvest.hymlc(abc_result = abc_result, outpath = outpath, prefix = prefix)
cbind(hymlc_result$INDLABEL, hymlc_result$Source, hymlc_result$rf, hymlc_result$abc)

#plot Hy-MLC result
hymlc_plot <- plot.hymlc(hymlc_result = hymlc_result,
            source0 = source0, source1 = source1, cutoff = 0.99)

ggexport(hymlc_plot, filename=paste0(outpath,prefix, "_HyMLC_figure.svg"))
