library(adegenet)
library(dartR)
library(gghybrid)
library(abc)
library(randomForest)
library(ggplot2)
library(ggpubr)



simulate.hybrids.gl<- function(gl=gl, s0=character(), s1=character(), Nhybrids=15){
    #subset a gl object for each of the parental species (non-admixed individuals)
    gl0<-gl.keep.pop(gl, pop.list = s0, recalc = TRUE)
    gl1<-gl.keep.pop(gl, pop.list = s1, recalc = TRUE)
    #filter loci that ended up having no observations within each species as these are uninformative
    gl0<- gl.filter.allna(gl0)
    gl1<- gl.filter.allna(gl1)
    # now keep loci that are present in both parental species
    loci_to_keep<-intersect(gl0$loc.names, gl1$loc.names)
    gl0<- gl.keep.loc(gl0, loc.list = loci_to_keep)
    gl1<- gl.keep.loc(gl1, loc.list = loci_to_keep)
    # convert to genind object to then use hybridize
    source0  <- gl2gi(gl0); source0$pop <- factor(rep("source0", length(source0$pop)))
    source1  <- gl2gi(gl1); source1$pop <- factor(rep("source1", length(source1$pop)))
    #Simulate 15 hybrids on each category
    F1  <- hybridize(source0,source1, n=Nhybrids, pop= "F1", hyb.label = "F1_")
    F2  <- hybridize(F1, F1,  n=Nhybrids, pop ="F2", hyb.label = "F2_")
    F1BX0 <- hybridize(source0,F1,  n=Nhybrids, pop ="BX0_1", hyb.label = "F1BX0_")
    F2BX0 <- hybridize(source0,F2,  n=Nhybrids, pop ="BX0_2", hyb.label = "F2BX0_")
    F1BX1 <- hybridize(source1,F1,  n=Nhybrids, pop ="BX1_1", hyb.label = "F1BX1_")
    F2BX1 <- hybridize(source1,F2,  n=Nhybrids, pop ="BX1_2", hyb.label = "F2BX1_")
    # re-pool all simulated hybrid and non-admixed (source) specimens into one assembly
    gi_repooled <- repool(source0, source1, F1, F2, F1BX0, F2BX0, F1BX1,F2BX1)
    # convert back to genlight
    gl_sim <- gi2gl(gi_repooled)
    # remove monomorphic loci
    gl_sim <- gl.filter.monomorphs(gl_sim)
    return(gl_sim)
}

prepare.data <- function(
    gl=gl, 
    source0=character(), source1=character(),
    prefix="/tmp", outpath=tempdir(),
    simulated=FALSE){
    #first generate a structure file
    gl2structure(gl ,addcolumns = gl@pop,
        outfile = paste0(prefix,".str"),
        outpath = outpath)
    #IMPORTANT TO ADJUST THE STRUCTURE FILE FORMAT
    # Before reading in the structrue file into gghybrid add the collumns header for "INDLABEL" and "POPID". 
    ## Also replaced all 9's with NA in the file for missing data. 
    pre_str<-readLines(paste0(outpath, prefix, ".str"))
    pre_str[1]<-paste("INDLABEL\tPOPID\t",pre_str[1], sep = "")
    pre_str<-gsub("\t-9", "\tNA", pre_str)
    conection<-file(paste0(outpath, prefix, ".str"))
    writeLines(text=pre_str,con=conection)

    dat = read.data(file=paste0(outpath, prefix, ".str"),
        nprecol=2, NUMLOCI= gl$n.loc, NUMINDS= length(gl$ind.names),
        INDLABEL=1,ONEROW=0, POPID=1, MARKERNAME=1, MISSINGVAL=NA)

    prepdata =data.prep(data=dat$data,
        loci=dat$loci,
        alleles=dat$alleles,
        S0=ifelse(simulated,"source0",source0), #POPID names for the first parental reference set#
        S1=ifelse(simulated,"source1",source1), #POPID names for the second parental reference set#
        precols=dat$precols,
        return.genotype.table=F,
        return.locus.table=F)
    return(prepdata)
}


read.hi.data <- function(path_to_hi= character()){
    hi <- read.table(path_to_hi, header = T)
    hi_data<- vector("list",4)
    names(hi_data) <- c("all_variables", "observed_data", "simulated_data", "model_index")
    hi_data$all_variables <- c("h_posterior_mode","h_cred_int_lower","h_cred_int_upper", "beta_mean", "beta_var", "beta_shape1","beta_shape2")
    #for convenience split simulated (training) and observed data
    hi_obs <- hi[hi$dataset =="observed",]
    hi_data$observed_data <- hi_obs
    #Prepare data and variables (a.k.a.,summary statistics) for ABC and RF
    #for the training dataset
    hi_sim <- hi[hi$dataset =="simulated",]
    hi_sim$INDLABEL[hi_sim$Source == "S0"] <- "source0"
    hi_sim$INDLABEL[hi_sim$Source == "S1"] <- "source1"
    hi_sim$INDLABEL <- sapply(strsplit(hi_sim$INDLABEL,"_"), function(x)x[1])
    #lump F1 and F2 categories since these cannot be diferentiated based on hybrid index..
    hi_sim$INDLABEL[hi_sim$INDLABEL %in% c("F1","F2")] <- "F1_F2"
    hi_sim$INDLABEL[hi_sim$INDLABEL %in% c("F1BX0","F2BX0")] <- "BX0"
    hi_sim$INDLABEL[hi_sim$INDLABEL %in% c("F1BX1","F2BX1")] <- "BX1"
    hi_data$model_index <- hi_sim$INDLABEL
    hi_data$simulated_data <- hi_sim
    return(hi_data)
}


rf.varImp <- function(hi_data=hi_data, ntree=10000, impVars=4) {
    data <- hi_data$simulated_data[,hi_data$all_variables]
    model_index <- factor(hi_data$model_index)
    rf_result <- hi_data
    #build random forest model with variable importance analysis
    rf <- randomForest(model_index ~., data=data , ntree=ntree, importance=T)
    #before predicting data with ABC we select important variables
    #select important variables
    importance<-data.frame(rf$importance)
    vars_ranked_mda<-row.names(importance)[order(importance$MeanDecreaseAccuracy, decreasing = T)][1:impVars]
    vars_ranked_mdg<-row.names(importance)[order(importance$MeanDecreaseGini, decreasing = T)][1:impVars]
    important_variables<-union(vars_ranked_mda, vars_ranked_mdg)
    #lets predict hybrid categories in the data now
    vote <- predict(rf, hi_data$observed_data[,hi_data$all_variables], type = "prob")
    vote<-data.frame(as.matrix(vote))
    vote$INDLABEL <- hi_data$observed_data$INDLABEL
    rf_result$rf <- rf
    rf_result$important_variables <- important_variables
    rf_result$rf_votes <- vote
    return(rf_result)
}

build.predict.abc <- function(rf_result=rf_result, nval = numeric(), tols=numeric()){
    abc_result <- rf_result
    #create the predictive posterior probs
    cv <-cv4postpr(rf_result$model_index, sumstat = rf_result$simulated_data[,rf_result$important_variables],
                nval = nval, method = "rejection", tols = tols)
    #determine posterior model probabilities for empirical data using all data
    Ne <-nrow(rf_result$observed_data) #number of empirical observations
    ABC.summary <- list()
    ABC.table <- matrix(NA, Ne, 5)
    for (i in 1:Ne){
        abc.emp <- postpr(target = rf_result$observed_data[i,rf_result$important_variables],
                            index = rf_result$model_index,
                            sumstat = rf_result$simulated_data[,rf_result$important_variables],
                            tol = tols, method="rejection")
        ABC.summary[[rf_result$observed_data$INDLABEL[i]]] <- summary(abc.emp)
        ABC.table[i,]<- ABC.summary[[rf_result$observed_data$INDLABEL[i]]]$Prob
    }
    ABC.table<-data.frame(ABC.table)
    colnames(ABC.table) <- names(ABC.summary[[rf_result$observed_data$INDLABEL[i]]]$Prob)
    ABC.table$INDLABEL <- rf_result$observed_data$INDLABEL
    abc_result$ABC_cv <- cv
    abc_result$ABC_posterior_predictions <- ABC.table
    return(abc_result)
}


harvest.hymlc <- function(abc_result=abc_result, export_summary=TRUE, outpath=outpath, prefix=prefix ){
    ## Categorize RESULTS FROM ABC ANC RF and append sample label to each
    categories<-c("BX0", "BX1", "F1_F2", "source0", "source1")
    inferred_hybrid_generation <- data.frame(
    INDLABEL=abc_result$observed_data$INDLABEL,
    rf=categories[unlist(apply(abc_result$rf_votes[,categories],MARGIN = 1, FUN=which.max))],
    max_rf_votes =apply(abc_result$rf_votes[,categories],MARGIN = 1, FUN=max),
    abc=categories[unlist(apply(abc_result$ABC_posterior_predictions[,categories],MARGIN = 1, FUN=which.max))],
    max_abc_prob=apply(abc_result$ABC_posterior_predictions[,categories],MARGIN = 1, FUN=max)
    )

    hymlc_full_result <- left_join(left_join(abc_result$observed_data, 
    left_join(abc_result$rf_votes,abc_result$ABC_posterior_predictions,by="INDLABEL", suffix = c(".rf_votes", ".abc_prob")),
    by= "INDLABEL"), inferred_hybrid_generation, by = "INDLABEL")

    #export results
    if(export_summary){
    write.table(hymlc_full_result,
        file=paste0(outpath,prefix, ".hymlc_results.txt"), 
        col.names = TRUE,row.names = TRUE, quote = FALSE, sep="\t")
    #export RF and ABC main diagnostic results more detail in the abc_result object
    sink(paste0(outpath,prefix, ".RF_ABC_diagnostics.txt"))
    write.table("#Random Forest confusion matrix:", col.names = F,row.names = F, quote = F)
    write.table(round(abc_result$rf$confusion,3),col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
    write.table("\n#Random Forest variable importance:", col.names = F,row.names = F, quote = F)
    write.table(abc_result$rf$importance,col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
    write.table("\n#Random Forest variable importanceSD:", col.names = F,row.names = F, quote = F)
    write.table(abc_result$rf$importanceSD,col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
    write.table("\n#Approximate Bayesian computation cross-validation error:", col.names = F,row.names = F, quote = F)
    summary(abc_result$ABC_cv, probs=T)
    sink()
    
    svg(paste0(outpath,prefix, ".RF_variable_importance.svg"))
    varImpPlot(abc_result$rf,
    main="Random Forest variable importance", 
    color= rev(ifelse(abc_result$all_variables %in% abc_result$important_variables, "blue","black"))) #this does not make sense but works
    dev.off()
    }
    
    return(hymlc_full_result)

}





plot.hymlc <- function(hymlc_result=hymlc_result, source0=character(), source1=character(),cutoff=0.75){
    #curate those ambiguous categories with posteerior below the selected cutoff
    hymlc_result$category <- ifelse(hymlc_result$max_abc_prob >= cutoff, hymlc_result$abc, "ambiguous")
    #change source(a priori) names to something more meaningful
    hymlc_result$Source
    hymlc_result$Source[hymlc_result$Source=="S0"] <- source0
    hymlc_result$Source[hymlc_result$Source=="TEST"] <- "putative_hybrid"
    hymlc_result$Source[hymlc_result$Source=="S1"] <- source1

    #define colums to use in each subplot
    gg_cols<-c("INDLABEL","category", "Source","h_posterior_mode" ,"h_cred_int_lower" ,"h_cred_int_upper")
    rf_cols<-c("INDLABEL",names(hymlc_result)[grep("\\.rf_votes", names(hymlc_result))])
    abc_cols<-c("INDLABEL",names(hymlc_result)[grep("\\.abc_prob", names(hymlc_result))])

    #define color palletes
    colors <- c(RColorBrewer::brewer.pal(5, "Dark2"), "#A9A9A9")  # Use Dark2 for sober colors
    names(colors)<- c("source0","BX0","F1_F2","BX1","source1", "ambiguous") 
    #pie(rep(1,6),col=colors)
    #define names of s0 and s1 species, and order of these
    POPID_order<- c(source0, "putative_hybrid", source1)

    #make a legend for the plots
    custom_legend <- get_legend(ggplot(data.frame(x=1:6, y=1:6,Source=c(source0, rep("putative_hybrid",3), source1,"putative_hybrid" ),category=names(colors),color=colors))+
    geom_bar(aes(x=x, fill=category))+
    geom_point(aes(x=x,y=y, shape=Source))+
    scale_shape_manual(values = setNames(c(15,16,17), POPID_order))+
    scale_fill_manual(values=colors)+theme_classic()+
    theme(legend.position = "bottom", legend.text = element_text(size=5),
    legend.title = element_text(size = 5)))

    #subset hi data from gghybrids
    tmp_gg<- hymlc_result[, gg_cols]
    # Reorder Sequence_ID based on the levels of Source
    INDLABEL_order <- tmp_gg$INDLABEL[order(match(tmp_gg$Source, POPID_order))]
    tmp_gg$INDLABEL <- factor(tmp_gg$INDLABEL, levels =INDLABEL_order )
    tmp_gg$Source <- factor(tmp_gg$Source, levels = POPID_order)
    tmp_gg$category <- factor(tmp_gg$category, levels = names(colors))

    #
    s0_line<- max(tmp_gg[tmp_gg$Source == source0,c("h_cred_int_lower","h_cred_int_upper") ])
    s1_line<- min(tmp_gg[tmp_gg$Source == source1,c("h_cred_int_lower","h_cred_int_upper") ])

    pgg<-ggplot(tmp_gg) +
    geom_point(aes(y = h_posterior_mode,x = INDLABEL, shape=Source, color=category), size=1)+
    geom_segment(aes(y = h_cred_int_lower,yend=h_cred_int_upper, x = INDLABEL,xend = INDLABEL,color=category))+
    scale_color_manual(values = colors[names(colors) %in% levels(tmp_gg$category)],guide="none")+
    scale_shape_manual(values = setNames(c(15,16,17), POPID_order))+
    geom_hline(yintercept = s0_line, linewidth=0.5, lty=2, color=colors[1])+
    geom_hline(yintercept = s1_line, linewidth=0.5, lty=2, color=colors[5])+
    ylab("Hibrid Index")+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90, vjust =0, size = 4),#, vjust = 0.5, size=5),
            axis.text.y = element_text( size=5),
            axis.title = element_text(size=5),
            legend.position = "none")
    #RF
    tmp_rf<- melt(hymlc_result[,rf_cols])
    tmp_rf$variable<-gsub("\\.","/",gsub("\\.rf_votes","",tmp_rf$variable))
    tmp_rf$variable<- factor(tmp_rf$variable, levels = names(colors))
    tmp_rf$INDLABEL <- factor(tmp_rf$INDLABEL, levels = INDLABEL_order)

    prf <- ggplot(tmp_rf, aes(fill=variable, y=value, x=INDLABEL)) + 
            geom_bar(position="fill", stat="identity")+
            scale_fill_manual(values = colors5)+
            theme_minimal()+
            ylab("Votes (RF)")+
            theme(axis.text.x = element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_text(size = 5),
            axis.text = element_text(size = 5),
            legend.position="none")

    #ABC
    tmp_abc<- melt(hymlc_result[,abc_cols])
    tmp_abc$variable<-gsub("\\.","/",gsub("\\.abc_prob","",tmp_abc$variable))
    tmp_abc$variable<- factor(tmp_abc$variable, levels = names(colors))
    tmp_abc$INDLABEL <- factor(tmp_abc$INDLABEL, levels = INDLABEL_order)

    pabc <- ggplot(tmp_abc, aes(fill=variable, y=value, x=INDLABEL)) + 
            geom_bar(position="fill", stat="identity")+
            scale_fill_manual(values = colors5)+
            theme_minimal()+
            ylab("Bayes Factor (ABC)")+
            theme_classic()+
            theme(axis.text.x = element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_text(size = 5),
            axis.text = element_text(size = 5),
            legend.position="none")
            
    p<-ggarrange(ggarrange(prf,pabc, ncol = 1, nrow = 2), pgg, ncol=1, nrow=2, heights = c(2,3),
    legend.grob = custom_legend, legend = "bottom" )

    return(p)
}
