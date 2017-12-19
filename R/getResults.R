##########################################
###### AUTOMATIZE READING RESULTS ########
##########################################

## INPUT :
#         1) studyType : name of the IMPACT directory where the results/data are saved
#         2) method : Method implemented to be analysed (LASSO,RF,GBM,CF)
#         3) LT : is the data left truncated : TRUE v FALSE

## OUTPUT :
#         For each of the methods to be analyzed will give
#         1) CI distribution bar plot
#         2) An influence plot of top hit genes
#         3) Kaplan meier plot based on 4 risk groups
#         4) Proportion of mutated genes in each of the 4 risk groups, per top hit gene


####### Make description file :

#' getResults
#'
#' This functions let's the user study the output of the VariableSelection function. This function takes as input
#' one of the (or the one) objects returned from the different machine learning algorithms chosen previously.
#' Only one such object can be inputted everytime in the getResults function.
#' @param VarSelectObject A list object outputed by the VariableSelection function.
#' @param numGroups The number of groups to be made when stratifying by risk groups. Options are 2,3 and 4 (for now implementing
#' broader version). Default is 2.
#' @param cuts Numeric vector of the points in the distribution of risk scores where groups will be splitted. Needs to be of length
#' numGroups - 1. eg : c(0.25,0.5,0.75) when numgroups is 4. Default is 0.5.
#' @param geneList Optional character vector of gene names to use to generate the pie charts.
#' @return ciSummary Summary of the distribution of the concordance index accross all runs.
#' @return inflPlot Bar plot of frequency of the 20 most selected features.
#' @return topHits Character vector of the top 10 most selected features.
#' @return average.risk Average predicted risk score for each patient in the data set.
#' @return data.out The data that was used for the analysis.
#' @return selectInflPlot Volcano plot of the selection frequency against the average mean coefficient of each feature accross all runs.
#' Note that this plot is interactive.
#' @return RiskRefit Refitted cox proportional hazard model with the predicted average risk score as the continuous covariate.
#' @return RiskHistogram Histogram of the density distribution of the average predicted risk score. Note it has been rescaled from 0 to 10
#' for simplicity.
#' @return Fits Data frame reporting the coefficients found for each feature at each run.
#' @return time.type Time unit used. Options are Days or Months.
#' @return RiskScoreSummary Distribution summary of the average predicted risk score.
#' @return KM_Plot Kaplan-Meier plot stratified by risk group.
#' @return SurvSum Summary of survival metrics per risk group.
#' @return mut_Plot Mutation distribution by features bar plot.
#' @return PieChart Interactive pie chart of the mutation dsitribution using either the most frequently selected features or
#' the manually inputted gene list. Each of the pies represent one of the risk groups.
#' @return GenesUsed Character vector of the features used to make the pie charts.
#' @keywords Results
#' @export
#' @examples library(OncoCast)
#' test <- OncoCast(data=survData,formula=Surv(time,status)~.,
#'                          method=c("LASSO"),
#'                          runs = 5,cores = 1,sampling = "cv",
#'                          pathResults = "./Test/",studyType = "ExampleRun",save=F)
#' out <- getResults_OC(test$LASSO,numGroups=2,cuts=0.5,geneList=NULL)
#' @import survival
#' @import ggplot2
#' @import plotly
#' @import plyr
#' @import stats
#' @import reshape2
#' @import scales
#' @import survminer
#' @import data.table


getResults_OC <- function(VarSelectObject,numGroups=2,cuts=0.5,geneList=NULL,mut.data=F){


  ############## CHECKS
  if(!(numGroups %in% 2:4)) {stop("ERROR : You must use a number of groups between 2 and 4. We are working on implementing a broader version")}
  if(max(cuts) >= 1 || min(cuts) <= 0){stop("ERROR : You must select cuts that are between 0 and 1 (Not included).")}
  if(!(VarSelectObject[[1]]$method %in% c("LASSO","RIDGE","ENET"))){"ERROR : The inputted object is NOT an output
    of the VariableSelection function."}

  try(library(plotly,lib.loc="/usr/local/lib/R/site-library/"),silent=T)
  try(library(plotly),silent=T)
  library(reshape2)
  library(ggplot2)
  library(scales)
  library(survminer)
  library(data.table)
  library(survival)
  library(plyr)

  ##### COX ######
  data <- VarSelectObject[[1]]$data
  ## determine if left truncated
  if(length(grep("time",colnames(data)))  == 1) {LT = FALSE}
  if(length(grep("time",colnames(data)))  == 2) {LT = TRUE}

  if( LT && max(data$time2) > 1000){MD = 365
  time.type = "Days"}
  if( LT && max(data$time2) < 1000){MD = 12
  time.type = "Months"}
  if( !LT && max(data$time) > 1000){MD = 365
  time.type = "Days"}
  if( !LT && max(data$time) < 1000){MD = 12
  time.type = "Months"}
  ### LASSO ANALYSIS ###

  #try(setwd("./results/"))

  if(VarSelectObject[[1]]$method %in% c("LASSO","RIDGE","ENET")) {

    ConcordanceIndex <- as.data.frame(as.vector(unlist(sapply(VarSelectObject, "[[", "CI"))))
    summary.CI <- round(as.data.frame(c(quantile(ConcordanceIndex[,1],c(0.1,0.25,0.5,0.75,0.9),na.rm = T))),digits = 2)
    colnames(summary.CI) <- "Concordance Index"
    rownames(summary.CI) <- c("Lower 10%","1st Quarter","Median","3rd Quarter","Upper 10%")
    CI.BP <- as.data.frame(t(summary.CI))

    if(LT){genes <- colnames(data)[-c(1:3)]}
    if(!LT){genes <- colnames(data)[-c(1:2)]}
    selected.genes.lasso <- as.data.frame(matrix(0L,nrow=1,ncol = length(genes)))
    colnames(selected.genes.lasso) <- genes

    for(i in 1:length(VarSelectObject)){
      temp <-VarSelectObject[[i]]$fit
      if(!is.null(temp)){
        selected.temp <- rownames(temp)
        selected.genes.lasso[,match(selected.temp,colnames(selected.genes.lasso))] <-
          selected.genes.lasso[,match(selected.temp,colnames(selected.genes.lasso))] +1
      }
    }

    topHits <- names(sort(selected.genes.lasso,decreasing = TRUE))[1:10]
    if(length(selected.genes.lasso)  >= 20){ranked.lasso <- sort(selected.genes.lasso,decreasing = TRUE)[1:20]}
    if(length(selected.genes.lasso)  < 20){ranked.lasso <- sort(selected.genes.lasso,decreasing = TRUE)[1:length(selected.genes.lasso)]}
    melt.rank.lasso <- melt(ranked.lasso)
    melt.rank.lasso$value <- melt.rank.lasso$value/length(VarSelectObject)
    colnames(melt.rank.lasso) <- c("Gene","Frequency")

    influencePlot <- ggplot(melt.rank.lasso,aes(x=Gene,y=Frequency,fill=Gene))+geom_col()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      #scale_fill_discrete(guide = guide_legend(title = "Category")) +
      #scale_fill_manual(values=rainbow(length(unique(melt.rank.lasso$Gene)))) +
      labs(title = paste("Most recurrent selected genes out of",length(VarSelectObject),"runs")) +
      theme(legend.position="none")

    if(LT) Variables <- colnames(data)[-c(1:3)]
    if(!LT) Variables <- colnames(data)[-c(1:2)]
    allCoefs <- as.data.frame(matrix(nrow=length(VarSelectObject),ncol=length(Variables)))
    colnames(allCoefs) <- Variables

    for(x in 1:length(VarSelectObject)){
      if(!is.na(VarSelectObject[[x]]$fit[1])){
        coefsValues <- VarSelectObject[[x]]$fit[,1]
        allCoefs[x,match(names(coefsValues),colnames(allCoefs))] <- as.numeric(coefsValues)}
    }
    allCoefs[is.na(allCoefs)] <- 0

    meanCoefs <- apply(allCoefs,2,function(x){mean(x,na.rm = TRUE)})
    selectFreq <- apply(allCoefs,2,function(x){
      length(which(x!=0))/length(x)
    })

    ## get mu freq
    if(LT) data.temp <- data[,-c(1:3)]
    if(!LT) data.temp <- data[,-c(1:2)]
    MutationFrequency <- apply(data.temp,2,function(x){
      sum(x)/length(x)
    })

    resultsAll <- as.data.frame(cbind(meanCoefs,selectFreq,MutationFrequency))
    colnames(resultsAll) <- c("MeanCoefficient","SelectionFrequency","MutationFrequency")
    rownames(resultsAll) <- names(meanCoefs)
    resultsAll <- resultsAll[complete.cases(resultsAll),]
    resultsAll$GeneName <- rownames(resultsAll)
    resultsAll$MutationFrequency2 <- cut(resultsAll$MutationFrequency, c(0,0.10,0.20,0.40))

    selectInflPlot <- plot_ly(data = resultsAll, x = ~MeanCoefficient, y = ~SelectionFrequency,
                              text = ~paste('Gene :',GeneName,
                                            '</br> Hazard Ratio :',round(exp(MeanCoefficient),digits=2)),
                              mode = "markers",size = ~MutationFrequency,color = ~MutationFrequency) %>%
      layout(title ="Volcano Plot")

    final.pred <- as.data.frame(matrix(nrow= nrow(data),ncol = length(VarSelectObject)))
    rownames(final.pred) <- rownames(data)
    for(i in 1:length(VarSelectObject)){
      temp <- VarSelectObject[[i]]$predicted
      final.pred[match(names(temp),rownames(final.pred)),i] <- as.numeric(temp)
    }

    average.risk <- apply(final.pred,1,function(x){
      mean(as.numeric(x),na.rm = TRUE)
    })
    average.risk[which(is.na(average.risk))] <- NA
    RiskScore <- rescale(as.numeric(average.risk), to = c(0, 10), from = range(average.risk, na.rm = TRUE, finite = TRUE))
    summary.RiskScore <- round(as.data.frame(c(quantile(RiskScore,c(0.1,0.25,0.5,0.75,0.9),na.rm = TRUE))),digits = 2)
    colnames(summary.RiskScore) <- "Risk Score"
    rownames(summary.RiskScore) <- c("Lower 10%","1st Quarter","Median","3rd Quarter","Upper 10%")
    ## refit coxph model with average risk as covariate
    meanRS <- mean(RiskScore)
    #RiskScore <- average.risk #- meanRS
    if(LT) refit.risk <- coxph(Surv(data$time1,data$time2,data$status)~RiskScore)
    if(!LT) refit.risk <- coxph(Surv(data$time,data$status)~RiskScore)
    Risk <- as.data.frame(RiskScore)

    RiskHistogram <- ggplot(Risk, aes(x = RiskScore, y = ..density..)) +
      geom_histogram(show.legend = FALSE, aes(fill=..x..),
                     breaks=seq(min(Risk$RiskScore,na.rm = T), max(Risk$RiskScore,na.rm = T), by=0.05)) +
      geom_density(show.legend = FALSE) +
      theme_minimal() +
      labs(x = "Average risk score", y = "Density") +
      scale_fill_gradient(high = "red", low = "green")
  }


  #####################################
  ####### RISK STRATIFICATION #########
  #####################################


  data$lvl4Groups <- rep(NA,nrow(data))
  qts <- as.numeric(quantile(average.risk,cuts,na.rm = T))

  if(numGroups ==2){
    for(i in 1:nrow(data)){

      if(is.character(try(average.risk[i] < qts,silent=T))){
        stop("ERROR : Increase the number of runs. Some patients were never attributed to the testing set. Recommended value is above 20")
      }

      if(average.risk[i] < qts) {
        data$lvl4Groups[i] <- "Low"
      }

      if(average.risk[i] >= qts) {
        data$lvl4Groups[i] <- "High"
      }
    }
  }

  if(numGroups == 3){

    for(i in 1:nrow(data)){

      if(is.character(try(average.risk[i] < qts,silent=T))){
        stop("ERROR : Increase the number of runs. Some patients were never attributed to the testing set. Recommended value is above 20")
      }

      if(average.risk[i] < qts[1]) {
        data$lvl4Groups[i] <- "Low"
      }
      if(average.risk[i] >= qts[1] && average.risk[i] < qts[2]){
        data$lvl4Groups[i] <- "Intermediate"
      }
      if(average.risk[i] >= qts[2]) {
        data$lvl4Groups[i] <- "High"
      }
    }
  }

  if(numGroups == 4){
    for(i in 1:nrow(data)){

      if(is.character(try(average.risk[i] < qts,silent=T))){
        stop("ERROR : Increase the number of runs. Some patients were never attributed to the testing set. Recommended value is above 20")
      }

      if(average.risk[i] < qts[1]) {
        data$lvl4Groups[i] <- "Low"
      }
      if(average.risk[i] >= qts[1] && average.risk[i] < qts[2]){
        data$lvl4Groups[i] <- "Low-intermediate"
      }
      if(average.risk[i] >= qts[2] && average.risk[i] < qts[3]){
        data$lvl4Groups[i] <- "High-intermediate"
      }
      if(average.risk[i] >= qts[3]) {
        data$lvl4Groups[i] <- "High"
      }
    }
  }


  ##################

  ### Is this left truncated ?
  # Will be LT if second column is binary
  if(LT == TRUE) {
    colnames(data)[1:3] <- c("time1","time2","status")
    survObj <<- with(data, Surv(time=time1, time2=time2, event=status))
    if(max(data$time2) > 1000){
      timeType = "Days"
      intercept = 5*365}
    else{
      timeType = "Months"
      intercept = 5*12}
  }
  if(LT == FALSE) {
    colnames(data)[1:2] <- c("time","status")
    survObj <<- with(data, Surv(time=time, event=status))
    if(max(data$time) > 1000){
      timeType = "Days"
      intercept = 5*365}
    else{
      timeType = "Months"
      intercept = 5*12}
  }
  ## for 4 levels
  # make KM for Lasso 1 :
  if(numGroups == 2){data$lvl4Groups <- factor(data$lvl4Groups, levels =c("Low","High") )}
  if(numGroups == 3){data$lvl4Groups <- factor(data$lvl4Groups, levels =c("Low","Intermediate","High") )}
  if(numGroups == 4){data$lvl4Groups <- factor(data$lvl4Groups, levels =c("Low","Low-intermediate","High-intermediate","High") )}
  data$lvl4Groups <- as.factor(data$lvl4Groups)
  RiskGroup <<- as.factor(data$lvl4Groups)
  #km.lvl4 <<- survfit(survObj ~ RiskGroup,data=data, conf.type = "log-log")
  ### Get pvalue
  fit0 <- coxph(survObj ~ RiskGroup,data=data,
                na.action=na.exclude)
  log.test.pval <- as.vector(summary(fit0)[10][[1]])[3]
  CI <- as.numeric(as.vector(summary(fit0)[14])[[1]][1])
  if(LT) limit <- as.numeric(quantile(data$time2,0.95))
  if(!LT) limit <- as.numeric(quantile(data$time,0.95))
  KM_2LVLS <- ggsurvplot(survfit(survObj ~ RiskGroup,data=data, conf.type = "log-log"),conf.int  = TRUE,surv.median.line = "hv",
                         data = data,xlim=c(0,limit),break.time.by = 6) + xlab("Time (Months)") +
    labs(title = paste("Kaplan Meier Plot (p-value : " ,round(log.test.pval,digits =4)," and CI : ",round(CI,digits=4), ")",sep=""))+
    geom_vline(xintercept=intercept,col="red", lty = 2)


  ### MAKE TABLE
  if(numGroups == 2) {Groups <- c("Low","High")}
  if(numGroups == 3) {Groups <- c("Low","Intermediate","High")}
  if(numGroups == 4) {Groups <- c("Low","Low-intermediate","High-intermediate","High")}
  survivalGroup <- as.data.frame(matrix(nrow=length(Groups),ncol=4))
  rownames(survivalGroup) <- Groups
  colnames(survivalGroup) <- c("MedianOS","95%CI","1Ysurvival","3Ysurvival")
  # for each group find closest value to median
  if(timeType == "Months"){YR1 <- 1*12;YR3 <- 3*12}
  if(timeType == "Days"){YR1 <- 1*365;YR3 <- 3*365}
  for(i in 1:length(Groups)){
    if(LT == TRUE){NewObject <- with(data[data$lvl4Groups == Groups[i],],Surv(time1,time2,status))}
    if(LT == FALSE){NewObject <- with(data[data$lvl4Groups == Groups[i],],Surv(time,status))}
    Fit <- survfit(NewObject ~ 1,data=data[data$lvl4Groups == Groups[i],], conf.type = "log-log")
    med.index <- which.min(abs(Fit$surv-0.5))
    YR3.index <- which.min(abs(Fit$time-YR1))
    YR5.index <- which.min(abs(Fit$time-YR3))
    survivalGroup[i,] <- c(round(Fit$time[med.index],digits=2),paste0("(",round(Fit$time[which.min(abs(Fit$lower-0.5))],digits=2),",",
                                                                      round(Fit$time[which.min(abs(Fit$upper-0.5))],digits=2),")"),
                           round(Fit$surv[YR3.index],digits=2),round(Fit$surv[YR5.index],digits=2))
  }

  if(mut.data){

    ### MAKE CORRESPONDING MUTATION PER GROUP PLOT ###

    mutDistrib <- as.data.frame(matrix(nrow = numGroups,ncol = length(topHits)))
    rownames(mutDistrib) <- Groups
    colnames(mutDistrib) <- topHits
    for( gene in 1:length(topHits)){
      if(length(data[data$lvl4Groups == "Low",match(topHits[gene],colnames(data))]) != 0 ){
        mutDistrib[match("Low",rownames(mutDistrib)),gene] <- sum(data[data$lvl4Groups == "Low",match(topHits[gene],colnames(data))])/length(data[data$lvl4Groups == "Low",match(topHits[gene],colnames(data))])
      }
      if(length(data[data$lvl4Groups == "Low-intermediate",match(topHits[gene],colnames(data))]) != 0 ){
        mutDistrib[match("Low-intermediate",rownames(mutDistrib)),gene] <- sum(data[data$lvl4Groups == "Low-intermediate",match(topHits[gene],colnames(data))])/length(data[data$lvl4Groups == "Low-intermediate",match(topHits[gene],colnames(data))])
      }
      if(length(data[data$lvl4Groups == "Intermediate",match(topHits[gene],colnames(data))]) != 0 ){
        mutDistrib[match("Intermediate",rownames(mutDistrib)),gene] <- sum(data[data$lvl4Groups == "Intermediate",match(topHits[gene],colnames(data))])/length(data[data$lvl4Groups == "Intermediate",match(topHits[gene],colnames(data))])
      }
      if(length(data[data$lvl4Groups == "High-intermediate",match(topHits[gene],colnames(data))]) != 0  ){
        mutDistrib[match("High-intermediate",rownames(mutDistrib)),gene] <- sum(data[data$lvl4Groups == "High-intermediate",match(topHits[gene],colnames(data))])/length(data[data$lvl4Groups == "High-intermediate",match(topHits[gene],colnames(data))])
      }
      if(length(data[data$lvl4Groups == "High",match(topHits[gene],colnames(data))]) != 0 ){
        mutDistrib[match("High",rownames(mutDistrib)),gene] <- sum(data[data$lvl4Groups == "High",match(topHits[gene],colnames(data))])/length(data[data$lvl4Groups == "High",match(topHits[gene],colnames(data))])
      }
    }
    mutDistrib$Risk <- Groups
    melted.prop <- melt(mutDistrib)
    colnames(melted.prop) <- c("Risk","Gene","Proportion")
    #melted.prop$Risk <- as.factor(melted.prop$Risk)#c("Low","High")
    melted.prop$Risk <- factor(melted.prop$Risk, levels =Groups )
    mut_2LVLS <- ggplot(melted.prop, aes(x = Gene, y = Proportion, fill = Risk)) +
      geom_bar(stat = "identity", position=position_dodge()) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(title = "Proportion of mutations for each genes per risk and method", subtitle = "Using genetic data only")

    ### PIE CHARTS ###
    # count.dups <- function(DF){
    #   DT <- data.table(DF)
    #   DT[,.N, by = names(DT)]
    # }

    ### FIT topHits ###
    if(LT) {
      fit.data <- data[,match(c("time1","time2","status",topHits[1:5]),colnames(data))]
      fit.topHits <- coxph(Surv(time1,time2,status)~.,data= fit.data)}
    if(!LT) {
      fit.data <- data[,match(c("time","status",topHits[1:5]),colnames(data))]
      fit.topHits <- coxph(Surv(time,status)~.,data= fit.data)
    }

    if( LT && max(data$time2) > 1000){MD = 365
    time.type = "Days"}
    if( LT && max(data$time2) < 1000){MD = 12
    time.type = "Months"}
    if( !LT && max(data$time) > 1000){MD = 365
    time.type = "Days"}
    if( !LT && max(data$time) < 1000){MD = 12
    time.type = "Months"}

    ####
    if(length(geneList) == 0) {useGenes <- topHits[1:5]}
    else{
      geneList <- gsub(" ","",geneList)
      useGenes <- geneList}
    pie.data <- data[,match(useGenes,colnames(data))]
    Groups <- as.character(data[,match("lvl4Groups",colnames(data))])
    profile <- apply(pie.data,1,function(x){
      return(paste0(paste0(names(x),"=",as.numeric(x),collapse =",")))
    })
    make.pie.data <- as.data.frame(cbind(Groups,profile))

    ## Save top 2 scenarios for KM
    top.scenarios <- as.data.frame(matrix(nrow = 8,ncol = 5))
    colnames(top.scenarios) <- useGenes
    start <-1

    # high
    high.pie <- try(filter(make.pie.data,Groups == "High"))
    high.dups <- ddply(high.pie,.(Groups,profile),nrow)
    colnames(high.dups)[ncol(high.dups)] <- "N"
    profiles.high <- cbind(high.dups[order(-high.dups$N),],
                           paste("Profile",1:nrow(high.dups)))
    colnames(profiles.high) <- c("Groups","profile","N","Profile")
    scenarios <- lapply(1:nrow(profiles.high),function(x){
      geneList <- as.numeric(gsub(".*=","",unlist(strsplit(as.character(profiles.high$profile[x]),split=","))))
    })
    profiles.high$profile <- unlist(lapply(1:nrow(profiles.high),function(x){
      temp <- unlist(strsplit(as.character(profiles.high$profile[x]),split = ","))
      prof.temp <- temp[grep("=1",temp)]
      prof <- toString(gsub("=1","",prof.temp))
    }) )
    scenarios <- as.data.frame(do.call(rbind, scenarios))
    colnames(scenarios) <- useGenes
    #Prob3YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (3*MD))$surv)
    #Prob5YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (5*MD))$surv)
    #profiles.high <- profiles.high #cbind(profiles.high,Prob3YSurvival,Prob5YSurvival)
    profiles.high$profile[which(profiles.high$profile =="")] <- "No mutants"
    top.scenarios[start:(start+2),] <- scenarios[1:3,]
    rownames(top.scenarios)[start:(start+2)] <- c("High : Profile 1","High : Profile 2","High : Profile 3")
    start <- start+3
    ## Inter high
    if(numGroups==4){
      high.inter.pie <- try(filter(make.pie.data,Groups == "High-intermediate"))
      high.inter.dups <- ddply(high.inter.pie,.(Groups,profile),nrow)
      colnames(high.inter.dups)[ncol(high.inter.dups)] <- "N"
      profiles.inter.high <- cbind(high.inter.dups[order(-high.inter.dups$N),],
                                   paste("Profile",1:nrow(high.inter.dups)))
      colnames(profiles.inter.high) <- c("Groups","profile","N","Profile")
      scenarios <- lapply(1:nrow(profiles.inter.high),function(x){
        geneList <- as.numeric(gsub(".*=","",unlist(strsplit(as.character(profiles.inter.high$profile[x]),split=","))))
      })
      profiles.inter.high$profile <- unlist(lapply(1:nrow(profiles.inter.high),function(x){
        temp <- unlist(strsplit(as.character(profiles.inter.high$profile[x]),split = ","))
        prof.temp <- temp[grep("=1",temp)]
        prof <- toString(gsub("=1","",prof.temp))
      }) )
      scenarios <- as.data.frame(do.call(rbind, scenarios))
      colnames(scenarios) <- useGenes
      #Prob3YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (3*MD))$surv)
      #Prob5YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (5*MD))$surv)
      #profiles.inter.high <- profiles.inter.high #cbind(profiles.inter.high,Prob3YSurvival,Prob5YSurvival)
      profiles.inter.high$profile[which(profiles.inter.high$profile =="")] <- "No mutants"
      top.scenarios[start:(start+2),] <- scenarios[1:3,]
      rownames(top.scenarios)[start:(start+2)] <- c("High-Intermediate : Profile 1","High-Intermediate : Profile 2","High-Intermediate : Profile 3")
      start <- start+3}

    ## INTERMEDIATE
    if(numGroups==3){
      inter.pie <- try(filter(make.pie.data,Groups == "Intermediate"))
      inter.dups <- ddply(inter.pie,.(Groups,profile),nrow)
      colnames(inter.dups)[ncol(inter.dups)] <- "N"
      profiles.inter <- cbind(inter.dups[order(-inter.dups$N),],
                              paste("Profile",1:nrow(inter.dups)))
      scenarios <- lapply(1:nrow(profiles.inter),function(x){
        geneList <- as.numeric(gsub(".*=","",unlist(strsplit(as.character(profiles.inter$profile[x]),split=","))))
      })
      profiles.inter$profile <- unlist(lapply(1:nrow(profiles.inter),function(x){
        temp <- unlist(strsplit(as.character(profiles.inter$profile[x]),split = ","))
        prof.temp <- temp[grep("=1",temp)]
        prof <- toString(gsub("=1","",prof.temp))
      }) )
      scenarios <- as.data.frame(do.call(rbind, scenarios))
      colnames(scenarios) <- useGenes
      #Prob3YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (3*MD))$surv)
      #Prob5YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (5*MD))$surv)
      #profiles.inter <- profiles.inter #cbind(profiles.inter,Prob3YSurvival,Prob5YSurvival)
      profiles.inter$profile[which(profiles.inter$profile =="")] <- "No mutants"
      top.scenarios[start:(start+2),] <- scenarios[1:3,]
      rownames(top.scenarios)[start:(start+2)] <- c("Intermediate : Profile 1","Intermediate : Profile 2","Intermediate : Profile 3")
      start <- start+3}


    ## LOW INTER
    if(numGroups==4){
      low.inter.pie <- try(filter(make.pie.data,Groups == "Low-intermediate"))
      low.inter.dups <- ddply(low.inter.pie,.(Groups,profile),nrow)
      colnames(low.inter.dups)[ncol(low.inter.dups)] <- "N"
      profiles.inter.low <- cbind(low.inter.dups[order(-low.inter.dups$N),],
                                  paste("Profile",1:nrow(low.inter.dups)))
      colnames(profiles.inter.low) <- c("Groups","profile","N","Profile")
      scenarios <- lapply(1:nrow(profiles.inter.low),function(x){
        geneList <- as.numeric(gsub(".*=","",unlist(strsplit(as.character(profiles.inter.low$profile[x]),split=","))))
      })
      profiles.inter.low$profile <- unlist(lapply(1:nrow(profiles.inter.low),function(x){
        temp <- unlist(strsplit(as.character(profiles.inter.low$profile[x]),split = ","))
        prof.temp <- temp[grep("=1",temp)]
        prof <- toString(gsub("=1","",prof.temp))
      }) )
      scenarios <- as.data.frame(do.call(rbind, scenarios))
      colnames(scenarios) <- useGenes
      #Prob3YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (3*MD))$surv)
      #Prob5YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (5*MD))$surv)
      profiles.inter.low <- profiles.inter.low #cbind(profiles.inter.low,Prob3YSurvival,Prob5YSurvival)
      profiles.inter.low$profile[which(profiles.inter.low$profile =="")] <- "No mutants"
      top.scenarios[start:(start+2),] <- scenarios[1:3,]
      rownames(top.scenarios)[start:(start+2)] <- c("Low-Intermediate : Profile 1","Low-Intermediate : Profile 2","Low-Intermediate : Profile 3")
      start <- start+3
    }

    ## LOW
    low.pie <- try(filter(make.pie.data,Groups == "Low"))
    low.dups <- ddply(low.pie,.(Groups,profile),nrow)
    colnames(low.dups)[ncol(low.dups)] <- "N"
    profiles.low <- cbind(low.dups[order(-low.dups$N),],
                          paste("Profile",1:nrow(low.dups)))
    colnames(profiles.low) <- c("Groups","profile","N","Profile")
    scenarios <- lapply(1:nrow(profiles.low),function(x){
      geneList <- as.numeric(gsub(".*=","",unlist(strsplit(as.character(profiles.low$profile[x]),split=","))))
    })
    profiles.low$profile <- unlist(lapply(1:nrow(profiles.low),function(x){
      temp <- unlist(strsplit(as.character(profiles.low$profile[x]),split = ","))
      prof.temp <- temp[grep("=1",temp)]
      prof <- toString(gsub("=1","",prof.temp))
    }) )
    scenarios <- as.data.frame(do.call(rbind, scenarios))
    colnames(scenarios) <- useGenes
    #Prob3YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (3*MD))$surv)
    #Prob5YSurvival <- as.numeric(summary(survfit(fit.topHits, newdata = scenarios, se.fit = F, conf.int = F), times = (5*MD))$surv)
    #profiles.low <- cbind(profiles.low,Prob3YSurvival,Prob5YSurvival)
    profiles.low$profile[which(profiles.low$profile =="")] <- "No mutants"
    top.scenarios[start:(start+2),] <- scenarios[1:3,]
    rownames(top.scenarios)[start:(start+2)] <- c("Low : Profile 1","Low : Profile 2","Low : Profile 3")
    start <- start+3



    if(numGroups == 2){
      pie.chart <- plot_ly() %>%
        add_pie(data = profiles.low, labels = ~profile, values = ~N,
                name = "Low Group", domain = list(x = c(0, 0.4), y = c(0.3, 0.8)),
                hoverinfo="hovertext",
                hovertext = ~paste('',Profile
                                   #'</br> 3 year survival',round(Prob3YSurvival,digits = 2),
                                   #'</br> 5 year survival',round(Prob5YSurvival,digits = 2)
                ),
                type="pie") %>%
        add_pie(data = profiles.high, labels = ~profile, values = ~N,
                name = "High Group", domain = list(x = c(0.45, 0.95), y = c(0.3, 0.8)),hoverinfo="hovertext",
                hovertext = ~paste('',Profile
                                   #'</br> 3 year survival',round(Prob3YSurvival,digits = 2),
                                   #'</br> 5 year survival',round(Prob5YSurvival,digits = 2)
                ),
                type="pie") %>%
        layout(title = "Genetic Profile Distribution", showlegend = F,
               xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
    }

    if(numGroups == 3){
      pie.chart <- plot_ly() %>%
        add_pie(data = profiles.low, labels = ~profile, values = ~N,
                name = "Low Group", domain = list(x = c(0, 0.4), y = c(0.5, 1)),hoverinfo="hovertext",
                hovertext = ~paste('',Profile
                                   #'</br> 3 year survival',round(Prob3YSurvival,digits = 2),
                                   #'</br> 5 year survival',round(Prob5YSurvival,digits = 2)
                ),
                type="pie") %>%
        add_pie(data = profiles.inter, labels = ~profile, values = ~N,
                name = "Intermediate Group", domain = list(x = c(0.45, 0.95), y = c(0.5, 1)),hoverinfo="hovertext",
                hovertext = ~paste('',Profile
                                   #'</br> 3 year survival',round(Prob3YSurvival,digits = 2),
                                   #'</br> 5 year survival',round(Prob5YSurvival,digits = 2)
                ),
                type="pie") %>%
        add_pie(data = profiles.high, labels = ~profile, values = ~N,
                name = "High Group", domain = list(x = c(0.2, 0.6), y = c(0, 0.5)),hoverinfo="hovertext",
                hovertext = ~paste('',Profile
                                   #'</br> 3 year survival',round(Prob3YSurvival,digits = 2),
                                   #'</br> 5 year survival',round(Prob5YSurvival,digits = 2)
                ),
                type="pie") %>%
        layout(title = "Genetic Profile Distribution", showlegend = F,
               xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
    }

    if(numGroups == 4){
      pie.chart <- plot_ly() %>%
        add_pie(data = profiles.low, labels = ~profile, values = ~N,
                name = "Lo", domain = list(x = c(0, 0.22), y = c(0, 0.9)),hoverinfo="hovertext",
                hovertext = ~paste('',Profile
                                   #'</br> Mutant genes',profile
                                   #'</br> 3 year survival',round(Prob3YSurvival,digits = 2),
                                   #'</br> 5 year survival',round(Prob5YSurvival,digits = 2)
                ),
                type="pie") %>%
        add_pie(data = profiles.inter.low, labels = ~profile, values = ~N,
                name = "ILo", domain = list(x = c(0.26, 0.48), y = c(0, 0.9)),hoverinfo="hovertext",
                hovertext = ~paste('',Profile
                                   #'</br> Mutant genes',profile
                                   #'</br> 3 year survival',round(Prob3YSurvival,digits = 2),
                                   #'</br> 5 year survival',round(Prob5YSurvival,digits = 2)
                ),
                type="pie") %>%
        add_pie(data = profiles.inter.high, labels = ~profile, values = ~N,
                name = "IHi", domain = list(x = c(0.52, 0.74), y = c(0, 0.9)),hoverinfo="hovertext",
                hovertext = ~paste('',Profile
                                   #'</br> Mutant genes',profile
                                   #'</br> 3 year survival',round(Prob3YSurvival,digits = 2),
                                   #'</br> 5 year survival',round(Prob5YSurvival,digits = 2)
                ),
                type="pie") %>%
        add_pie(data = profiles.high, labels = ~profile, values = ~N,
                name = "Hi", domain = list(x = c(0.78, 1), y = c(0, 0.9)),hoverinfo="hovertext",
                hovertext = ~paste('',Profile
                                   #'</br> Mutant genes',profile
                                   #'</br> 3 year survival',round(Prob3YSurvival,digits = 2),
                                   #'</br> 5 year survival',round(Prob5YSurvival,digits = 2)
                ),
                type="pie") %>%
        add_annotations(x= 0.05, y= 1, xref = "paper", yref = "paper", text = "<b>Low risk</b>", showarrow = F) %>%
        add_annotations(x= 0.3, y= 1, xref = "paper", yref = "paper", text = "<b>Intermediate low</b>", showarrow = F) %>%
        add_annotations(x= 0.61, y= 1, xref = "paper", yref = "paper", text = "<b>Intermediate high</b>", showarrow = F) %>%
        add_annotations(x= 0.9, y= 1, xref = "paper", yref = "paper", text = "<b>High risk</b>", showarrow = F) %>%
        layout(title = "Genetic Profile Distribution", showlegend = F,
               xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
    }
    useGenes <- paste0(useGenes,collapse = ",")

  }
  #####################################
  else{
    pie.chart <- NULL
    mut_2LVLS <- NULL
    useGenes <- NULL
  }

  return(list("ciSummary" = CI.BP,"inflPlot" = influencePlot,"topHits" = topHits,"average.risk"=average.risk,"data.out"= data,
              "selectInflPlot" = selectInflPlot,"RiskRefit"=refit.risk,
              "RiskHistogram"=RiskHistogram,"Fits"=allCoefs,"time.type"=time.type,
              "RiskScoreSummary"=as.data.frame(t(summary.RiskScore)),"KM_Plot" = KM_2LVLS , "mut_Plot" = mut_2LVLS,
              "SurvSum" = survivalGroup,"PieChart"=pie.chart,"GenesUsed"=useGenes))

}

