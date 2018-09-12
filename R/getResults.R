#' getResults_OC
#'
#' This functions let's the user study the output of the OncoCast function. This function takes as input
#' one of the (or the one) objects returned from the different machine learning algorithms chosen previously.
#' Only one such object can be inputted at a time in the getResults_OC function.
#' @param OC_object A list object outputed by the OncoCast function.
#' @param data A dataframe that corresponds to the data used to generate the OncoCast output.
#' @param numGroups The number of groups to be made when stratifying by risk groups. Options are 2,3 and 4 (for now implementing
#' broader version). Default is 2.
#' @param cuts Numeric vector of the points in the distribution of risk scores where groups will be splitted. Needs to be of length
#' numGroups - 1. eg : c(0.25,0.5,0.75) when numgroups is 4. Default is 0.5.
#' @param geneList Optional character vector of gene names to use to generate the pie charts. Default is NULL, which
#' leads to using the 5 most frequently selected features.
#' @param mut.data Boolean argument indicating if the user is using mutation predictors (binary data). Default is FALSE.
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
#'                           method=c("LASSO"),
#'                           runs = 25,cores = 1,sampling = "cv",
#'                           pathResults = "./Test/",studyType = "ExampleRun",save=F)
#' out <- getResults_OC(test$LASSO,numGroups=2,cuts=0.5,geneList=NULL,mut.data=T)
#' @import survival
#' @import ggplot2
#' @import plotly
#' @import plyr
#' @import stats
#' @import reshape2
#' @import scales
#' @import survminer
#' @import data.table
#' @import gplots


getResults_OC <- function(OC_object,data,numGroups=2,cuts=0.5,geneList=NULL,mut.data=F){


  ############## CHECKS
  if(!(numGroups %in% 2:4)) {stop("ERROR : You must use a number of groups between 2 and 4. We are working on implementing a broader version")}
  if(max(cuts) >= 1 || min(cuts) <= 0){stop("ERROR : You must select cuts that are between 0 and 1 (Not included).")}

  OC_object <- Filter(Negate(is.null), OC_object)

  ## determine if left truncated
  if(length(grep("time",colnames(data)))  == 1) {LT = FALSE}
  if(length(grep("time",colnames(data)))  == 2) {LT = TRUE}

  MD <- 12

  #################################################################
  ### get all the basic results from the output of the pipeline ###
  #################################################################

  basic.results <- outputSummary(OC_object,data)

  average.risk <- basic.results$average.risk
  #####################################
  ####### RISK STRATIFICATION #########
  #####################################

  strat.results <- riskStrat(data,average.risk,numGroups,cuts)
  data <- strat.results$data.out
  topHits <- basic.results$topHits

  ######################################
  ####### MUTATION DISTRIBUTION ########
  ######################################
  if(mut.data){
    mut.results <- mutSummary(data,average.risk,topHits,numGroups,geneList)
    pie.chart <- mut.results$PieChart
    mut_2LVLS <- mut.results$mut_Plot
    useGenes <- mut.results$GenesUsed
    #MutProfiles <- mut.results$MutProfiles
  }
  #####################################
  else{
    pie.chart <- NULL
    mut_2LVLS <- NULL
    useGenes <- NULL
    #MutProfiles <- NULL
  }

  return(list("ciSummary" = basic.results$ciSummary,"inflPlot" = basic.results$inflPlot,
              "topHits" = basic.results$topHits,"average.risk"=basic.results$average.risk,
              "selectInflPlot" = basic.results$selectInflPlot,"RiskRefit"=basic.results$RiskRefit,
              "scaled.risk"=basic.results$scaled.risk,
              "RiskHistogram"=basic.results$RiskHistogram,"Fits"=basic.results$Fits,"time.type"=basic.results$time.type,
              "RiskScoreSummary"=basic.results$RiskScoreSummary,
              "KM_Plot" = strat.results$KM_Plot,"SurvSum" = strat.results$SurvSum,
              "mut_Plot" = mut_2LVLS,"PieChart"=pie.chart,"GenesUsed"=useGenes))
              # ,
              # "MutProfiles"=MutProfiles))

}

