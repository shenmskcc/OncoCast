#' riskStrat
#'
#' The riskStrat function let's the user explore the possible stratification created by the risk score
#' attributed to patients by the getResults_OC function.
#' @param data The data used to perform the cross-validated statistical learning run.
#' @param average.risk The average risk score assigned to each patient. This vector can be retrieved
#' by using the outputSummary's function output.
#' @param numGroups The number of groups to be made when stratifying by risk groups. Options are 2,3 and 4 (for now implementing
#' broader version). Default is 2.
#' @param cuts Numeric vector of the points in the distribution of risk scores where groups will be splitted. Needs to be of length
#' numGroups - 1. eg : c(0.25,0.5,0.75) when numgroups is 4. Default is 0.5.
#'
#' @return KM_Plot : A Kaplan-Meier plot with patients stratified by the cuts made in the risk score.
#' @return SurvSum : A summary table of the survival for each group generated
#' @return data.out : The data used to generate the plots with an extra column giving the respective groups of patients.
#'
#' @keywords Results, stratification, survival
#' @export
#' @examples library(OncoCast)
#' test <- OncoCast(data=survData,formula=Surv(time,status)~.,
#'                           method=c("LASSO"),
#'                           runs = 25,cores = 1,sampling = "cv",
#'                           pathResults = "./Test/",studyType = "ExampleRun",save=F)
#' out <- outputSummary(test$LASSO)
#' riskout <- riskStrat(survData,out$average.risk,numGroups = 2,cuts=0.5)
#' @import survival
#' @import ggplot2
#' @import plotly
#' @import plyr
#' @import stats
#' @import reshape2
#' @import scales
#' @import survminer
#' @import data.table

riskStrat <- function(data,average.risk,numGroups,cuts){


  MD <- 12
  if(length(grep("time",colnames(data)))  == 1) {LT = FALSE}
  if(length(grep("time",colnames(data)))  == 2) {LT = TRUE}
  data$RiskGroup <- rep(NA,nrow(data))
  qts <- as.numeric(quantile(average.risk,cuts,na.rm = T))

  if(numGroups ==2){
    for(i in 1:nrow(data)){

      if(is.character(try(average.risk[i] < qts,silent=T))){
        stop("ERROR : Increase the number of runs. Some patients were never attributed to the testing set. Recommended value is above 20")
      }

      if(average.risk[i] < qts) {
        data$RiskGroup[i] <- "Low"
      }

      if(average.risk[i] >= qts) {
        data$RiskGroup[i] <- "High"
      }
    }
  }

  if(numGroups == 3){

    for(i in 1:nrow(data)){

      if(is.character(try(average.risk[i] < qts,silent=T))){
        stop("ERROR : Increase the number of runs. Some patients were never attributed to the testing set. Recommended value is above 20")
      }

      if(average.risk[i] < qts[1]) {
        data$RiskGroup[i] <- "Low"
      }
      if(average.risk[i] >= qts[1] && average.risk[i] < qts[2]){
        data$RiskGroup[i] <- "Intermediate"
      }
      if(average.risk[i] >= qts[2]) {
        data$RiskGroup[i] <- "High"
      }
    }
  }

  if(numGroups == 4){
    for(i in 1:nrow(data)){

      if(is.character(try(average.risk[i] < qts,silent=T))){
        stop("ERROR : Increase the number of runs. Some patients were never attributed to the testing set. Recommended value is above 20")
      }

      if(average.risk[i] < qts[1]) {
        data$RiskGroup[i] <- "Low"
      }
      if(average.risk[i] >= qts[1] && average.risk[i] < qts[2]){
        data$RiskGroup[i] <- "Low-intermediate"
      }
      if(average.risk[i] >= qts[2] && average.risk[i] < qts[3]){
        data$RiskGroup[i] <- "High-intermediate"
      }
      if(average.risk[i] >= qts[3]) {
        data$RiskGroup[i] <- "High"
      }
    }
  }


  ##################
  if(numGroups == 2){data$RiskGroup <- factor(data$RiskGroup, levels =c("Low","High") )}
  if(numGroups == 3){data$RiskGroup <- factor(data$RiskGroup, levels =c("Low","Intermediate","High") )}
  if(numGroups == 4){data$RiskGroup <- factor(data$RiskGroup, levels =c("Low","Low-intermediate","High-intermediate","High") )}
  data$RiskGroup <- as.factor(data$RiskGroup)

  ### Is this left truncated ?
  # Will be LT if second column is binary
  if(LT == TRUE) {
    fit0 <- coxph(Surv(time1,time2,status) ~ RiskGroup,data=data,
                  na.action=na.exclude)
    if(max(data$time2) > 1000){
      timeType = "Days"
      intercept = 5*365}
    else{
      timeType = "Months"
      intercept = 5*12}
  }
  if(LT == FALSE) {
    fit0 <- coxph(Surv(time,status) ~ RiskGroup,data=data,
                  na.action=na.exclude)
    if(max(data$time) > 1000){
      timeType = "Days"
      intercept = 5*365}
    else{
      timeType = "Months"
      intercept = 5*12}
  }

  log.test.pval <- as.vector(summary(fit0)[10][[1]])[3]
  CI <- as.numeric(as.vector(summary(fit0)[14])[[1]][1])
  if(LT) limit <- as.numeric(quantile(data$time2,1))
  if(!LT) limit <- as.numeric(quantile(data$time,1))

  if(LT) {KM_2LVLS <- ggsurvplot(survfit(Surv(time1,time2,status) ~ RiskGroup,data=data, conf.type = "log-log"),conf.int  = TRUE,surv.median.line = "hv",
                         data = data,xlim=c(0,limit),break.time.by = 6) + xlab("Time (Months)") +
    labs(title = paste("Kaplan Meier Plot (p-value : " ,round(log.test.pval,digits =4),")",sep=""))}
  if(!LT){KM_2LVLS <- ggsurvplot(survfit(Surv(time,status) ~ RiskGroup,data=data, conf.type = "log-log"),conf.int  = TRUE,surv.median.line = "hv",
                                 data = data,xlim=c(0,limit),break.time.by = 6) + xlab("Time (Months)") +
    labs(title = paste("Kaplan Meier Plot (p-value : " ,round(log.test.pval,digits =4),")",sep=""))}


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
    if(LT == TRUE){NewObject <- with(data[data$RiskGroup == Groups[i],],Surv(time1,time2,status))}
    if(LT == FALSE){NewObject <- with(data[data$RiskGroup == Groups[i],],Surv(time,status))}
    Fit <- survfit(NewObject ~ 1,data=data[data$RiskGroup == Groups[i],], conf.type = "log-log")
    # med.index <- which.min(abs(Fit$surv-0.5))
    YR3.index <- which.min(abs(Fit$time-YR1))
    YR5.index <- which.min(abs(Fit$time-YR3))
    survivalGroup[i,] <- c(as.numeric(round(summary(Fit)$table[7],digits=2)),
                           paste0("(",as.numeric(round(summary(Fit)$table[8],digits=2)),",",
                                  as.numeric(round(summary(Fit)$table[9],digits=2)),")"),
                           round(Fit$surv[YR3.index],digits=2),round(Fit$surv[YR5.index],digits=2))
  }

return(list("KM_Plot" = KM_2LVLS,"SurvSum" = survivalGroup,"data.out"=data))

}
