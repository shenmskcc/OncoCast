#' outputSummary
#'
#' This functions let's the user study the basic output of the OncoCast function. This function takes as input
#' one of the (or the one) objects returned from the different machine learning algorithms chosen previously.
#' Only one such object can be inputted everytime in the outputSummary function.
#' @param OC_object A list object outputed by the VariableSelection function.
#' @param data A dataframe that corresponds to the data used to generate the OncoCast output.
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
#' @keywords Results
#' @export
#' @examples library(OncoCast)
#' test <- OncoCast(data=survData,formula=Surv(time,status)~.,
#'                           method=c("LASSO"),
#'                           runs = 25,cores = 1,sampling = "cv",
#'                           pathResults = "./Test/",studyType = "ExampleRun",save=F)
#' out <- outputSummary(test$LASSO)
#' @import survival
#' @import ggplot2
#' @import plotly
#' @import plyr
#' @import stats
#' @import reshape2
#' @import scales

outputSummary <- function(OC_object,data){


  OC_object <- Filter(Negate(is.null), OC_object)

  ## determine if left truncated
  if(length(grep("time",colnames(data)))  == 1) {LT = FALSE}
  if(length(grep("time",colnames(data)))  == 2) {LT = TRUE}

  MD <- 12

  ConcordanceIndex <- as.data.frame(as.vector(unlist(sapply(OC_object, "[[", "CI"))))
  summary.CI <- round(as.data.frame(c(quantile(ConcordanceIndex[,1],c(0.1,0.25,0.5,0.75,0.9),na.rm = T))),digits = 2)
  colnames(summary.CI) <- "Concordance Index"
  rownames(summary.CI) <- c("Lower 10%","1st Quarter","Median","3rd Quarter","Upper 10%")
  CI.BP <- as.data.frame(t(summary.CI))

  allCoefs <- t(sapply(OC_object,"[[","fit"))
  allCoefs[is.na(allCoefs)] <- 0

  selected.genes.lasso <- apply(allCoefs,2,function(x){sum(x!=0)})

  #topHits <- names(sort(selected.genes.lasso[selected.genes.lasso >0.5],decreasing = TRUE))
  if(length(selected.genes.lasso)  >= 20){ranked.lasso <- sort(selected.genes.lasso,decreasing = TRUE)[1:20]}
  if(length(selected.genes.lasso)  < 20){ranked.lasso <- sort(selected.genes.lasso,decreasing = TRUE)[1:length(selected.genes.lasso)]}
  melt.rank.lasso <- melt(ranked.lasso)
  melt.rank.lasso$Gene <- factor(rownames(melt.rank.lasso), levels = rownames(melt.rank.lasso))
  melt.rank.lasso$value <- melt.rank.lasso$value/length(OC_object)
  colnames(melt.rank.lasso) <- c("Frequency","Gene")

  influencePlot <- ggplot(melt.rank.lasso,aes(x=Gene,y=Frequency,fill=Gene))+geom_col()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = paste("Most recurrent selected genes out of",length(OC_object),"runs")) +
    theme(legend.position="none")

  if(LT) Variables <- colnames(data)[-c(1:3)]
  if(!LT) Variables <- colnames(data)[-c(1:2)]

  meanCoefs <- apply(allCoefs,2,function(x){mean(x,na.rm = TRUE)})
  selectFreq <- apply(allCoefs,2,function(x){
    length(which(x!=0))/length(x)
  })

  if(length(selectFreq[selectFreq > 0.5]) > 2) topHits <- names(selectFreq[selectFreq > 0.5])
  else topHits <- names(selectFreq[order(selectFreq,decreasing = T)])[1:10]
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
                                          '<br> Hazard Ratio :',round(exp(MeanCoefficient),digits=2)),
                            mode = "markers",size = ~MutationFrequency,color = ~MutationFrequency) %>%
    layout(title ="Volcano Plot")

  # final.pred <- as.data.frame(matrix(nrow= nrow(data),ncol = length(OC_object)))
  # rownames(final.pred) <- rownames(data)
  # for(i in 1:length(OC_object)){
  #   temp <- OC_object[[i]]$predicted
  #   final.pred[match(names(temp),rownames(final.pred)),i] <- as.numeric(temp)
  # }

  final.pred <- sapply(OC_object,"[[","predicted")
  average.risk <- apply(final.pred,1,function(x){
    mean(as.numeric(x),na.rm = TRUE)
  })
  average.risk[which(is.na(average.risk))] <- NA
  to <- c(0,10)
  from <- range(average.risk, na.rm = TRUE, finite = TRUE)
  RiskScore <- (as.numeric(average.risk)-from[1])/diff(from)*diff(to)+to[1]
  #RiskScore <- rescale(as.numeric(average.risk), to = c(0, 10), from = range(average.risk, na.rm = TRUE, finite = TRUE))
  summary.RiskScore <- round(as.data.frame(c(quantile(RiskScore,c(0.1,0.25,0.33,0.5,0.66,0.75,0.9),na.rm = TRUE))),digits = 2)
  colnames(summary.RiskScore) <- "Risk Score"
  rownames(summary.RiskScore) <- c("Lower 10%","1st Quarter","1st Tertile","Median","2nd Tertile","3rd Quarter","Upper 10%")
  ## refit coxph model with average risk as covariate
  meanRS <- mean(RiskScore)
  #RiskScore <- average.risk #- meanRS
  if(LT) refit.risk <- coxph(Surv(data$time1,data$time2,data$status)~RiskScore)
  if(!LT) refit.risk <- coxph(Surv(data$time,data$status)~RiskScore)
  Risk <- as.data.frame(RiskScore)

  RiskHistogram <- ggplot(Risk, aes(x = RiskScore, y = ..density..)) +
    geom_histogram(show.legend = FALSE, aes(fill=..x..),
                   breaks=seq(min(Risk$RiskScore,na.rm = T), max(Risk$RiskScore,na.rm = T), by=0.25)) +
    geom_density(show.legend = FALSE) +
    theme_minimal() +
    labs(x = "Average risk score", y = "Density") +
    scale_fill_gradient(high = "red", low = "green")


  return(list("ciSummary" = CI.BP,"inflPlot" = influencePlot,"topHits" = topHits,"average.risk"=average.risk,
              "selectInflPlot" = selectInflPlot,"RiskRefit"=refit.risk,"scaled.risk"=RiskScore,
              "RiskHistogram"=RiskHistogram,"Fits"=allCoefs,
              "RiskScoreSummary"=as.data.frame(t(summary.RiskScore))))

}
