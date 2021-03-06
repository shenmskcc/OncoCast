% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/riskStrat.R
\name{riskStrat}
\alias{riskStrat}
\title{riskStrat}
\usage{
riskStrat(data, average.risk, numGroups, cuts)
}
\arguments{
\item{data}{The data used to perform the cross-validated statistical learning run.}

\item{average.risk}{The average risk score assigned to each patient. This vector can be retrieved
by using the outputSummary's function output.}

\item{numGroups}{The number of groups to be made when stratifying by risk groups. Options are 2,3 and 4 (for now implementing
broader version). Default is 2.}

\item{cuts}{Numeric vector of the points in the distribution of risk scores where groups will be splitted. Needs to be of length
numGroups - 1. eg : c(0.25,0.5,0.75) when numgroups is 4. Default is 0.5.}
}
\value{
KM_Plot : A Kaplan-Meier plot with patients stratified by the cuts made in the risk score.

SurvSum : A summary table of the survival for each group generated

data.out : The data used to generate the plots with an extra column giving the respective groups of patients.
}
\description{
The riskStrat function let's the user explore the possible stratification created by the risk score
attributed to patients by the getResults_OC function.
}
\examples{
library(OncoCast)
test <- OncoCast(data=survData,formula=Surv(time,status)~.,
                          method=c("LASSO"),
                          runs = 25,cores = 1,sampling = "cv",
                          pathResults = "./Test/",studyType = "ExampleRun",save=F)
out <- outputSummary(test$LASSO)
riskout <- riskStrat(survData,out$average.risk,numGroups = 2,cuts=0.5)
}
\keyword{Results,}
\keyword{stratification,}
\keyword{survival}
