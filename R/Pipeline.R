##########################################################################################################
##################### FULL VARIABLE SELECTION ANALYSIS FOR IMPACT GENETIC DATA ###########################
##########################################################################################################

##########################################################################################################
##################### DATA :
### This function takes as input data a processed data frame containing :
# 1) IMPACT mutation data in binary form (mutated = 1, non-mutated = 0)
# 2) Survival Times and outcome status
# 3) No missing data is allowed (maybe allow for imputation later)

##################### ARGUMENTS :
# 1) data (specified in the way above)
# 2) A survival formula (eg : Surv(time,status)~. or Surv(time1,time2,status)~.)
# 3) Variable selection method to be used possible are : "LASSO","RF","GBM" --> default all
# 4) In the case where "RF or "GBM" are chosen can specify :
# i) number of trees --> default 1000
# ii) nsplits ("RF" only) --> default 10
# 5) Number of runs to be performed (number of splits) --> default 100
# 6) Where to save the results directory --> default working directory
# 7) Name of the output (begining) --> default none
# 8) number of cores to be used "cores" --> default 1
# 9) bootControl number of CV for the GBM parameters --> default 50

#################### NOTES :
# 1) In the case where interaction depth of shrinkage is in fact a vector
# Cross validation will be automatically performed with 50 folds


#' OncoCast
#'
#' This functions let's the user select one or multiple machine learning algorithms. This is intended for survival data, and
#' the mehods can handle left truncated survival data. The inputed data should be a data frame with columns representing the
#' variables of interest while the rows should correspond to patients.
#' The output will vary depending on the methods. However they will all return a list of length equal to the number of
#' crossvalidation performed, including the predicted risk score at each cross validation for all the patients falling in the
#' test set.
#' @param data name of the data frame with variables as columns and patients as rows.
#' @param formula Must match the family chosen in the previous parameter. For the "gaussian" and "binomial" families the
#' formula should be of the form 'y~.'. For the 'cox' family, a survival formula with the names of the variables to be used in the data frame provided in the first argument.
#'  eg : Surv(time,status)~. or Surv(time1,time2,status)~. (Note all the variable available will be used regardless of the right
#'  side of the formula).
#' @param method Character vector of the names of the method(s) to be used, options are : 1) LASSO ("LASSO") 2) Ridge ("RIDGE")
#'  3) Elastic Net ("ENET"). Default is all.
#' @param runs Number of cross validation iterations to be performed. Default is 100.
#' @param sampling The method use for sampling, options are bootstrapping ("boot") and cross-validation ("cv").
#' Default is cross-validation.
#' @param penalizeCol Name of variables you do not with to penalize (available only for LASSO, RIDGE and ENET methods). Default is NULL.
#' @param cores If you wish to run this function in parallel set the number of cores to be used to be greater than 1. Default is 1.
#' CAUTION : Overloading your computer can lead to serious issues, please check how many cores are available on your machine
#' before selecting an option !
#' @param pathResults String of where you wish to output the results. Default is current directory.
#' @param studyType String that will be the prefix to the name of the outputed results. Default is empty.
#' @param save Boolean value : Default is TRUE, the results will be saved with the specified name in the specified path. If FALSE the results
#' will be returned directly from the function and won't be saved.
#' @return CI : For each iteration the concordance index of the generated model will be calculated on the testing set
#' @return predicted : The predicted risk for all the patients in the testing set
#' @return fit : For LASSO, RIDGE and ENET methods return the refitted cox proportional hazard model with the beta coefficients found
#' in the penalized regression model.
#' @return Selection : The number of times each variables is selected accross all the trees. For the CF method only.
#' @keywords Selection, penalized regression
#' @export
#' @examples library(OncoCast)
#' test <- OncoCast(data=survData,formula=Surv(time,status)~.,
#'                           method=c("LASSO"),
#'                           runs = 5,cores = 1,sampling = "cv",
#'                           pathResults = "./Test/",studyType = "ExampleRun",save=F)
#' @import survival
#' @import penalized
#' @import foreach
#' @import doParallel
#' @import party
#' @import stats


OncoCast <- function(data,formula, method = c("LASSO","RIDGE","ENET"),
                               runs = 100,penalizeCol = NULL,
                               sampling = "cv",
                               cores = 1,
                               pathResults = "",studyType = "",save = T){

  # load libraries
  library(foreach)
  library(doParallel)

  # Missingness
  if(anyNA(data)){
    stop("ERROR : Missing data is not allowed at this time, please remove or impute missing data.")
  }

  if( anyNA(match(method,c("LASSO","ENET","RIDGE"))) ){stop("ERROR : The method you have selected is not available.")}

  if(!(sampling %in% c("cv","boot"))){
    stop("ERROR : This sampling method is not available. OPTIONS : 'cv' or 'boot'.")
  }


  ### Prepare cores
  cl <- makeCluster(cores) #not to overload your computer
  registerDoParallel(cl)

  # appropriate formula
  survFormula <- as.formula(formula)
  survResponse <- survFormula[[2]]

  if(!length(as.list(survResponse)) %in% c(3,4)){
    stop("ERROR : Response must be a 'survival' object with 'Surv(time, event)' or 'Surv(time1, time2, event)'.")
  }
  ### reprocess data
  if(length(as.list(survResponse)) == 3){
    colnames(data)[match(as.list(survResponse)[2:3],colnames(data))] <- c("time","status")
    LT = FALSE
  }
  if(length(as.list(survResponse)) == 4){
    colnames(data)[match(as.list(survResponse)[2:4],colnames(data))] <- c("time1","time2","status")
    LT = TRUE
  }


  print("Data check performed, ready for analysis.")

  ##### SURVIVAL (INCLUDING LEFT TRUNCATION) #####
  LASSO <- NULL
  ENET <- NULL
  RIDGE <- NULL


  ########## LASSO #############
  # returns and saves the CI and the refitted model
  final.lasso <- list()
  if("LASSO" %in% method) {
    print("LASSO SELECTED")
    LASSO <- foreach(run=1:runs) %dopar% {

      #load into for each
      library(survival)
      #library(plyr, lib.loc="/usr/local/lib/R/site-library/")
      try(library(penalized, lib.loc="/usr/local/lib/R/site-library/"),silent = TRUE)
      try(library(penalized),silent = TRUE)
      library(foreach)
      library(doParallel)
      library(stats)

      ### BUILD TRAINING AND TESTING SET ###
      set.seed(run)
      print(paste("Run : ", run,sep=""))
      # split data
      if(sampling == "cv"){
        rm.samples <- sample(1:nrow(data), ceiling(nrow(data)*1/3),replace = FALSE)
        train <- data[-rm.samples,]
        test <- data[rm.samples,]}

      if(sampling == "boot"){
        rm.samples <- sample(1:nrow(data), nrow(data),replace = TRUE)
        train <- data[rm.samples,]
        test <- data[-rm.samples,]
      }


      # make new survival objects
      if(LT) {
        trainSurv <- with(train,Surv(time1,time2,status))
        testSurv <- with(test,Surv(time1,time2,status))
        if(is.null(penalizeCol)){
          opt <- try(optL1(trainSurv,data = train,penalized = train[,4:ncol(train)],
                           unpenalized = ~0,fold = 5,trace=FALSE))
        }

        else{
          noPen.index <- match(penalizeCol,colnames(train))
          noPen.index.pen <-c(1:3,noPen.index)
          opt <- try(optL1(trainSurv,data = train,penalized = train[,-noPen.index.pen],
                           unpenalized = train[,noPen.index],fold = 5,trace=FALSE))}
      }

      else {
        trainSurv <- with(train,Surv(time,status))
        testSurv <- with(test,Surv(time,status))

        if(is.null(penalizeCol)){
          opt <- try(optL1(trainSurv,data = train,penalized = train[,3:ncol(train)],
                           unpenalized = ~0,fold = 5,trace=FALSE))}

        else{
          noPen.index <- match(penalizeCol,colnames(train))
          noPen.index.pen <-c(1,2,noPen.index)
          opt <- try(optL1(trainSurv,data = train,penalized = train[,-noPen.index.pen],
                           unpenalized = train[,noPen.index],fold = 5,trace=FALSE))
        }

      }


      if(typeof(opt) == "list" && length(coefficients(opt$fullfit)) != 0){
        # get optimal coefficients
        optimal.coefs <- coefficients(opt$fullfit)
        coefs.left <- names(coefficients(opt$fullfit))
        if(LT) {lasso.formula <- as.formula(paste("Surv(time1,time2,status) ~ ",paste(coefs.left, collapse= "+")))}
        else {lasso.formula <- as.formula(paste("Surv(time,status) ~ ",paste(coefs.left, collapse= "+")))}

        # Refit model with the corresponding coefficients
        lasso.fit <- coxph(lasso.formula, data=train,init=optimal.coefs,iter=0)

        # save what you need to save
        if(sampling == "cv"){
          CI <- as.numeric(survConcordance(testSurv ~predict(lasso.fit, newdata=test),
                                           test)$concordance)}
        if(sampling=="boot"){
          CI.0632 <- 0.632*as.numeric(survConcordance(testSurv ~predict(lasso.fit, newdata=test),
                                                      test)$concordance) +
            0.368*as.numeric(survConcordance(trainSurv ~predict(lasso.fit, newdata=train),
                                             train)$concordance)}
        CI <- as.numeric(survConcordance(testSurv ~predict(lasso.fit, newdata=test),
                                         test)$concordance)
        fit <- summary(lasso.fit)$coefficients
        predicted <- predict(lasso.fit, newdata=test)
      }
      else{
        CI <- NA
        fit <- NA
        predicted <- NA
        CI.0632 <- NA
        if(run == 1){
          final.lasso$method <- "LASSO"
          final.lasso$data <- data}
      }
      if(run == 1){
        final.lasso$method <- "LASSO"
        final.lasso$data <- data
      }
      final.lasso$CI <- CI
      if(sampling == "boot"){ final.lasso$CI.0632 <- CI.0632 }
      final.lasso$fit <- fit
      final.lasso$predicted <- predicted

      return(final.lasso)
    }

    if(save){
      save(LASSO,file = paste0(pathResults,studyType,"_",sampling,"_LASSO.Rdata"))
      LASSO <- NULL}
  }


  ##### RIDGE #####

  final.ridge <- list()
  if("RIDGE" %in% method) {
    print("RIDGE SELECTED")
    RIDGE <- foreach(run=1:runs) %dopar% {

      #load into for each
      library(survival)
      #library(plyr, lib.loc="/usr/local/lib/R/site-library/")
      try(library(penalized, lib.loc="/usr/local/lib/R/site-library/"),silent = TRUE)
      try(library(penalized),silent = TRUE)
      library(foreach)
      library(doParallel)
      library(stats)

      ### BUILD TRAINING AND TESTING SET ###
      set.seed(run)
      print(paste("Run : ", run,sep=""))
      # split data
      if(sampling == "cv"){
        rm.samples <- sample(1:nrow(data), ceiling(nrow(data)*1/3),replace = FALSE)
        train <- data[-rm.samples,]
        test <- data[rm.samples,]}

      if(sampling == "boot"){
        rm.samples <- sample(1:nrow(data), nrow(data),replace = TRUE)
        train <- data[rm.samples,]
        test <- data[-rm.samples,]
      }


      # make new survival objects
      if(LT) {
        trainSurv <- with(train,Surv(time1,time2,status))
        testSurv <- with(test,Surv(time1,time2,status))
        opt <- try(optL2(trainSurv,data = train,penalized = train[,4:ncol(train)],fold = 5,trace=FALSE))
      }

      else {
        trainSurv <- with(train,Surv(time,status))
        testSurv <- with(test,Surv(time,status))
        opt <- try(optL1(trainSurv,data = train,penalized = train[,3:ncol(train)],fold = 5,trace=FALSE))
      }

      if(typeof(opt) == "list" && length(coefficients(opt$fullfit)) != 0){
        # get optimal coefficients
        optimal.coefs <- coefficients(opt$fullfit)
        coefs.left <- names(coefficients(opt$fullfit))
        if(LT) {ridge.formula <- as.formula(paste("Surv(time1,time2,status) ~ ",paste(coefs.left, collapse= "+")))}
        else {ridge.formula <- as.formula(paste("Surv(time,status) ~ ",paste(coefs.left, collapse= "+")))}

        # Refit model with the corresponding coefficients
        ridge.fit <- coxph(ridge.formula, data=train,init=optimal.coefs,iter=0)

        # save what you need to save
        if(sampling == "cv"){
          CI <- as.numeric(survConcordance(testSurv ~predict(ridge.fit, newdata=test),
                                           test)$concordance)}
        if(sampling=="boot"){
          CI <- 0.632*as.numeric(survConcordance(testSurv ~predict(ridge.fit, newdata=test),
                                                 test)$concordance) +
            0.368*as.numeric(survConcordance(trainSurv ~predict(ridge.fit, newdata=train),
                                             train)$concordance)}
        final.ridge$CI <- CI
        final.ridge$fit <- summary(ridge.fit)$coefficients
        final.ridge$predicted <- predict(ridge.fit, newdata=test)
        if(run == 1){
          final.ridge$method <- "RIDGE"
          final.ridge$data <- data}

        return(final.ridge)
      }

    }
    if(save){
      save(RIDGE,file = paste0(pathResults,studyType,"_",sampling,"_RIDGE.Rdata"))
      RIDGE <- NULL}
  }

  #### LASSO ENET #####

  final.enet <- list()
  if("ENET" %in% method) {
    print("ENET SELECTED")
    ENET <- foreach(run=1:runs) %dopar% {
      library(survival)
      try(library(penalized, lib.loc="/usr/local/lib/R/site-library/"),silent = TRUE)
      try(library(penalized),silent = TRUE)
      library(foreach)
      library(doParallel)
      library(stats)

      ### BUILD TRAINING AND TESTING SET ###
      set.seed(run)
      print(paste("Run : ", run,sep=""))
      # split data
      if(sampling == "cv"){
        rm.samples <- sample(1:nrow(data), ceiling(nrow(data)*1/3),replace = FALSE)
        train <- data[-rm.samples,]
        test <- data[rm.samples,]}

      if(sampling == "boot"){
        rm.samples <- sample(1:nrow(data), nrow(data),replace = TRUE)
        train <- data[rm.samples,]
        test <- data[-unique(rm.samples),]
      }

      if(LT){
        trainSurv <- with(train,Surv(time1,time2,status))
        testSurv <- with(test,Surv(time1,time2,status))}
      if(!LT){
        trainSurv <- with(train,Surv(time,status))
        testSurv <- with(test,Surv(time,status))
      }

      alpha.final <- 0.01
      opt <- optL1(trainSurv,data = train,penalized = train[,4:ncol(train)],
                   unpenalized = ~0,fold=5,lambda2 = alpha.final,model = "cox",trace = FALSE)

      if(typeof(opt) == "list" && length(coefficients(opt$fullfit)) != 0){
        optimal.coefs <- coefficients(opt$fullfit)
        coefs.left <- names(coefficients(opt$fullfit))
        if(LT) lasso.formula <- as.formula(paste("Surv(time1,time2,status) ~ ",paste(coefs.left, collapse= "+")))
        if(!LT) lasso.formula <- as.formula(paste("Surv(time,status) ~ ",paste(coefs.left, collapse= "+")))
        # Refit model with the corresponding coefficients
        lasso.fit <- coxph(lasso.formula, data=train,init=optimal.coefs,iter=0)
        # save what you need to save
        if(sampling == "cv"){
          CI <- as.numeric(survConcordance(testSurv ~predict(lasso.fit, newdata=test),
                                           test)$concordance)}
        if(sampling=="boot"){
          CI <- 0.632*as.numeric(survConcordance(testSurv ~predict(lasso.fit, newdata=test),
                                                 test)$concordance) +
            0.368*as.numeric(survConcordance(trainSurv ~predict(lasso.fit, newdata=train),
                                             train)$concordance)}

        final.enet$CI <- CI
        final.enet$fit <- summary(lasso.fit)$coefficients
        final.enet$predicted <- predict(lasso.fit, newdata=test)
        final.enet$alphas <- alpha.final
        final.enet$data <- NULL
        if(run == 1){
          final.enet$method <- "ENET"
          final.enet$data <- data
        }
        return(final.enet)
      }

    }
    if(save){
      save(ENET,file = paste0(pathResults,studyType,"_",sampling,"_ENET.Rdata"))
      ENET <- NULL}
  }


  ################

  stopCluster(cl)
  OUTPUT <- list()
  if(!is.null(LASSO)){OUTPUT$LASSO <- LASSO}
  if(!is.null(RIDGE)){OUTPUT$RIDGE <- RIDGE}
  if(!is.null(ENET)){OUTPUT$ENET <- ENET}

  if(!is.null(OUTPUT)){return(OUTPUT)}
  else{return(0)}

}


# library(OncoCast)
# test <- OncoCast(data=survData,formula=Surv(time,status)~.,
#                   method=c("LASSO"),
#                   runs = 50,cores = 1,sampling = "cv",
#                   pathResults = "./",studyType = "ExampleRun",save=F)
#
# out <- getResults_OC(test$LASSO,numGroups=2,cuts=0.5,geneList=NULL,mut.data = T)


