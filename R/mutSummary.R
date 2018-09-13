#' mutSummary
#'
#' Let's the user explore the distribution of the mutational profiles by risk groups. Only works with binary predictors.
#' @param data The data used for the upstream analysis with an extract column indicating the risk group of each patients.
#' Can be obtained from the output of the riskStrat function.
#' @topHist The name of the most frequently selected genes (note will be overthrown by the geneList argument).
#' @param numGroups The number of groups to be made when stratifying by risk groups. Options are 2,3 and 4 (for now implementing
#' broader version). Default is 2.
#' @param cuts Numeric vector of the points in the distribution of risk scores where groups will be splitted. Needs to be of length
#' numGroups - 1. eg : c(0.25,0.5,0.75) when numgroups is 4. Default is 0.5.
#' @param geneList Optional character vector of gene names to use to generate the pie charts. Default is NULL, which
#' leads to using the 5 most frequently selected features.
#'
#' @return mut_Plot Mutation distribution by features bar plot.
#' @return PieChart Interactive pie chart of the mutation dsitribution using either the most frequently selected features or
#' the manually inputted gene list. Each of the pies represent one of the risk groups.
#' @return GenesUsed Character vector of the features used to make the pie charts.
#'
#' @examples library(OncoCast)
#' test <- OncoCast(data=survData,formula=Surv(time,status)~.,
#'                           method=c("LASSO"),
#'                           runs = 25,cores = 1,sampling = "cv",
#'                           pathResults = "./Test/",studyType = "ExampleRun",save=F)
#' out <- outputSummary(test$LASSO)
#' riskout <- riskStrat(survData,out$average.risk,numGroups = 2,cuts=0.5)
#' data <- riskout$data.out
#' topHits <- out$topHits
#' mut.results <- mutSummary(data,topHits,numGroups = 2,geneList=NULL)
#'
#' @import survival
#' @import ggplot2
#' @import plotly
#' @import plyr
#' @import stats
#' @import reshape2
#' @import scales
#' @import survminer
#' @import data.table


mutSummary <- function(data,average.risk,topHits,numGroups,geneList){

  if(length(grep("time",colnames(data)))  == 1) {LT = FALSE}
  if(length(grep("time",colnames(data)))  == 2) {LT = TRUE}

  ### MAKE CORRESPONDING MUTATION PER GROUP PLOT ###

  mutDistrib <- as.data.frame(matrix(nrow = numGroups,ncol = length(topHits)))
  Groups <- unique(data$RiskGroup)
  rownames(mutDistrib) <- Groups
  colnames(mutDistrib) <- topHits
  for( gene in 1:length(topHits)){
    if(length(data[data$RiskGroup == "Low",match(topHits[gene],colnames(data))]) != 0 ){
      mutDistrib[match("Low",rownames(mutDistrib)),gene] <- sum(data[data$RiskGroup == "Low",match(topHits[gene],colnames(data))])/length(data[data$RiskGroup == "Low",match(topHits[gene],colnames(data))])
    }
    if(length(data[data$RiskGroup == "Low-intermediate",match(topHits[gene],colnames(data))]) != 0 ){
      mutDistrib[match("Low-intermediate",rownames(mutDistrib)),gene] <- sum(data[data$RiskGroup == "Low-intermediate",match(topHits[gene],colnames(data))])/length(data[data$RiskGroup == "Low-intermediate",match(topHits[gene],colnames(data))])
    }
    if(length(data[data$RiskGroup == "Intermediate",match(topHits[gene],colnames(data))]) != 0 ){
      mutDistrib[match("Intermediate",rownames(mutDistrib)),gene] <- sum(data[data$RiskGroup == "Intermediate",match(topHits[gene],colnames(data))])/length(data[data$RiskGroup == "Intermediate",match(topHits[gene],colnames(data))])
    }
    if(length(data[data$RiskGroup == "High-intermediate",match(topHits[gene],colnames(data))]) != 0  ){
      mutDistrib[match("High-intermediate",rownames(mutDistrib)),gene] <- sum(data[data$RiskGroup == "High-intermediate",match(topHits[gene],colnames(data))])/length(data[data$RiskGroup == "High-intermediate",match(topHits[gene],colnames(data))])
    }
    if(length(data[data$RiskGroup == "High",match(topHits[gene],colnames(data))]) != 0 ){
      mutDistrib[match("High",rownames(mutDistrib)),gene] <- sum(data[data$RiskGroup == "High",match(topHits[gene],colnames(data))])/length(data[data$RiskGroup == "High",match(topHits[gene],colnames(data))])
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
    fit.data <- data[,match(c("time1","time2","status",topHits[1:min(5,length(topHits))]),colnames(data))]
    fit.topHits <- coxph(Surv(time1,time2,status)~.,data= fit.data)}
  if(!LT) {
    fit.data <- data[,match(c("time","status",topHits[1:min(5,length(topHits))]),colnames(data))]
    fit.topHits <- coxph(Surv(time,status)~.,data= fit.data)
  }

  ####
  if(length(geneList) == 0) {useGenes <- topHits[1:min(5,length(topHits))]}
  else{
    geneList <- gsub(" ","",geneList)
    useGenes <- geneList}
  pie.data <- data[,match(useGenes,colnames(data))]
  Groups <- as.character(data[,match("RiskGroup",colnames(data))])
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
    colnames(profiles.inter) <- c("Groups","profile","N","Profile")
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


  ###### MUTATION PROFILE ######
  # to <- c(0,10)
  # from <- range(average.risk, na.rm = TRUE, finite = TRUE)
  # RiskScore <- (as.numeric(average.risk)-from[1])/diff(from)*diff(to)+to[1]
  # names(RiskScore) <- names(average.risk)
  # oo=order(RiskScore)
  #
  # if(is.null(geneList)){
  #   geneList <- topHits
  # }
  # mut=data[,match(geneList,colnames(data))]
  # mut=mut[names(RiskScore), ]
  # mut.sorted=t(mut)[,oo]
  #
  # #bw=colorpanel(2, low="white", high="black")
  # n = nrow(data)
  # cols=colorpanel(n, low="green", high="red")
  # profiles <- heatmap.2(mut.sorted, dendrogram="none",  Rowv=NA, Colv=NA,
  #           trace="none", labCol = "",
  #           ColSideColors=cols,
  #           key=FALSE, margins = c(0,7),
  #           col=colorpanel(2, low="white", high="black"),
  #           lwid=c(0.1,4), lhei=c(0.1,4))

  useGenes <- paste0(useGenes,collapse = ",")
  return(list("mut_Plot" = mut_2LVLS,"PieChart"=pie.chart,
              "GenesUsed"=useGenes))#,"MutProfiles"=profiles))
}
