#######################################################################################
################### EFFECT OF taVNS ON HRV - A BAYSIAN META-ANALYSIS ##################
#######################################################################################

################### Helper functions #################################################
#----
#Load files
load(file = "df.RDa")
load(file = "screened.RDa")

#Generate function "shrinkage.ggplot"
  #This function shows the shrinkage of a studies estimate compared to the initial estimate
shrinkage.ggplot <- function(data, i) {
  effect <- seq(-1, 2, length = 200)
  colors <- c("initial estimate" = "blue", "shrinkage estimate" = "red")
  print(ggplot(data = NULL, aes(x=effect,y=bma$dposterior(theta = effect, individual = 1))) + 
    geom_line(aes(x = effect, y = dnorm(effect, mean = bma$y[i], sd = bma$sigma[i]), col = "initial estimate")) + 
    geom_line(aes(x = effect, y = bma$dposterior(theta = effect, individual = 1), col = "shrinkage estimate")) + 
    geom_vline(xintercept = 0, col = "gray") + 
    theme_minimal_hgrid(12) + 
    labs(x = "effect", y = "probability density", color = "legend") + 
    scale_color_manual(values = colors) + 
    ggtitle(bma$labels[i]))
}

#Generate function "posterior.ggPlot"
  #This function generates a posterior and posterior predictive probability density plot
  #xmin and xmax to scale the x-axis
  #bma has to be a bayesmeta object
posterior.ggplot <- function(xmin, xmax, bma) {
  effect <- seq(xmin, xmax, length = 200)
  colors <- c("posterior" = "red", "posterior predictive" = "blue")
  devAskNewPage(ask = FALSE)
  ggplot(data = NULL, aes(x=effect,y=bma$dposterior(mu = effect))) + 
    geom_line(aes(x = effect, y = bma$dposterior(mu = effect), col = "posterior (?)")) + 
    geom_line(aes(x = effect, y = bma$dposterior(mu = effect, predict = T), col = "posterior predictive (???k+1)")) + 
    geom_vline(xintercept = 0, col = "gray") + 
    theme_minimal_hgrid(12) + 
    labs(x = "effect", y = "probability density", color = "legend") + 
    scale_color_manual(values = colors)
}

#Generate function "priorposteriorlikelihood.ggplot"
priorposteriorlikelihood.ggplot <- function(bma) {
  effect <- seq(-2, 2, length = 200)
  colors <- c("posterior" = "#D55E00", "prior" = "#009E73", "likelihood" = "#0072B2")
  devAskNewPage(ask = FALSE)
  ggplot(data = NULL, aes(x=effect,y=bma$dposterior(mu = effect))) + 
    geom_line(aes(x = effect, y = bma$likelihood(mu = effect), col = "likelihood")) +
    geom_line(aes(x = effect, y = bma$dposterior(mu = effect), col = "posterior")) + 
    geom_line(aes(x = effect, y = bma$dprior(mu = effect), col = "prior")) + 
    geom_vline(xintercept = 0, col = "gray") + 
    theme_minimal_hgrid(12) + 
    labs(x = "effect", y = "probability density", color = "legend") + 
    scale_color_manual(values = colors)
}

#Generate function "tauprior.ggplot"
tauprior.ggplot <- function(bma) {
  effect <- seq(0, 2, length = 200)
  devAskNewPage(ask = FALSE)
  ggplot(data = NULL, aes(x=effect,y=bma$dprior(tau = effect))) + 
    geom_line(aes(x = effect, y = bma$dprior(tau = effect))) + 
    geom_vline(xintercept = 0, col = "gray") + 
    theme_minimal_hgrid(12) + 
    labs(x = "tau", y = "probability density") 
}

#Generate function "bayesmeta.master"
  #This function does the full analysis for a whole dataset/subset:
    #aggregation of effect sizes
    #summary of results
    #plots
bayesmeta.master <- function(data) {
  #Aggregate within study effect sizes
  aggES <<- agg(id     = Paper.No,
               es     = yi,
               var    = vi,
               data   = data,
               cor = .5,
               method = "BHHR")
  #Merging aggregated ES with original dataframe 
  MA <<- merge(x = aggES, y = data, by.x = "id", by.y = "Paper.No")
  MA <<- unique(setDT(MA) [sort.list(id)], by = "id")
  MA <<- with(MA, MA[order(MA$es)])
  #Generate bayesmeta-object "bma", which stores all relevant results
  bma <<- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                   tau.prior = function(t) dhalfnormal(t, scale = 0.5), 
                   mu.prior.mean = 0, mu.prior.sd = 4)
  #Show summary of "bma"
  summary(bma)
  #ggPlot: Prior, posterior, and likelihood
  print(priorposteriorlikelihood.ggplot(bma))
  #ggPlot: Original and shrinkage estimates for the all studys seperatly
  for (study in c(1:nrow(MA))){
    shrinkage.ggplot(data = bma, i=study)
  }
  #ggPlot: Posterior and posterior predictive probability density plot
  print(posterior.ggplot(-1,1,bma))
  #Generate basic plots
  forestplot.bayesmeta(bma)
  plot(bma, which=2)
  plot(bma, which=3)
  plot(bma, which=4)
}

# BF robustness check
robustness <- function(MA,SD, tauprior) {
  narrow <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                    tau.prior = tauprior, 
                    mu.prior = c("mean" = 0, "sd" = (SD/2)))
  default <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                    tau.prior = tauprior, 
                    mu.prior = c("mean" = 0, "sd" = 1.5))
  user <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                    tau.prior = tauprior, 
                    mu.prior = c("mean" = 0, "sd" = SD))
  wide <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                    tau.prior = tauprior, 
                    mu.prior = c("mean" = 0, "sd" = SD+1))
  ultrawide <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                         tau.prior = tauprior, 
                         mu.prior = c("mean" = 0, "sd" = SD+2))
  BFS <- c(narrow$bayesfactor[1,2], user$bayesfactor[1,2], wide$bayesfactor[1,2], ultrawide$bayesfactor[1,2])
  SDS <- c(SD/2,SD,SD+1,SD+2)
  names <- c("Narrow", "User","Wide","Ultrawide")
  deflabel <- data.frame(Ref = "Default", val = 1.5, stringsAsFactors = F)
  ggplot(data = NULL, aes(SDS,BFS)) +
    geom_hline(yintercept = 3, color = "grey", size = 0.3, alpha = .6) + 
    geom_hline(yintercept = 10, color = "grey", size = 0.3, alpha = .6) + 
    geom_hline(yintercept = 20, color = "grey", size = 0.3, alpha = .6) + 
    geom_hline(yintercept = 30, color = "grey", size = 0.3, alpha = .6) + 
    geom_vline(xintercept = 1.5, color = "#D55E00", size = 0.3, alpha = .4) +
    geom_point(aes(SDS,BFS), alpha = .75) + 
    geom_point(aes(1.5, default$bayesfactor[1,2]), alpha = .75) +
    scale_x_continuous(breaks = c(SDS, 1.5), limits = c(SDS[1]-.5,SDS[4]+3)) +
    scale_y_continuous(limits = c(0,BFS[4]+10)) +
    geom_line(aes(SDS,BFS), alpha = .4) +
    geom_text(aes(label = names),vjust = -1) +
    geom_text(mapping = aes(x = val, y = 0, label = Ref, hjust = -.2, vjust = -1), data = deflabel, color = "#D55E00") +
    labs(x = "standard deviations", y = "BF01") +
    theme_cowplot(12) +
    annotate("text", label = "Anecdotal evidence for H0", x = SDS[4] + 2, y = 6.5) +
    annotate("text", label = "Substantial evidence for H0", x = SDS[4] + 2, y = 15) +
    annotate("text", label = "Strong evidence for H0", x = SDS[4] + 2, y = 25) +
    annotate("text", label = "Strong evidence for H0", x = SDS[4] + 2, y = (BFS[4]+10+30)/2) +
    annotate("text", label = "Anecdotal evidence for H1", x = SDS[4] + 2, y = 0)
}

#Check for outliers
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}


################### Helper functions for frequentist analysis#########################
#Meta analysis regression with summary, forest and funnel plot
myMareg <- function(formula = es~1,
                    var = var,
                    data,
                    addfit_forest = T) {
  model <- mareg(formula = formula,
                 var = var,
                 data = data,
                 slab = data$study)
  
  result <- list(summary(model),
                 confint(model),
                 metafor::forest.rma(x = model, showweights = T, addfit = addfit_forest,
                                     order = "obs",
                                     xlim=c(-10,8)),
                 funnel(model, xlab = "Observed outcome"))
  listnames <- c("Modelsummary", "Model-Fit", "Forestplot", "Funnelplot")
  names(result) <- listnames
  return(result)
  
}

#Coefficients table
myCoef <- function(coef){
  coef<-round(coef,digits = 3)
  DT::datatable(coef,
                rownames=c("Intercept"),
                colnames=c("b","S.E.","z","lower CI","upper CI","p"))
}

#Model fit
myFit <- function(fit){
  fit<-round(fit,digits = 3)
  DT::datatable(fit)
}

#Model uniqueness
myUni<-function(rand){
  uni<-as.data.frame(rand)
  uni$estimate<-round(uni$estimate,digits = 2)
  uni$ci.lb<-round(uni$ci.lb,digits = 2)
  uni$ci.ub<-round(uni$ci.ub,digits = 2)
  DT::datatable(uni,
                colnames = c("estimate","lower CI","upper CI"))
}


