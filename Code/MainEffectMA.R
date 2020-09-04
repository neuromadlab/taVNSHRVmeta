#######################################################################################
################### EFFECT OF taVNS ON HRV - A BAYSIAN META-ANALYSIS ##################
#######################################################################################

################### Main Effect Meta Analysis #########################################
#----
#Clear global environment and load required packages
rm(list=ls())
dev.off()
if (!require("pacman")) install.packages("pacman")
pacman::p_load(bayesmeta, car, compute.es, cowplot, data.table, dplyr, esc, forestplot, 
               ggplot2, ggstatsplot, knitr, MAd, metafor, outliers, reactable, readr, readxl, 
               R.rsp, sjPlot, stringr, tidyr) 
                
#Source helper functions
source("HelperFunctions.R")

#Load "df.RDa"
load(file = "df.RDa")

#Filter sham and healthy samples
df <- df %>% filter(df$Sample %in% "healthy"
                    & df$Control.Type %in% "sham")

#Aggregate within study effect sizes
aggES <- agg(id     = df$Paper.No,
             es     = yi,
             var    = vi,
             data   = df,
             cor = .5,
             method = "BHHR")

#Merging aggregated ES with original dataframe 
MA <- merge(x = aggES, y = df, by.x = "id", by.y = "Paper.No")
MA <- unique(setDT(MA) [sort.list(id)], by = "id")
MA <- with(MA, MA[order(MA$es)])
  #write_delim(MA, "~/Documents/Uni Salzburg/Praktikum Master/MA_taVNS_HR/MA_R_Project/DF4jasp.csv", delim = ";")

#Generate bayesmeta-object "bma", which stores all relevant results
bma <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                tau.prior = function(t) dhalfcauchy(t, scale = 0.5), 
                mu.prior = c("mean" = 0, "sd" = 1.5))

#Robustness check
 # robustness(MA, 1.5, function(t) dhalfcauchy(t, scale = 0.5))

#Summary tables
summary(bma)

#Plot: Original and shrinkage estimates for the all studys separatly
  #for (study in c(1:nrow(MA))){
   # shrinkage.ggplot(data = bma, i=study)
    #}

#Plot: Prior, posterior and likelihood
# priorposteriorlikelihood.ggplot(bma)

#ggPlot: Posterior and posterior predictive probability density plot
# posterior.ggplot(-1,1,bma)

#Generate basic plots
forestplot.bayesmeta(bma)
plot(bma, which=2)
plot(bma, which=3)
plot(bma, which=4)

#Funnel plot
funnel.bayesmeta(bma, main = "")

#Study weights
bma$weights

#Posterior predicitve p values
#pppvalue(bma, n = 20, alternative = "two.sided")

#Overview table 1
overview<-data.frame(MA[,c(6, 67, 14, 17,20,24,25)])

#Boxplot outliers
MAo <- MA %>% tibble::rownames_to_column(var = "outlier") %>% mutate(is_outlier=ifelse(is_outlier(es), es, as.numeric(NA)))
MAo$study[which(is.na(MAo$is_outlier))] <- as.numeric(NA)
ggplot(MAo, aes(x = factor(0), es)) +
  geom_boxplot(outlier.size = 3, position = "dodge", fill = "lightgrey") +
  geom_text(aes(label=study),na.rm = T, nudge_y = 0.02, nudge_x = 0.05) +
  stat_boxplot(geom = "errorbar", width = 0.05) +
  scale_x_discrete(breaks = NULL) +
  xlab(NULL) + ylab("Hedges' g") +
  theme_minimal_hgrid(12)


a <- car::Boxplot(MA$es, ylab = "Hedges' g", id = T)
outliernames <- cbind(rownames(MA)[a],MA[a,6])
outliernames <- cbind(outliernames, MA[a,2])
outliernames <- outliernames %>% mutate_if(is.numeric, round, digits=3)
DT::datatable(outliernames, colnames = c("ID","Study", "Hedges' g"), rownames = F,
          options = list(dom = 't', 
                         pageLength = nrow(outliernames)))
MAoutliersremoved <- MA[-1,] # Gancheva2018
boxplot(MAoutliersremoved$es)

#Normality 
gghistostats(MA, es, xlab = "Hedges' g",
             ggtheme = cowplot::theme_minimal_hgrid(12))
# ggbetweenstats(MA, Design, es, plot.type = "box",
#                outlier.tagging = T,
#                outlier.label = study,
#                ylab = "Hedges' g")

#Results for outliers removed
bma2 <- bayesmeta(y = MAoutliersremoved$es,sigma = sqrt(MAoutliersremoved$var), labels = MAoutliersremoved$study, 
                 tau.prior = function(t) dhalfcauchy(t, scale = 0.5), 
                 mu.prior = c("mean" = 0, "sd" = 1.5))
summary(bma2)
forestplot.bayesmeta(bma2, main="")
plot(bma2, which=2)
plot(bma2, which=3)
plot(bma2, which=4)
funnel.bayesmeta(bma2,main="")

