#######################################################################################
################### EFFECT OF taVNS ON HRV - A BAYSIAN META-ANALYSIS ##################
#######################################################################################

################### Frequentist Meta Analysis #########################################
#----
#Clear global environment and load required packages
rm(list=ls())
dev.off()
if (!require("pacman")) install.packages("pacman")
pacman::p_load(bayesmeta, compute.es, cowplot, data.table, dplyr, esc, forestplot, 
               ggplot2, knitr, MAd, metafor, reactable, readr, readxl, rmarkdown, 
               R.rsp, stringr, tidyr) 

#Source helper functions
source("HelperFunctions.R")

#Load "df.RDa"
load(file = "df.RDa")

#Filter sham and healthy samples
df <- df %>% filter(df$Sample %in% "healthy"
                  & df$Control.Type %in% "sham")
                    
#Aggregate within study effect sizes
aggES <- agg(id     = Paper.No,
             es     = yi,
             var    = vi,
             data   = df,
             cor = .5,
             method = "BHHR")

#Merging aggregated ES with original dataframe 
MA <- merge(x = aggES, y = df, by.x = "id", by.y = "Paper.No")
MA <- unique(setDT(MA) [sort.list(id)], by = "id")
MA <- with(MA, MA[order(MA$es)])

#Meta regression
m0<-myMareg(data = MA)[1:2]

#Model summary
myCoef(m0$Modelsummary$coef)

#Model fit
myFit(m0$Modelsummary$fit)

#Model uniqueness
myUni(m0$`Model-Fit`$random)

