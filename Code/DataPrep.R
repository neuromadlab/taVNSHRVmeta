#######################################################################################
################### EFFECT OF taVNS ON HRV - A BAYSIAN META-ANALYSIS ##################
#######################################################################################

################### Data preparation ##################################################
#---- 
#Clear global environment and load required packages
rm(list=ls())
dev.off()
if (!require("pacman")) install.packages("pacman")
pacman::p_load(bayesmeta, compute.es, cowplot, data.table, dplyr, esc, forestplot, 
               ggplot2, knitr, MAd, metafor, reactable, readr, readxl, rmarkdown, 
               R.rsp, stringr, tidyr) 

#Import masterfile, in which articles were coded
Master <- read_excel("Master.xlsx")

#Generate column "study" serving as publication identifier (last name of 1st author + publication year)
firstauthor <- word(string = Master$Author, sep = ",") #Grab last name of the 1st author
study<-as.data.frame(cbind(firstauthor,Master$Publication.Year)) #Dataframe with 1st authors last name and publication year
study <- study %>% unite(study, sep = ", ") #Unite it into 1 column
Master <- cbind(study,Master)

#Table showing reasons for exclusions
excluded <- Master %>% filter(Master$Included == 0)
excluded <- unique(setDT(excluded) [sort.list(study)], by = "DOI")
reasons_exc <- table(excluded$Reason.for.exclusion)

#Generate data frame "df" consisting only of included studys 
df <- Master %>% filter(Master$Included==1)

#Delete irrelevant columns
df <- base::subset(df, select = -c(Reason.for.exclusion, Included, Pages, Publication.Title, 
                                   Issue, Volume, HRV.controlled.for.respiration, 
                                   HRV.measurement.duration..min.))

#Calculate effect sizes (Hedges' g) and store them in "es"
es <- as.data.frame(es<-escalc(m1i = df$Mean_HRV_intervention,
            sd1i = df$SD_HRV_intervention,
            n1i = df$N_Intervention,
            m2i = df$Mean_HRV_control,
            sd2i = df$SD_HRV_control,
            n2i = df$N_Control,
            measure = "SMD",
            yi, sei, vi))
   #Add missing ES manually (Burger, 2019a; Burger, 2019b)
      #Burger2019b
         esc_t(0.81,0.42,58,29,29,"g")
         es[16,1:2] <- c(0.2099,0.2634)
      #Burger2019a
         esc_t(p = 0.99,totaln = 97, grp1n = 48,grp2n = 49, es.type = "g")
         es[13,1:2] <- c(0.0025,0.2031)
         esc_t(p = 0.44,totaln = 97, grp1n = 48,grp2n = 49, es.type = "g")
         es[14,1:2] <- c(0.1562,0.2034)
         esc_t(p = 0.83,totaln = 97, grp1n = 48,grp2n = 49, es.type = "g")
         es[15,1:2] <- c(0.0434,0.2031)
   #Bind "es" to "df"
   df<-cbind(es,df)

#Generate factor variables
df$Design<-factor(df$Design,levels = c(0,1),labels = c("between","within"))
df$Sample<-factor(df$Sample,levels = c(0,4,5,7,8,9,10),labels = c("healthy","high worriers","Depression patients","Hypertension patients", "Pancreatitis patients", "Diastolic dysfunction patients", "Parkinson's disease patients"))
df$Blindness<-factor(df$Blindness,levels = c(0:2),labels = c("double blind","single blind","transparent"))
df$Age.Category<-factor(df$Age.Category,levels = c(0,2,3),labels = c("adults","elderly", "adolescents"))
df$Country<-factor(df$Country,levels = c(0,2:6,8),labels = c("US","Brazil","UK","Netherlands","Belgium","Germany","Denmark"))
df$Stimulation.side<-factor(df$Stimulation.side,levels = c(0,1,3),labels = c("right","left", "bilateral"))
df$Part.of.the.ear.stimulated<-factor(df$Part.of.the.ear.stimulated,levels = c(0:1),labels = c("cymba concha","tragus"))
df$Control.Type<-factor(df$Control.Type,levels = c(0:1),labels = c("sham","pre stimulation baseline"))
df$taVNS.Device<-factor(df$taVNS.Device,levels = c(0:2),labels = c("Cerbomed","TENS","other"))
df$Intensity.fixed.or.mean<-factor(df$Intensity.fixed.or.mean,levels = c(0:1),labels = c("individually","fixed"))
df$Respiratory.Gated<-factor(df$Respiratory.Gated,levels = c(0,1),labels = c("no","yes"))
df$HRV.Recording.Method<-factor(df$HRV.Recording.Method,levels = c(0,1,3),labels = c("ECG (glued)", "ECG (chest strap)", "Photoplethysmography"))
df$HRV.Parameter<-factor(df$HRV.Parameter,levels = c(0:5),labels = c("RMSSD","HF-HRV","pNN50","LF/HF Ratio","SDNN","CVT"))
df$HRV.baseline<-factor(df$HRV.baseline,levels = c(0:1), labels = c("task baseline", "resting baseline"))
df$Direction.of.effect<-factor(df$Direction.of.effect,levels = c(0:2),labels = c("null","positive","negative"))
df$Inverse<-factor(df$Inverse,levels = c(0,1),labels = c("no","yes"))
df$study<-as.factor(df$study)
df$RMSSD<-factor(df$RMSSD, levels = c(0,1), labels = c(" ", "x"))
df$`HF-HRV`<-factor(df$`HF-HRV`, levels = c(0,1), labels = c(" ", "x"))
df$pNN50<-factor(df$pNN50, levels = c(0,1), labels = c(" ", "x"))
df$`LF/HF Ratio`<-factor(df$`LF/HF Ratio`, levels = c(0,1), labels = c(" ", "x"))
df$SDNN<-factor(df$SDNN, levels = c(0,1), labels = c(" ", "x"))
df$CVT<-factor(df$CVT, levels = c(0,1), labels = c(" ", "x"))

#Change publication year for Burger 2019a and 2019b, respectively
df[13:16,4] <- 2019
df$Publication.Year <- as.numeric(df$Publication.Year)

#Inverse LF/HF Ratio ES directions (because lower LF/HF-ratio means higher HRV)
df$yi <- ifelse(df$Inverse=="yes", df$yi*(-1), df$yi*1)

#Add "https://doi.org/" infront of dois
Master$DOI <- paste0("https://doi.org/", Master$DOI)
Master$DOI <- paste0("<a href='", Master$DOI,"'target='_blank'>",Master$DOI,"</a>")


#Write "df" as .csv (currently deactivated)
   write_delim(df, "df.csv", delim = ";")

#Save "df" as .RDa
save(df, file = "df.RDa")

#All screened studies
screened <- Master[,c(1,7,9,14,15)]
screened$Included <- factor(screened$Included,levels = c(0,1), labels = c("no","yes"))
screened <- unique(setDT(screened) [sort.list(study)], by = "DOI")
colnames(screened) <- c("Study", "Title", "DOI", "Included", "Reason for exclusion")
save(screened, file = "screened.RDa")

