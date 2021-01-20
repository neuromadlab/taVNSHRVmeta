#######################################################################################
################### EFFECT OF taVNS ON HRV - A BAYESIAN META-ANALYSIS #################
#######################################################################################

################### Shiny App v.2 21.12.2020 UI #######################################

# Load required packages and source helper functions #----
library(shiny)
library(bayesmeta)
library(cowplot)
library(dplyr)
library(DT)
library(data.table)
library(esc)
library(ggplot2)
library(MAd)
library(readr)
library(R.rsp)
library(shinycssloaders)
library(shinythemes)
library(stringr)
library(tidyr)
library(xtable)
source("HelperFunctions.R")
#----
# Define UI
ui <- fluidPage(theme = shinytheme("cosmo"),
                titlePanel(title = div("The influence of taVNS on HRV - a Bayesian living & interactive meta-analysis", 
                                       img(src='ukt.png',style="width: 200px", align = "right")),
                           windowTitle = "taVNS on HRV | Bayesian living & interactive MA"),
                sidebarLayout(
                  sidebarPanel(fluidRow(
                    tabsetPanel(
                      tabPanel("Study criteria",
                               checkboxGroupInput(inputId = "design", label = "Study design", choices = levels(df$Design), selected = levels(df$Design)),
                               checkboxGroupInput(inputId = "control", label = "Control type", choices = levels(df$Control.Type), selected = "sham"),
                               checkboxGroupInput(inputId = "hrvparam", label = "Heart rate variability parameter", choices = levels(df$HRV.Parameter), selected = c("RMSSD", "HF-HRV", "pNN50", "CVT", "RSA")),
                               checkboxGroupInput(inputId = "sample", label = "Sample", choices = levels(df$Sample), selected = "healthy"),
                               sliderInput(inputId = "gender", label = "Percent females", min = min(df$Percent.Females), max = max(df$Percent.Females), value = c(min(df$Percent.Females), max(df$Percent.Females)), step = 1, sep = "", ticks = F),
                               sliderInput(inputId = "age", label = "Age", min = min(df$Mean.Age), max = max(df$Mean.Age), value = c(min(df$Mean.Age), max(df$Mean.Age)), step = 1, sep = "", ticks = F),
                               checkboxGroupInput(inputId = "blind", label = "Blindness of Study", choices = levels(df$Blindness), selected = c("single blind", "double blind")),
                               sliderInput(inputId = "pubyear", label = "Publication year", min = min(df$Publication.Year), max = max(df$Publication.Year), value = c(min(df$Publication.Year), max(df$Publication.Year)), step = 1, sep = "", ticks = F),
                               checkboxGroupInput(inputId = "included", label = "Include/exclude specific studies", choices = levels(df$study), selected = levels(df$study))),
                      tabPanel("taVNS specifications",
                               checkboxGroupInput(inputId = "side", label = "Stimulation side", choices = levels(df$Stimulation.side), selected = levels(df$Stimulation.side)),
                               checkboxGroupInput(inputId = "part", label = "Stimulation site", choices = levels(df$Part.of.the.ear.stimulated), selected = levels(df$Part.of.the.ear.stimulated)),
                               sliderInput(inputId = "stimduration", label = "Stimulation duration (in sec.)", min = min(df$Stimulation.duration.sec.), max = max(df$Stimulation.duration.sec.), value = c(min(df$Stimulation.duration.sec.), max(df$Stimulation.duration.sec.)), step = 10, sep = "", ticks = F),
                               checkboxGroupInput(inputId = "intensity1", label = "Intensity chosen", choices = levels(df$Intensity.fixed.or.mean), selected = levels(df$Intensity.fixed.or.mean)),
                               checkboxGroupInput(inputId = "device", label = "taVNS device", choices = levels(df$taVNS.Device), selected = levels(df$taVNS.Device)),
                               sliderInput(inputId = "intensity2", label = "Intensity (in mA)", min = min(df$Mean.Intensity..mA.), max = max(df$Mean.Intensity..mA.), value = c(min(df$Mean.Intensity..mA.), max(df$Mean.Intensity..mA.)), step = 0.1, sep = "", ticks = F)),
                      tabPanel("Prior specifications",
                               numericInput(inputId = "mupriormean", label = "µ prior mean", value = 0, step = 0.1),
                               numericInput(inputId = "mupriorsd", label = "µ prior standard deviation", value = 1.5, step = 0.1, min = 0),
                               radioButtons(inputId = "robust", label = "µ Bayes Factor robustness check", choices = c(No = "No", Yes = "Yes")),
                               p("Selecting 'Yes' will lead to an increase in computation time. Bayes Factors over a variety of prior standard deviations will be plotted."), hr(),
                               radioButtons(inputId = "tauprior", label = "τ prior", choices = c("Half cauchy", "Half student t","uniform", "sqrt", "Jeffreys", "BergerDeely", "conventional", "DuMouchel", "shrinkage", "I2")),
                               numericInput(inputId="scaletau", label="τ prior scale (for half cauchy or half student t)", value=0.5, step=0.05),
                               a("Further information on choosing an appropriate τ prior.", href="https://cran.r-project.org/web/packages/bayesmeta/bayesmeta.pdf", target = "_blank")),
                      hr(),
                      submitButton("Re-Calculate Meta-Analysis", icon("sync"))
                    ))),
                  mainPanel(
                    fluidRow(
                      tabsetPanel(
                        tabPanel("Explanation", br(),
                                 h3("Welcome to the living & interactive meta-analysis on the influence of transcutaneous auricular vagus nerve stimulation (taVNS) on heart rate variablilty (HRV)."), br(),
                                 h4("Purpose:"),
                                 p("Non-invasive brain stimulation techniques, such as transcutaneous auricular vagus nerve stimulation (taVNS), have considerable potential for clinical use. Beneficial effects of taVNS have been demonstrated on symptoms in patients with mental or neurological disorders as well as transdiagnostic dimensions, including mood and motivation. However, since taVNS research is still an emerging field, the underlying neurophysiological processes are not yet fully understood, and the replicability of findings on biomarkers of taVNS effects has been questioned. Here, we perform a living Bayesian random effects meta-analysis to synthesize the current evidence concerning the effects of taVNS on heart rate variability (HRV), a candidate biomarker that has, so far, received most attention in the field. To keep the synthesis of evidence transparent and up to date as new studies are being published, we developed a Shiny web app that regularly incorporates new results and enables users to modify study selection criteria to evaluate the robustness of the inference across potential confounds. By increasing transparency and timeliness, we believe that the concept of living meta-analyses can lead to transformational benefits in emerging fields such as non-invasive brain stimulation."), br(),
                                 h4("Explanation:"),
                                 p("Select inclusion criteria and prior settings in the panel on the left. Click 'Re-Calculate Meta-Analysis' to rerun the analysis. Results are displayed in the central panels."), br(),
                                 h4("Paper:"),
                                 ("This app accompanies the following "),
                                 a("paper", href="https://doi.org/10.1101/2021.01.18.426704", target = "_blank"), 
                                 ("where a more detailed explanation and an exemplary analysis can be found. Inclusion criteria and prior settings of the exemplary analysis correspond to the default criteria in the app."), br(), br(),
                                 h4("Contact:"),
                                 p("The app is maintained by neuroMADLAB (Neuroscience of Motivation, Action, & Desire) at the University of Tuebingen, Germany."),
                                 p("Contact/Visit us:"),
                                   a(shiny::actionButton(inputId = "email1", 
                                                         label = "Mail", 
                                                         icon = icon("envelope", lib = "font-awesome")),
                                     href="mailto:neuromadlab@klinikum.uni-tuebingen.de"),
                                   a(shiny::actionButton(inputId = "website", 
                                                         label = "Web", 
                                                         icon = icon("globe", lib = "font-awesome")),
                                     href="http://www.neuromadlab.com"),
                                   a(shiny::actionButton(inputId = "website", 
                                                         label = "Twitter", 
                                                         icon = icon("twitter", lib = "font-awesome")),
                                     href="http://www.twitter.com/neuromadlab")
                                 ),
                        tabPanel("Study overview", br(),
                                 h4("This table provides an overview of the studies included in the analysis:"),
                                 textOutput("warning"), br(),
                                 DT::dataTableOutput("studies") %>% withSpinner(type = 6, color = "#3498DB"), br(),
                                 h4("Abbreviations:"), p("RMSSD: root mean square of successive differences; 
                                   HF-HRV: high frequency heart rate variability; 
                                   pNN50: percentage of successive normal sinus RR intervals more than 50ms;
                                   LF/HF Ratio: low frequency to high frequency ratio;
                                   SDNN: standard deviation of all R-to-R intervals;
                                   CVT: cardiac vagal tone calculated with 'phase shift demodulation';
                                   RSA: respiratory sinus arrhythmia.")),
                        tabPanel("Outlier check", br(),
                                 h4("Boxplot graph:"), plotOutput("boxplot") %>% withSpinner(type = 6, color = "#3498DB")),
                        tabPanel("Forest plot", br(),
                                 h4("Forest plot with 95% credible intervals:"),
                                 plotOutput("forest") %>% withSpinner(type = 6, color = "#3498DB")),
                        tabPanel("Funnel plot", br(),
                                 h4("Funnel plot to assess publication bias:"),
                                 plotOutput("funnel")  %>% withSpinner(type = 6, color = "#3498DB")),
                        tabPanel("Statistics", br(),
                                 h4("Bayes factors:"), verbatimTextOutput("bf") %>% withSpinner(type = 6, color = "#3498DB"),
                                 p("Bayes factors are only computed if the priors for τ and μ are proper."), br(),
                                 h4("Marginal posterior summary:"), verbatimTextOutput("summary") %>% withSpinner(type = 6, color = "#3498DB"), br(),
                                 h4("Maximum-likelihood:"), verbatimTextOutput("ML") %>% withSpinner(type = 6, color = "#3498DB"), br(),
                                 h4("Maximum-a-posteriori:"), verbatimTextOutput("MAP") %>% withSpinner(type = 6, color = "#3498DB")),
                        tabPanel("Full texts screened", br(),
                                 h4("This table provides an overview of all full-text we assessed for eligibility:"),
                                 DT::dataTableOutput("screened") %>% withSpinner(type = 6, color = "#3498DB"), br(),
                                 h4("Important note:"),
                                 ("If we missed a paper"),
                                 ("or you published a paper which is not yet shown on this list or"),
                                 ("if you have unpublished HRV data of your taVNS experiments, please"),
                                 a("contact us.", href="mailto:neuromadlab@klinikum.uni-tuebingen.de", target = "_blank")),
                        tabPanel("Additional plots", br(),
                                 h4("Joint posterior density of heterogeneity τ and effect μ:"), plotOutput("joint") %>%  withSpinner(type = 6, color = "#3498DB"),
                                 p("Darker shading corresponds to higher probability density."),
                                 p("Red lines indicate (approximate) 2-dimensional credible regions,"),
                                 p("green lines show marginal posterior medians and 95% credible intervals,"),
                                 p("blue lines show conditional posterior mean effect as a function of the heterogeneity along with a 95% interval."),
                                 p("Red cross (+): posterior mode"),
                                 p("Pink cross (x): ML estimate"), br(),
                                 h4("Prior, posterior, & likelihood:"), plotOutput("evupdate") %>% withSpinner(type = 6, color = "#3498DB"), br(),
                                 h4("τ prior distribution:"), plotOutput("taupriorplot") %>% withSpinner(type = 6, color = "#3498DB")),
                        tabPanel("Bayes factor robustness check", br(),
                                 h4("Bayes Factors over a variety of prior standard deviations:"),
                                 p("Will only be computed if 'Yes' is selected for 'µ Bayes Factor robustness check' and the priors for τ and μ are proper."),
                                 conditionalPanel(condition = "input.robust == 'Yes'", 
                                 plotOutput("robustplot") %>% withSpinner(type = 6, color = "#3498DB"),
                                 p("Default: 1.5 (orange horizontal line)"),
                                 p("Narrow: user selected standard deviation / 2"),
                                 p("User: user selected standard deviation"),
                                 p("Wide: user selected standard deviation + 1"),
                                 p("Ultrawide: user selected standard deviation + 2"),
                                 p("Interpretations of Bayes factors are based on Jeffreys (1961) und should be considered with caution.")))
                          )
                        )
                      )
                    )
                  )
                
