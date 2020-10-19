#######################################################################################
################### EFFECT OF taVNS ON HRV - A BAYSIAN META-ANALYSIS ##################
#######################################################################################

################### Shiny App v.2 20.08.2020 SERVER ###################################
# Initialize #----
# rm(list=ls())
# if (!require("pacman")) install.packages("pacman")
# library(pacman)
# pacman::p_load(bayesmeta, car, cowplot, dplyr, DT, esc, data.table, 
#                ggplot2, knitr, MAd, readr, readxl, R.rsp, shiny, 
#                shinycssloaders, shinythemes, stringr, tidyr, xtable) 
# load(file = "df.RDa")
# load(file = "screened.RDa")
# source("HelperFunctions.R")
#----
# Define server logic
server <- function(input, output) {
  # Create parameter reactive for HRV parameter columns in study overview
  parameters <- reactive({input$hrvparam})
  # Create MA reactive for all outputs
  MA <- reactive({
    # Create subset based on chosen inclusion criteria
    df_sub <- df %>% filter(df$Design %in% input$design
                            & df$HRV.Parameter %in% input$hrvparam
                            & df$Sample %in% input$sample
                            & df$Age.Category %in% input$age
                            & df$Blindness %in% input$blind
                            & df$Stimulation.side %in% input$side
                            & df$Publication.Year >= input$pubyear
                            & df$Control.Type %in% input$control
                            & df$Part.of.the.ear.stimulated %in% input$part
                            & df$Intensity.fixed.or.mean %in% input$intensity1
                            & df$Mean.Intensity..mA. >= input$intensity2
                            & df$Stimulation.duration.sec. >= input$stimduration
                            & df$taVNS.Device %in% input$device
                            & df$study %in% input$included)
    # Aggregate effect sizes
    aggES <- agg(id     = Paper.No,
                 es     = yi,
                 var    = vi,
                 data   = df_sub,
                 cor = .5,
                 method = "BHHR")
    # Merging aggregated ES with original dataframe 
    MA <- merge(x = aggES, y = df_sub, by.x = "id", by.y = "Paper.No")
    MA <- unique(setDT(MA) [sort.list(id)], by = "id")
    MA <- with(MA, MA[order(MA$es)])
  })
  # Create bma reactive for all outputs
  bma <- reactive({
    # Generate bayesmeta-object "bma", which stores all relevant results, depending on tau prior chosen
    if (input$tauprior == "Half cauchy") {
      bma <- bayesmeta(y = MA()$es,sigma = sqrt(MA()$var), labels = MA()$study, 
                       tau.prior = function(t) dhalfcauchy(t, scale = input$scaletau), 
                       mu.prior = c("mean" = input$mupriormean, "sd" = input$mupriorsd))
    } else if (input$tauprior == "Half student t") {
      bma <- bayesmeta(y = MA()$es,sigma = sqrt(MA()$var), labels = MA()$study, 
                       tau.prior = function(t) dhalfnormal(t, scale = input$scaletau), 
                       mu.prior = c("mean" = input$mupriormean, "sd" = input$mupriorsd))
    } else {
      bma <- bayesmeta(y = MA()$es,sigma = sqrt(MA()$var), labels = MA()$study, 
                       tau.prior = input$tauprior, 
                       mu.prior = c("mean" = input$mupriormean, "sd" = input$mupriorsd))
    }
  })
  # Study overview
  output$studies <- DT::renderDataTable({
    MA <- as.data.frame(MA())
    MA <- MA %>% mutate_if(is.numeric, round, digits=3)
    parameters <- base::subset(MA, select = c(parameters()))
    criteria <- base::subset(MA, select = c(study, es, var, Design, Sample, Age.Category, Blindness, Total.N, Stimulation.side, Part.of.the.ear.stimulated, taVNS.Device, Comment))
    colnames(criteria) <- c("Study", "Hedges' g", "Variance", "Design", "Sample", "Age group", "Blindness", "Total N", "Stimulation side", "Stimulation site", "taVNS device", "Comment")
    MAclean <- cbind(criteria,parameters)
    DT::datatable(MAclean,
                  options = list(dom = 't', 
                                 pageLength = nrow(MAclean)))
  })
  # Warning message if 3 or less studies are included
  output$warning <- renderPrint({
    MA <- as.data.frame(MA())
    if (nrow(MA) < 4) {print('WARNING: With the chosen inclusion criteria, 3 or less studies will be included in the analysis.')}
  })
  # Outliers
  output$boxplot <- renderPlot({
    MAo <- MA() %>% tibble::rownames_to_column(var = "outlier") %>% mutate(is_outlier=ifelse(is_outlier(es), es, as.numeric(NA)))
    MAo$study[which(is.na(MAo$is_outlier))] <- as.numeric(NA)
    ggplot(MAo, aes(x = factor(0), es)) +
      geom_boxplot(outlier.size = 3.5, outlier.colour = "#D55E00", outlier.shape = 18, fill = "lightgrey") +
      geom_text(aes(label=study),na.rm = T, nudge_y = 0.02, nudge_x = 0.05) +
      stat_boxplot(geom="errorbar", width = 0.05) +
      scale_x_discrete(breaks = NULL) +
      xlab(NULL) + ylab("Hedges' g") +
      theme_minimal_hgrid(12)
  }, width = 600, height = 600)
  
  # Forest Plot
  output$forest <- renderPlot({
    forestplot.bayesmeta(bma(), xlab = "Hedges' g")
  })
  # Funnel Plot
  output$funnel <- renderPlot({
    funnel.bayesmeta(bma(), main = "")
  })
  # Statistics
  output$bf <- renderPrint ({
    bma()$bayesfactor[1,]
  })
  output$summary <- renderPrint({
    bma()$summary
  })
  output$ML <- renderPrint({
    bma()$ML
  })
  output$MAP <- renderPrint({
    bma()$MAP
  })
  # Full texts screened
  output$screened <- DT::renderDataTable({
    datatable(screened, escape = F, options = list(columnDefs = list(list(
                                         targets = 2,
                                         render = JS(
                                           "function(data, type, row, meta) {",
                                           "return type === 'display' && data != null && data.length > 30 ?",
                                           "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
                                           "}")
                                       ))),
              class = "display")
  })
  # Additional plots
  output$evupdate <- renderPlot({
    priorposteriorlikelihood.ggplot(bma())
  }, width = 800)
  output$joint <- renderPlot({
    plot.bayesmeta(bma(), which=2, main = "")
  }, width = 800)
  output$taupriorplot <- renderPlot({
    tauprior.ggplot(bma())
  }, width = 800)
  output$robustplot <- renderPlot({
    if (input$robust == "Yes" &
        input$tauprior == "Half cauchy") {
      robustness(MA(),SD = input$mupriorsd, tauprior = function(t) dhalfcauchy(t, scale = input$scaletau))
    } else if (input$robust == "Yes" &
               input$tauprior == "Half student t") {
      robustness(MA(),SD = input$mupriorsd, tauprior = function(t) dhalfnormal(t, scale = input$scaletau))
    } 
  }, width = 800)
}
