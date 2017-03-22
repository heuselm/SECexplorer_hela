library(shiny)
library(ggplot2)
library(plotly)
library(data.table)


# data preparation
load("data.rda")

annotations <- names(trace_annotation_cum)
for (i in seq_along(annotations)){
  assign(annotations[i], trace_annotation_cum[[i]])
}

shinyUI(pageWithSidebar(

  
  # Application title
  headerPanel("SEC-SWATH-Explorer Hela CCL2 cell-cycle arrest experiment Interphase vs. Mitosis"),
  
  # Sidebar with a slider input for number of observations
  sidebarPanel(
    width = 2,
    selectInput("fcolumn",
                "Choose Annotation Column for filtering",
                annotations,
                selected = "Entry_name"),
    
    uiOutput("fcolumnvalues"), #The possible choices of this field are calculated on the server side and passed over by uiOutput
    
    checkboxInput("remove_legend", label = "Remove plot legends (aligns x axes)",
              value = FALSE),
    
    checkboxInput("logscale", "LOG10 Y-Axis", value = FALSE),
    checkboxInput("show_monomers", "Indicate monomer expected fraction", value = TRUE),
    selectizeInput("pepplot_replicate", label = "Select replicate for peptide level plot",
                   choices = c(1:3), selected = 1, multiple = FALSE),
    
    p("Select annotation column and value to display matching SEC-SWATH chromatograms"),
    p("Chromatograms obtained from mitotic cells are red, those from cells in interphase black"),
    p("The lower 2 plots allow visualization of all the peptide chromatograms obtained.
      Choose the replicate for which you want to see all peptide level info.
      The more peptides support the shift, the higher the confidence."),
    p("Uniprot information for the current protein is summarized in the lower table."),
    p("To align the plots, select 'none' for the legend position in the box above."),
    p("DISCLAIMER: THIS DATA IS UNPUBLISHED WORK - share only with permission of the creators."),
    p("Amon S., Koehler, M., Heusel M., Kutay U., Aebersold R."),
    p("Contact: heuselm@imsb.biol.ethz.ch")
    
    ),
  
  # Show a plot of the generated distribution
  mainPanel(
    tabsetPanel(tabPanel("Chromatogram View",
                plotlyOutput("protplot"),
                plotlyOutput("pepplot"),
                dataTableOutput("table")),
                tabPanel("Shift Score Map",
                         plotlyOutput("shiftmap_combined"),
                         dataTableOutput("shift_scores")),
                tabPanel("Pearson Shift Score Map", plotlyOutput("shiftmap_pearson")),
                tabPanel("FracDiff Shift Score Map", plotlyOutput("shiftmap_fracdiff"))
                )
    
  )
))