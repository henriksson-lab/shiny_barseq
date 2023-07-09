library(plotly)
library(shiny)
library(ggplot2)

################################################################################
########### Samplemeta #########################################################
################################################################################

tab_samplemeta <- sidebarLayout(
  
  sidebarPanel(
    
    selectInput(
      inputId = "samplemeta_pool",
      label = "Pool:",
      selectize = FALSE,
      multiple = FALSE,
      choices = names(all_samplemeta), 
      selected = names(all_samplemeta)[1]
    )

  ),
  mainPanel(
    h3("UMAPs pool overview"),
    plotOutput(outputId = "plotSamplemetaUmap", height = "700px")
  )
)


################################################################################
########### Gene viewer ########################################################
################################################################################

tab_grstats <- sidebarLayout(
  
  sidebarPanel(
    
    selectInput(
      inputId = "grstats_pool",
      label = "Pool:",
      selectize = FALSE,
      multiple = FALSE,
      choices = names(all_grstats), 
      selected = names(all_grstats)[1]
    ),

    selectInput(
      inputId = "grstats_volcano",
      label = "Volcano:",
      selectize = TRUE,
      multiple = FALSE,
      choices = c(""), 
      selected = NULL
    ),
    
    
    selectInput(
      inputId = "grstats_scatter",
      label = "Scatter:",
      selectize = TRUE,
      multiple = FALSE,
      choices = c(""), 
      selected = NULL
    ),
    
    
    selectInput(
      inputId = "grstats_gene",
      label = "Gene:",
      selectize = TRUE,
      multiple = FALSE,
      choices = c(""), 
      selected = NULL
    ),
    
    selectInput(
      inputId = "grstats_genesummary",
      label = "Show gene summaries:",
      selectize = FALSE,
      multiple = FALSE,
      choices = c("Per condition","Per mouse"), 
      selected = "Per condition"
    )

  ),
  mainPanel(
    
    h3("Volcano plot for one pool"),
    plotlyOutput("plot_grstats_volcano", height = "400px"),

    h3("Pool comparisons"),
    plotlyOutput("plot_grstats_scatterplot", height = "400px"),
    
    h3("Relative abundance over time"),
    plotlyOutput("plot_grstats_tcplot", height = "400px"),
    
  )
)

################################################################################
########### About page #########################################################
################################################################################

tab_about <- verbatimTextOutput("todo description of project here")




################################################################################
########### Total page #########################################################
################################################################################


ui <- fluidPage(
  
  titlePanel("Bushell lab barseq viewer"),

  tabsetPanel(type = "tabs",
              tabPanel("GRstats", tab_grstats),
              tabPanel("Samplemeta", tab_samplemeta),
              tabPanel("About", tab_about),
  )
  
)



