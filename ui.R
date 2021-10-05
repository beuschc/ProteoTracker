library(shiny)
library(shinydashboard)
library(tidyverse)
library(ggalluvial)
library(ggrepel)
library(shinycssloaders)
library(DT)
library(Cairo)
library(broom)
library(gridExtra)
library(cowplot)
library(metap)
library(svglite)

options(shiny.usecairo = T)
options(dplyr.summarise.inform = F)

cell.lines <- NULL

ui <- dashboardPage(
    title = 'ProteoTracker',
  dashboardHeader(
    title = img(src = 'ProteoTracker.jpg', height = 100, align = 'left'),
    titleWidth = 350
  ),
  
  dashboardSidebar(
    width = 350,
    
    #h3('Upload own data'),
    radioButtons('checkGroup', 
                       h4('Please specify the data'),
                       choices  = list('Use ProteTracker data' = 1,
                                       'Upload own data' = 2),
                       selected = character(0)),
    
    conditionalPanel(condition = 'input.checkGroup == 2',
                     fileInput('upload.file', 'Please choose tsv file',
                               accept = c('text/tsv', '.tsv'))),
    
    selectInput('in.t1.denominator', 'Please specify the denominator for T1', cell.lines, multiple = F),
    selectInput('in.t1.numerator', 'Please specify the numerator for T1', cell.lines, multiple = F),
    
    selectInput('in.t2.denominator', 'Please specify the denominator for T2', cell.lines, multiple = F),
    selectInput('in.t2.numerator', 'Please specify the numerator for T2', cell.lines, multiple = F),
    
    selectInput('in.t3.denominator', 'Please specify the denominator for T3', cell.lines, multiple = F),
    selectInput('in.t3.numerator', 'Please specify the numerator for T3', cell.lines, multiple = F),
    
    sliderInput("fisher.cutoff", "Specify cutoff for fisher significance", min = 0.001, max = 0.5, value = 0.05, step = 0.001),
    
    uiOutput('filterGenes'),
    
    uiOutput('action'),

    p(strong(uiOutput('citation')))
  ),
  dashboardBody(
    tags$head(tags$style(HTML('
                                .skin-blue .main-header .navbar {
                                background-color: #870052;
                                }
                                /* logo */
                                .skin-blue .main-header .logo {
                                background-color: #ffffff ;
                                }
                                /* logo when hovered */
                                .skin-blue .main-header .logo:hover {
                                background-color: #ffffff ;
                                }
                                /* main sidebar */
                                .skin-blue .main-sidebar {
                                background-color: #808080;
                                }                         
                                /* active selected tab in the sidebarmenu */
                                .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                                background-color: #808080;
                                }                          
                                /* other links in the sidebarmenu */
                                .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                                background-color: #808080;
                                }
                                /* other links in the sidebarmenu when hovered */
                                .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                                background-color: #808080;
                                }
                                /* toggle button when hovered  */                    
                                .skin-blue .main-header .navbar .sidebar-toggle:hover{
                                background-color: #870052;
                                }
                                .main-header { max-height: 100px; 
                                font-size:24px; 
                                font-weight:bold; 
                                line-height:24px;
                                }
                                .main-header .logo {
                                height: 100px;
                                font-size:24px; 
                                font-weight:bold; 
                                line-height:24px;align:
                                }
                                .skin-blue .sidebar a {
                                 color: #444;
                                 }
                                .main-sidebar {
                                float:top; margin-top:40px; padding-left:15px; padding-right:15px
                                }
                                '))),
    
    #result section      
    tabsetPanel(type = "tabs",
                tabPanel('Trajectory analysis',
                         h1(""),
                         
                         plotOutput('scatter.plot') %>% withSpinner(),
                         
                         h1(""),
                         
                         uiOutput('download.scatter.data'),
                         uiOutput('download.scatter.plot.pdf'),
                         uiOutput('download.scatter.plot.svg'),
                         
                         h1(""),
                         h1(""),
                         h1(""),
                         
                         plotOutput('sankey.plot') %>% withSpinner(),
                         
                         h1(""),
                         
                         uiOutput('download.sankey.data'),
                         uiOutput('download.sankey.plot.pdf'),
                         uiOutput('download.sankey.plot.svg'),
                         
                         h1(""),
                         h1(""),
                         h1(""),
                         
                         plotOutput('bar.plot') %>% withSpinner(),
                         
                         h1(""),
                         
                         uiOutput('download.bar.data'),
                         uiOutput('download.bar.plot.pdf'),
                         uiOutput('download.bar.plot.svg')),
                
                tabPanel('Sanky data',
                         DT::dataTableOutput("protein.trajectory") %>% withSpinner()),
                tabPanel('Instructions',
                         img(src = 'ProteoTracker_instructions.jpg', height = 1500, align = 'middle')))
  )
)
