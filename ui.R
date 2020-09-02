cell.lines <- c('' ,'hFF', 'RKO', 'iPSC', 'EB', 'ESC')

ui <- dashboardPage(
  dashboardHeader(
    title = img(src = 'ProteoTracker.jpg', height = 100, align = 'left'),
    titleWidth = 350
  ),
  
  dashboardSidebar(
    width = 350,
    
    selectInput("in.t1.denominator", "Please specify the denominator for T1", cell.lines, multiple = F),
    selectInput("in.t1.numerator", "Please specify the numerator for T1", cell.lines, multiple = F),
    
    selectInput("in.t2.denominator", "Please specify the denominator for T2", cell.lines, multiple = F),
    selectInput("in.t2.numerator", "Please specify the numerator for T2", cell.lines, multiple = F),
    
    selectInput("in.t3.denominator", "Please specify the denominator for T3", cell.lines, multiple = F),
    selectInput("in.t3.numerator", "Please specify the numerator for T3", cell.lines, multiple = F),
    
    sliderInput("fisher.cutoff", "Specify cutoff for fisher significance", min = 0, max = 1, value = 0.05, step = 0.01),
    
    uiOutput('filterGenes'),

    actionButton('start.analysis', 'Start Analysis'),

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
                        uiOutput('download.scatter.plot'),
                        
                        h1(""),
                        h1(""),
                        h1(""),
                        
                        plotOutput('sankey.plot') %>% withSpinner(),
                        
                        h1(""),
                        
                        uiOutput('download.sankey.data'),
                        uiOutput('download.sankey.plot'),
                        
                        h1(""),
                        h1(""),
                        h1(""),
                        
                        plotOutput('bar.plot') %>% withSpinner(),
                        
                        h1(""),
                        
                        uiOutput('download.bar.data'),
                        uiOutput('download.bar.plot')),
                        
                tabPanel('Sanky data',
                        dataTableOutput("protein.trajectory") %>% withSpinner()))
    

  )
)