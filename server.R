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

calculate.pval.and.fisher <- function(dat = NULL){
  cell.lines = unique(dat$cell)
  res <- list()
  n = 1
  for(i in 1:length(cell.lines)){
    for(j in 1:length(cell.lines)){
      if (j == i){next}
      denominator = c(cell.lines[i])
      numerator = c(cell.lines[j])
      name = paste0(numerator, '_vs_', denominator)
      
      test.Exp <- dat %>%
        filter(experiment == 'Exp') %>%
        filter(cell == denominator | cell == numerator) %>%
        mutate('category' = ifelse(cell == denominator, 'denominator', 'numerator')) %>%
        group_by(`ID`, `Gene names`, Peptides) %>%
        do({
          s <- .
          p <- try(t.test(value ~ category, data = s, var.equal = T), silent = T)
          tidy(p)
        }) %>%
        mutate('trajectory.change.Exp' = estimate2 - estimate1) %>%
        dplyr::select(ID:Peptides, trajectory.change.Exp, p.value) %>%
        rename('p.value.trajectory.change.Exp' = 'p.value')
      
      test.Sm <- dat %>%
        filter(experiment == 'Sm') %>%
        filter(cell == denominator | cell == numerator) %>%
        mutate('category' = ifelse(cell == denominator, 'denominator', 'numerator')) %>%
        group_by(`ID`, `Gene names`, Peptides) %>%
        do({
          s <- .
          p <- try(t.test(value ~ category, data = s, var.equal = T), silent = T)
          tidy(p)
        }) %>%
        mutate('trajectory.change.Sm' = estimate2 - estimate1) %>%
        dplyr::select(ID:Peptides, trajectory.change.Sm, p.value) %>%
        rename('p.value.trajectory.change.Sm' = 'p.value')
      
      res[[n]] <- test.Exp %>%
        full_join(test.Sm, by = c('ID', 'Gene names', 'Peptides')) %>%
        mutate('comparison' = name)
      n = n+1
    }
  }
  dat <- plyr::rbind.fill(res)
  
  res <- rep(NA, nrow(dat))
  for(i in 1:nrow(dat)){
    res[i] <-  metap::sumlog(c(dat$p.value.trajectory.change.Exp[i], dat$p.value.trajectory.change.Sm[i]))$p
  }
  dat$fisher.score = res
  
  return(dat)
}


server <- function(input, output, session){
  
  #extending max file size
  options(shiny.maxRequestSize = 50*1024^2)

  #file location for ProteTracker data
  data.set <- eventReactive(c(input$checkGroup, input$upload.file),{
    if(!is.null(input$checkGroup)){
      if(input$checkGroup == 1){
        ds <- 'data/ProteoTracker_data.tsv'
      }
      if(input$checkGroup == 2){
        req(input$upload.file)
        inFile <- input$upload.file
        if(!is.null(inFile)){
          ds <- inFile$datapath
        }
      } 
    }
    
    return(ds)
  })
  
  #read data into shiny
  dat <- reactive({
    ds <- data.set()
    if(is.character(ds)){
      df <- suppressMessages(read_tsv(ds)) %>%
          gather(key, value, -(ID:Peptides)) %>%
          separate(key, sep = '_', into = c('experiment', 'cell', 'replica'))
    }
    else{
      df <- NULL
    }
    
    return(df)
  })

  #determine denominator and numerator
  cell.states <- reactive({
    o <- NULL
    if(!is.null(dat())){
      d <- dat()
      o <- append('', unique(d$cell))
    }
    else{
      o <- NULL
    }
  
    return(o)
  })
  
  #determine all genes in
  output$filterGenes <- renderUI({
    dat <- dat()
    dat.choices = unique(dat$`Gene names`)
    
    selectizeInput('in.filterGenes', 'Please select here the genes you want to highlight in your analysis',
                   choices = dat.choices, selected = NULL, multiple = T, options = NULL)
  })
  
  #add start button
  output$action <- renderUI({
    req(dat())
    
    actionButton('start.analysis', 'Start Analysis')
  })
  
  #in.t1.denominator
  observe({
    cell.lines <- cell.states()
    x <- input$in.t1.denominator
    x.exclude1 <- input$in.t1.numerator
    x.exclude2 <- input$in.t2.numerator
    x.exclud <- append(x.exclude1, x.exclude2)
    cell.lines.out <- cell.lines[!cell.lines %in% x.exclud]
    
    updateSelectInput(session, "in.t3.denominator",
                      choices = cell.lines.out,
                      selected = x
    )
  })
  
  #in.t1.numerator
  observe({
    cell.lines <- cell.states()
    x <- input$in.t1.numerator
    x.exclude1 <- input$in.t1.denominator
    x.exclude2 <- input$in.t2.numerator
    x.exclud <- append(x.exclude1, x.exclude2)
    cell.lines.out <- cell.lines[!cell.lines %in% x.exclud]
    
    updateSelectInput(session, "in.t2.denominator",
                      choices = cell.lines.out,
                      selected = x
    )
  })
  
  #in.t2.denominator
  observe({
    cell.lines <- cell.states()
    x <- input$in.t2.denominator
    x.exclude1 <- input$in.t1.denominator
    x.exclude2 <- input$in.t2.numerator
    x.exclud <- append(x.exclude1, x.exclude2)
    cell.lines.out <- cell.lines[!cell.lines %in% x.exclud]
    
    updateSelectInput(session, "in.t1.numerator",
                      choices = cell.lines.out,
                      selected = x
    )
  })
  
  #in.t2.numerator
  observe({
    cell.lines <- cell.states()
    x <- input$in.t2.numerator
    x.exclude1 <- input$in.t1.denominator
    x.exclude2 <- input$in.t1.numerator
    x.exclud <- append(x.exclude1, x.exclude2)
    cell.lines.out <- cell.lines[!cell.lines %in% x.exclud]
    
    updateSelectInput(session, "in.t3.numerator",
                      choices = cell.lines.out,
                      selected = x
    )
  })
  
  #in.t3.denominator
  observe({
    cell.lines <- cell.states()
    x <- input$in.t3.denominator
    x.exclude1 <- input$in.t1.numerator
    x.exclude2 <- input$in.t2.numerator
    x.exclud <- append(x.exclude1, x.exclude2)
    cell.lines.out <- cell.lines[!cell.lines %in% x.exclud]
    
    updateSelectInput(session, "in.t1.denominator",
                      choices = cell.lines.out,
                      selected = x
    )
  })
  
  #in.t3.dnumerator
  observe({
    cell.lines <- cell.states()
    x <- input$in.t3.numerator
    x.exclude1 <- input$in.t1.denominator
    x.exclude2 <- input$in.t1.numerator
    x.exclud <- append(x.exclude1, x.exclude2)
    cell.lines.out <- cell.lines[!cell.lines %in% x.exclud]
    
    updateSelectInput(session, "in.t2.numerator",
                      choices = cell.lines.out,
                      selected = x
    )
  })
  
  #popup message
  observeEvent(input$start.analysis,{
    if(input$checkGroup == 2){
      showModal(modalDialog(
        title = '',
        'Your data is getting processed. Please have patience! \n
        Depending on the size of your input data this could take up to a few minutes',
        easyClose = T
      ))
    }
  })
  
  #calculate pvale and fisher score of all possible combinations
  calculate.dat <- reactive({
    if(!is.null(dat())){
      dat <- dat()
      if(input$checkGroup == 1){
        dat <- suppressMessages(read_tsv('data/ProteoTracker_data_processed.tsv'))
      }
      if(input$checkGroup == 2){
        dat <- calculate.pval.and.fisher(dat = dat)
      }
    }
    else{
      dat <- NULL
    }
    
    return(dat)    
  })
  
  #make fisher score cutoff criteria
  scatter.df <- eventReactive(c(input$start.analysis, input$in.filterGenes, input$fisher.cutoff,
                                input$in.t1.denominator, input$in.t1.numerator, input$in.t2.denominator,
                                input$in.t2.numerator, input$in.t3.denominator, input$in.t3.numerator), {
    req(input$in.t1.denominator)
    req(input$in.t1.numerator)
    req(input$in.t2.denominator)
    req(input$in.t2.numerator)
    req(input$in.t3.denominator)
    req(input$in.t3.numerator)
    req(input$start.analysis)
    
    dat <- calculate.dat()
    
    dat <- dat %>%
      mutate('trajectory.significance' = ifelse(fisher.score < input$fisher.cutoff, 'passed fisher significance', 'failed fisher significance'))
    
    comparisons = c(paste0(input$in.t1.numerator, '_vs_', input$in.t1.denominator), 
                    paste0(input$in.t2.numerator, '_vs_', input$in.t2.denominator),
                    paste0(input$in.t3.numerator, '_vs_', input$in.t3.denominator))
    
    res <- list()
    for(i in 1:length(comparisons)){
      res[[i]] <- dat %>%
        filter(comparison == comparisons[i]) %>%
        mutate('trajectory' = paste0('T', i))
    }
    dat <- plyr::rbind.fill(res)
    
    dat$poi <- 'bg'
    if(length(input$in.filterGenes) > 0){
      x  <- input$in.filterGenes
      for(i in 1:length(x)){
        w <- which(dat$`Gene names` == input$in.filterGenes[i])
        dat$poi[w] <- 'poi'
      }
    }
    
    return(dat)
  })
  
  #prepare data for scatter plot
  scatter.plot <- reactive({
    dat <- scatter.df()
    
    col.palette <- c('#9adbf5', '#1a1818', '#f20c18')
    names(col.palette) = c('failed fisher significance', 'passed fisher significance', 'selected protein')
    
    dat <- dat %>%
      mutate('poi' = as.factor(poi)) %>%
      mutate('highlight' = ifelse(poi == 'poi', 'selected protein', trajectory.significance)) %>%
      mutate('transparency' = ifelse(poi == 'poi', 1, 1))
    
    max.x = max(dat$trajectory.change.Sm)
    max.y = max(dat$trajectory.change.Exp)
    
    gg <- ggplot(dat, aes(y = trajectory.change.Exp, x = trajectory.change.Sm, colour = highlight)) +
      geom_point(aes(alpha = transparency), shape = 16) +
      geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.5) +
      theme_bw() +
      xlim(-max.x, max.x) +
      ylim(-max.y, max.y) +
      theme(legend.position = 'bottom') +
      geom_point(data = subset(dat, poi == 'poi'), aes(alpha = transparency), shape = 16) + 
      scale_color_manual(values = col.palette) +
      facet_wrap(~trajectory + comparison, nrow = 1) +
      geom_text_repel(data = subset(dat, highlight == 'selected protein'),
                      aes(label = `Gene names`), colour = 'red', size = 4, box.padding = unit(0.35, 'lines'),
                      point.padding = unit(0.3, 'lines')) +
      guides(alpha = F) +
      labs(fill = 'legend') +
      xlab('log2 mean Sm') +
      ylab('log2 mean Exp')
    
    return(gg)
  })
  
  #print data for scatter plot
  output$scatter.plot <- renderPlot({
    req(scatter.plot())
    scatter.plot()
  })
  
  #prepare data for sankey plot
  sankey.protein.df <- eventReactive(c(input$start.analysis, input$in.filterGenes, input$fisher.cutoff,
                                       input$in.t1.denominator, input$in.t1.numerator, input$in.t2.denominator,
                                       input$in.t2.numerator, input$in.t3.denominator, input$in.t3.numerator), {
    dat <- scatter.df()
    if(length(input$in.filterGenes) == 0){
      dat <- dat %>%
        filter(trajectory != 't3') %>%
        mutate('sector' = NA) %>%
        mutate('sector' = ifelse(trajectory.change.Sm > 0 & trajectory.change.Exp > 0, 'A', sector)) %>%
        mutate('sector' = ifelse(trajectory.change.Sm > 0 & trajectory.change.Exp < 0, 'B', sector)) %>%
        mutate('sector' = ifelse(trajectory.change.Sm < 0 & trajectory.change.Exp < 0, 'C', sector)) %>%
        mutate('sector' = ifelse(trajectory.change.Sm < 0 & trajectory.change.Exp > 0, 'D', sector)) %>%
        mutate('sector' = ifelse(trajectory.significance == 'failed fisher significance', 'E', sector))
      
    }
    
    if(length(input$in.filterGenes) > 0){
      dat <- dat %>%
        filter(trajectory != 't3') %>%
        mutate('sector' = NA) %>%
        mutate('sector' = ifelse(trajectory.change.Sm > 0 & trajectory.change.Exp > 0, 'A', sector)) %>%
        mutate('sector' = ifelse(trajectory.change.Sm > 0 & trajectory.change.Exp < 0, 'B', sector)) %>%
        mutate('sector' = ifelse(trajectory.change.Sm < 0 & trajectory.change.Exp < 0, 'C', sector)) %>%
        mutate('sector' = ifelse(trajectory.change.Sm < 0 & trajectory.change.Exp > 0, 'D', sector)) %>%
        mutate('sector' = ifelse(trajectory.significance == 'failed fisher significance', 'E', sector))
      
      dat$poi <- 0
      x  <- input$in.filterGenes
      for(i in 1:length(x)){
        w <- which(dat$`Gene names` == input$in.filterGenes[i])
        dat$poi[w] <- 1
      }
      dat$poi <- as.factor(dat$poi)
    }
    
    return(dat)
  })
  
  output$protein.trajectory <- DT::renderDataTable(server = T,{
    dat <- sankey.protein.df() %>%
      dplyr::select(ID:Peptides, trajectory, sector) %>%
      spread(trajectory, sector) %>%
      mutate('transition' = paste0(T1, '_to_', T2))
    
    protein.trajectory <- DT::datatable(data = dat,
                                        extensions = 'Buttons',
                                        options = list(dom = 'Bfrtip',
                                                       scrollX = T,
                                                       pageLength = 50,
                                                       buttons = list(list(extend = 'csv', filename = 'protein_trajectory'), 
                                                                      list(extend = 'excel', filename = 'protein_trajectory')),
                                                       fontSize = '100%'),
                                        class = "display")
    
    return(protein.trajectory)
  })
  
  
  sankey.df <- eventReactive(c(input$start.analysis,input$in.filterGenes, input$fisher.cutoff,
                               input$in.t1.denominator, input$in.t1.numerator, input$in.t2.denominator,
                               input$in.t2.numerator, input$in.t3.denominator, input$in.t3.numerator),{
    dat <- sankey.protein.df()
    
    if(length(input$in.filterGenes) == 0){
      s <- dat %>%
        dplyr::select(ID:Peptides, trajectory, sector) %>%
        spread(trajectory, sector) %>%
        unite('transistion', c(T1, T2), remove = F)
      
      s1 <- s %>%
        dplyr::select(-c(T2, T3)) %>%
        mutate('T' = 'T1')  %>%
        rename('sector' = 'T1')
      
      s2 <- s %>%
        dplyr::select(-c(T1, T3)) %>%
        mutate('T' = 'T2') %>%
        rename('sector' = 'T2')
      
      ss1 <- s1 %>%
        group_by(T, sector, transistion) %>%
        summarise(n = n())
      
      ss2 <- s2 %>%
        group_by(T, sector, transistion) %>%
        summarise(n = n())
      
      x <- rbind(ss1, ss2)
    }
    
    if(length(input$in.filterGenes) > 0){
      s <- dat %>%
        dplyr::select(ID:Peptides, poi, trajectory, sector) %>%
        spread(trajectory, sector) %>%
        mutate('transistion' = ifelse(poi == 1, paste(T1, T2, `Gene names`, sep = '_'), paste(T1, T2, 'NULL', sep = '_')))
      
      s1 <- s %>%
        dplyr::select(-c(T2, T3)) %>%
        mutate('T' = 'T1')  %>%
        rename('sector' = 'T1')
      
      s2 <- s %>%
        dplyr::select(-c(T1, T3)) %>%
        mutate('T' = 'T2') %>%
        rename('sector' = 'T2')
      
      ss1 <- s1 %>%
        group_by(T, sector, transistion) %>%
        summarise(n = n())
      
      ss2 <- s2 %>%
        group_by(T, sector, transistion) %>%
        summarise(n = n())
      
      x <- rbind(ss1, ss2) %>%
        separate(transistion, sep = '_', into = c('T1', 'T2', 'poi'), remove = F) %>%
        mutate('highlight' = ifelse(poi != 'NULL', poi, 'bg')) %>%
        mutate('sector.new' = ifelse(T1 == 'T1', paste0(sector, '_', poi), paste0(sector, '_', poi))) %>%
        mutate('transparency' = ifelse(poi != 'NULL', 1, 0.4)) 
      
    }
    
    return(x)
  })
  
  sankey.plot <- reactive({
    dat <- sankey.df()
    
    if(ncol(dat) == 4){
      gg <- ggplot(dat,
                   aes(x = T, stratum = sector, alluvium = transistion,
                       y = n,
                       fill = sector, label = sector)) +
        scale_fill_brewer(type = 'qual', palette = 'Set2') +
        geom_flow(stat = 'alluvium', lode.guidance = 'frontback',
                  color = 'darkgray') +
        geom_stratum(alpha = .5) +
        geom_text(stat = "stratum", size = 3) +
        theme_bw() +
        theme(legend.position = 'none') +
        xlab('Proteome transition') +
        ylab('Number of proteins')
    }
    
    if(ncol(dat) > 4){
      gg <- ggplot(dat,
                   aes(x = T,
                       y = n,
                       stratum = sector.new,
                       alluvium = transistion)) +
        geom_text_repel(data = subset(dat, T == 'T2'), aes(label = ifelse(highlight != 'bg', highlight, NA), colour = sector),
                        stat = "stratum", size = 4, direction = "y", nudge_x = 0.5) +
        geom_text_repel(data = subset(dat, T == 'T1'), aes(label = ifelse(highlight != 'bg', highlight, NA), colour = sector),
                        stat = "stratum", size = 4, direction = "y", nudge_x = -0.5) +
        stat_stratum(aes(fill = sector, stratum = sector), alpha = 0.5) +
        geom_text(data = subset(dat, highlight == 'bg'), aes(label = sector), stat = "stratum", size = 3) +
        geom_flow(stat = 'alluvium', aes(colour = sector, alpha = transparency)) +
        theme_bw() +
        theme(legend.position = 'none') +
        scale_fill_brewer(type = 'qual', palette = 'Set2') +
        scale_colour_brewer(type = 'qual', palette = 'Set2') +
        xlab('Proteome transition') +
        ylab('Number of proteins')
    }
    
    return(gg)
  })
  
  output$sankey.plot <- renderPlot({
    req(sankey.plot())
    sankey.plot()
  })
  
  sankey.plot.df <- reactive({
    req(sankey.df())
    
    dat <- sankey.df() %>%
      ungroup() %>%
      filter('T' == 'T1') %>%
      dplyr::select(-`T`)
    
    s <- dat %>%
      group_by(sector) %>%
      summarise(sum.sector = sum(n)) %>%
      rename('n.sector' = 'sum.sector')
    
    dat <- dat %>%
      full_join(s, by = 'sector') %>%
      mutate('proportion.to.sector' = round(n/n.sector, 4))
    
    return(dat)
  })
  
  bar.df <- eventReactive(input$in.filterGenes, {
    req(input$in.filterGenes)
    dat <- dat()

    dat$poi <- 0
    x <- input$in.filterGenes
    for(i in 1:length(x)){
      w <- which(dat$`Gene names` == input$in.filterGenes[i])
      dat$poi[w] <- 1
    }
    
    dat <- dat %>%
      filter(poi == 1)
    
    dat.ref <- dat %>%
      filter(cell == input$in.t1.denominator) %>%
      group_by(`Gene names`, cell, experiment) %>%
      summarise(ref.mean.value = mean(value)) %>%
      ungroup() %>%
      dplyr::select(-cell)
    
    dat <- dat %>%
      full_join(dat.ref, by = c('Gene names', 'experiment')) %>%
      mutate('t1.denominator.norm.value' = value - ref.mean.value)
    
    return(dat)
  })
  
  bar.plot <- reactive({
    dat <- bar.df()
    
    dat.ref <- dat %>%
      filter(cell == input$in.t1.denominator) %>%
      dplyr::select(ID, `Gene names`, Peptides, experiment, replica, t1.denominator.norm.value) %>%
      rename('test.ref' = 't1.denominator.norm.value')
    
    test <- dat %>%
      full_join(dat.ref, by = c('ID', 'Gene names', 'Peptides', 'experiment', 'replica')) %>%
      group_by(ID, `Gene names`, cell, experiment) %>%
      do({
        s <- .
        
        ss <- s %>%
          dplyr::select(ID, experiment, cell, t1.denominator.norm.value, test.ref) %>%
          gather(key, value, t1.denominator.norm.value:test.ref)
        
        p <- try(t.test(value ~ key, data = ss, var.equal = T), silent = T)
        
        tidy(p)
      }) %>%
      mutate('asterix' = 'ns') %>%
      mutate('asterix' = ifelse(p.value < 0.05, '*', asterix)) %>%
      mutate('asterix' = ifelse(p.value < 0.005, '**', asterix)) %>%
      mutate('asterix' = ifelse(p.value < 0.0005, '***', asterix)) 
    
    t <- test %>%
      ungroup() %>%
      dplyr::select(`Gene names`, cell, experiment, asterix) %>%
      filter(cell != input$in.t1.denominator)
    
    s <- dat %>%
      group_by(`Gene names`, cell, experiment) %>%
      summarise('mean.value' = mean(t1.denominator.norm.value),
                'sd.value' = sd(t1.denominator.norm.value)) %>%
      full_join(t, by = c('Gene names', 'cell', 'experiment')) %>%
      mutate('asterix.y' = ifelse(mean.value > 0, (mean.value+sd.value)+0.2, (mean.value-sd.value)-0.2))
    
    if(nrow(dat) > 0){
      gg1 <- ggplot(s, 
                    aes(x = cell, y = mean.value, fill = experiment)) +
        geom_bar(stat = 'identity', position = 'dodge') +
        geom_col(position = position_dodge(0.9)) +
        geom_errorbar(aes(ymin = mean.value - sd.value, ymax = mean.value + sd.value), width = 0.25, position = position_dodge(0.9)) +
        geom_hline(yintercept = 0, col = 'black', linetype = 'solid', alpha = 0.5) +
        geom_text(aes(y = asterix.y, label = asterix), position = position_dodge(0.9)) +
        facet_wrap(~`Gene names`) +
        scale_fill_manual(values = c('#f2bde1', '#f7bf5e')) +
        theme_bw() +
        theme(legend.position = 'bottom') +
        xlab('Cell line') +
        ylab('Foldchange in stability and expression
         compared to denominator in T1 +/- SD [log2(Fc)]')
      
      legend = data.frame('a' = 'ns  > 0.05    ',
                          'b' = '*  < 0.0    ',
                          'c' = '**  < 0.005    ',
                          'd' = '***  < 0.0005  ')
      
      gg2 <- tableGrob(legend, cols = NULL, rows = NULL, theme = ttheme_minimal())
      
      gg <- arrangeGrob(gg1, gg2, heights=c(15, 1))
      gg <- ggdraw(gg) + 
        theme(plot.background = element_rect(fill="white", color = NA))
    }
    
    return(gg)
  })
  
  output$bar.plot <- renderPlot({
    req(bar.plot())
    bar.plot()
  })
  
  ## DOWNLOADS
  #scatter plot
  output$download.scatter.data <- renderUI({
    req(scatter.df())
    downloadButton('out.download.scatter.data', 'Download scatter plot data')
  })
  
  output$out.download.scatter.data <- downloadHandler(
    filename = function() {'scatter_plot_data.tsv'},
    content = function(file){write.table(scatter.df(), file, row.names = F)
      })
  
  output$download.scatter.plot.svg <- renderUI({
    req(scatter.df())
    downloadButton('out.download.scatter.plot.svg', 'Download scatter plot as svg')
  })
  output$download.scatter.plot.pdf <- renderUI({
    req(scatter.df())
    downloadButton('out.download.scatter.plot.pdf', 'Download scatter plot as pdf')
  })
  
  output$out.download.scatter.plot.svg <- downloadHandler(
    filename = function(){'scatter_plot.svg'},
    content = function(file){ggsave(file, scatter.plot(), device = 'svg')
    }
  )
  output$out.download.scatter.plot.pdf <- downloadHandler(
    filename = function(){'scatter_plot.pdf'},
    content = function(file){ggsave(file, scatter.plot(), device = 'pdf')
    }
  )
  
  #sankey stuff
  output$download.sankey.data <- renderUI({
    req(sankey.plot.df())
    downloadButton('out.download.sankey.data', 'Download Sankey data')
  })
  
  output$out.download.sankey.data <- downloadHandler(
    filename = function() {'sankey_plot_data.tsv'},
    content = function(file){write.table(sankey.plot.df(), file, row.names = F)}
  )
  
  output$download.sankey.plot.svg <- renderUI({
    req(sankey.df())
    downloadButton('out.download.sankey.plot.svg', 'Download Sanky plot as svg')
  })
  output$download.sankey.plot.pdf <- renderUI({
    req(sankey.df())
    downloadButton('out.download.sankey.plot.pdf', 'Download Sanky plot as pdf')
  })
  
  output$out.download.sankey.plot.svg <- downloadHandler(
    filename = function(){'sankey_plot.svg'},
    content = function(file){ggsave(file, sankey.plot(), device = 'svg')
    }
  )
  output$out.download.sankey.plot.pdf <- downloadHandler(
    filename = function(){'sankey_plot.pdf'},
    content = function(file){ggsave(file, sankey.plot(), device = 'pdf')
    }
  )
  
  #barplot stuff
  output$download.bar.data <- renderUI({
    req(bar.df())
    downloadButton('out.download.bar.data', 'Download barplot data')
  })
  
  output$out.download.bar.data <- downloadHandler(
    filename = function() {'barplot_plot_data.tsv'},
    content = function(file){write.table(bar.df(), file, row.names = F)}
  )
  
  output$download.bar.plot.svg <- renderUI({
    req(bar.df())
    downloadButton('out.download.bar.plot.svg', 'Download barplot plot as svg')
  })
  output$download.bar.plot.pdf <- renderUI({
    req(bar.df())
    downloadButton('out.download.bar.plot.pdf', 'Download barplot plot as pdf')
  })
  
  output$out.download.bar.plot.svg <- downloadHandler(
    filename = function(){'barplot_plot.svg'},
    content = function(file){ggsave(file, bar.plot(), device = c('svg'))}
  )
  output$out.download.bar.plot.pdf <- downloadHandler(
    filename = function(){'barplot_plot.pdf'},
    content = function(file){ggsave(file, bar.plot(), device = c('pdf'))}
  )
  
  #link to paper
  url <- a('empty',
           href = 'github.com/RZlab')
  output$citation <- renderUI({
    tagList('Please cite:', url)
  })
}