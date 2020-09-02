server <- function(input, output, session){

  cell.lines <- c('' ,'hFF', 'RKO', 'iPSC', 'EB', 'ESC')
  
  comparison.df <- reactive({
    df <- read_tsv('data/finding_significant_proteins_fisher.tsv')
    return(df)
  })
  
  all.data.df <- reactive({
    df <- read_tsv('data/all_data.tsv')
    return(df)
  })
  
  output$filterGenes <- renderUI({
    dat <- comparison.df()
    dat.choices = unique(dat$`Gene names`)
   
    selectizeInput('in.filterGenes', 'Please select here the genes you want to highlight in your analysis',
                   choices = dat.choices, selected = NULL, multiple = T, options = NULL)
  })
  
  #in.t1.denominator
  observe({
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
  
  dummy.plot.df <- eventReactive(input$start.analysis, {
    req(input$in.t1.denominator)
    req(input$in.t1.numerator)
    req(input$in.t2.denominator)
    req(input$in.t2.numerator)
    req(input$in.t3.denominator)
    req(input$in.t3.numerator)
    
    il <- data.frame('x' = c(0, 2, 4, 0),
                     'y' = c(0, 2, 1, 0),
                     'denominator' = c('', '', '', ''),
                     'numerator' = c('', '', '', ''))
    
    il$denominator[1] <- input$in.t1.denominator
    il$numerator[1] <- input$in.t1.numerator
    
    il$denominator[2] <- input$in.t2.denominator
    il$numerator[2] <- input$in.t2.numerator
    
    il$denominator[3] <- input$in.t3.denominator
    il$numerator[3] <- input$in.t3.numerator
    
    il$denominator[4] <- input$in.t1.denominator
    il$numerator[4] <- input$in.t1.numerator
    
    return(il)
  })
  
    
  scatter.df <- eventReactive(c(input$start.analysis,input$in.filterGenes, input$fisher.cutoff), {
    req(input$in.t1.denominator)
    req(input$in.t1.numerator)
    req(input$in.t2.denominator)
    req(input$in.t2.numerator)
    req(input$in.t3.denominator)
    req(input$in.t3.numerator)
  
    dat <- comparison.df()
    
    dat <- dat %>%
      mutate('trajectory.significance' = ifelse(fisher < input$fisher.cutoff, 'passed fisher significance', 'failed fisher significance'))
    
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
  
  scatter.plot <- reactive({
    dat <- scatter.df()
    
    col.palette <- c('#9adbf5', '#1a1818', '#f20c18')
    names(col.palette) = c('failed fisher significance', 'passed fisher significance', 'selected protein')
    
    dat <- dat %>%
      mutate('poi' = as.factor(poi)) %>%
      mutate('highlight' = ifelse(poi == 'poi', 'selected protein', trajectory.significance)) %>%
      mutate('transparency' = ifelse(poi == 'poi', 1, 1))
    
    max.x = max(dat$trajectory.change.stability)
    max.y = max(dat$trajectory.change.expression)
    
    gg <- ggplot(dat, aes(y = trajectory.change.expression, x = trajectory.change.stability, colour = highlight)) +
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
  
  output$scatter.plot <- renderPlot({
    req(scatter.plot())
    scatter.plot()
  })
  
  sankey.protein.df <- eventReactive(c(input$start.analysis,input$in.filterGenes, input$fisher.cutoff), {
    dat <- scatter.df()
    
    if(length(input$in.filterGenes) == 0){
      dat <- dat %>%
        filter(trajectory != 't3') %>%
        mutate('sector' = NA) %>%
        mutate('sector' = ifelse(trajectory.change.stability > 0 & trajectory.change.expression > 0, 'A', sector)) %>%
        mutate('sector' = ifelse(trajectory.change.stability > 0 & trajectory.change.expression < 0, 'B', sector)) %>%
        mutate('sector' = ifelse(trajectory.change.stability < 0 & trajectory.change.expression < 0, 'C', sector)) %>%
        mutate('sector' = ifelse(trajectory.change.stability < 0 & trajectory.change.expression > 0, 'D', sector)) %>%
        mutate('sector' = ifelse(trajectory.significance == 'failed fisher significance', 'E', sector))
      
    }
    
    if(length(input$in.filterGenes) > 0){
      dat <- dat %>%
        filter(trajectory != 't3') %>%
        mutate('sector' = NA) %>%
        mutate('sector' = ifelse(trajectory.change.stability > 0 & trajectory.change.expression > 0, 'A', sector)) %>%
        mutate('sector' = ifelse(trajectory.change.stability > 0 & trajectory.change.expression < 0, 'B', sector)) %>%
        mutate('sector' = ifelse(trajectory.change.stability < 0 & trajectory.change.expression < 0, 'C', sector)) %>%
        mutate('sector' = ifelse(trajectory.change.stability < 0 & trajectory.change.expression > 0, 'D', sector)) %>%
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
  
  output$protein.trajectory <- DT::renderDataTable({
    dat <- sankey.protein.df() %>%
      select(ID:Peptides, trajectory, sector) %>%
      spread(trajectory, sector) %>%
      mutate('transition' = paste0(T1, '_to_', T2))

    return(DT::datatable(dat, extensions = 'Buttons',
                         options = list(dom = 'Bfrtip', scrollX = T, buttons = c('csv', 'excel'), fontSize = '100%'),
                         class = "display"))
  })
  
  sankey.df <- eventReactive(c(input$start.analysis,input$in.filterGenes, input$fisher.cutoff),{
    dat <- sankey.protein.df()
    
    if(length(input$in.filterGenes) == 0){
      s <- dat %>%
        dplyr::select(ID:Peptides, trajectory, sector) %>%
        spread(trajectory, sector) %>%
        unite('transistion', c(T1, T2), remove = F)
      
      s1 <- s %>%
        select(-c(T2, T2)) %>%
        mutate('T' = 'T1')  %>%
        rename('sector' = 'T1')
      
      s2 <- s %>%
        select(-c(T1, T1)) %>%
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
        select(-c(T2, T2)) %>%
        mutate('T' = 'T1')  %>%
        rename('sector' = 'T1')
      
      s2 <- s %>%
        select(-c(T1, T1)) %>%
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
      filter(T == 'T1') %>%
      select(-T)
    
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
    dat <- all.data.df()
  
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
      summarise(ref.mean.value = mean(norm.value)) %>%
      ungroup() %>%
      select(-cell)
    
    dat <- dat %>%
      full_join(dat.ref, by = c('Gene names', 'experiment')) %>%
      mutate('t1.denominator.norm.value' = norm.value - ref.mean.value)
  
    return(dat)
  })
  
  bar.plot <- reactive({
    dat <- bar.df()
    
    dat.ref <- dat %>%
      filter(cell == input$in.t1.denominator) %>%
      select(ID, `Gene names`, Peptides, experiment, replica, t1.denominator.norm.value) %>%
      rename('test.ref' = 't1.denominator.norm.value')
    
    test <- dat %>%
      full_join(dat.ref, by = c('ID', 'Gene names', 'Peptides', 'experiment', 'replica')) %>%
      group_by(ID, `Gene names`, cell, experiment) %>%
      do({
        s <- .
        
        ss <- s %>%
          select(ID, experiment, cell, t1.denominator.norm.value, test.ref) %>%
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
      select(`Gene names`, cell, experiment, asterix) %>%
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
    }

    return(gg1)
  })
  
  output$bar.plot <- renderPlot({
    req(bar.plot())

    legend = data.frame('a' = 'ns  > 0.05    ',
                    'b' = '*  < 0.0    ',
                    'c' = '**  < 0.005    ',
                    'd' = '***  < 0.0005  ')

    g2 <- tableGrob(legend, cols = NULL, rows = NULL)

    g <- arrangeGrob(bar.plot(), g2, heights=c(15, 1))

    return(plot(g))
  })
  
  ## DOWNLOADS
  #scatter plot
  output$download.scatter.data <- renderUI({
    req(scatter.df())
    downloadButton('out.download.scatter.data', 'Download scatter plot data')
  })
  
  output$out.download.scatter.data <- downloadHandler(
    filename = function() {'scatter_plot_data.csv'},
    content = function(file){write.csv(scatter.df(), file, row.names = F)}
  )
  
  output$download.scatter.plot <- renderUI({
    req(scatter.df())
    downloadButton('out.download.scatter.plot', 'Download scatter plot')
  })
  
  output$out.download.scatter.plot <- downloadHandler(
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
    filename = function() {'sankey_plot_data.csv'},
    content = function(file){write.csv(sankey.plot.df(), file, row.names = F)}
  )
  
  output$download.sankey.plot <- renderUI({
    req(sankey.df())
    downloadButton('out.download.sankey.plot', 'Download Sanky plot')
  })
  
  output$out.download.sankey.plot <- downloadHandler(
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
    filename = function() {'barplot_plot_data.csv'},
    content = function(file){write.csv(bar.df(), file, row.names = F)}
  )
  
  output$download.bar.plot <- renderUI({
    req(bar.df())
    downloadButton('out.download.bar.plot', 'Download barplot plot')
  })
  
   output$out.download.bar.plot <- downloadHandler(
    filename = function(){'barplot_plot.pdf'},
    content = function(file){ggsave(file, bar.plot(), device = 'pdf')}
  )
  
  #link to paper
  url <- a('empty',
           href = 'https://www.https://github.com/RZlab')
  output$citation <- renderUI({
    tagList('Please cite:', url)
  })
}
