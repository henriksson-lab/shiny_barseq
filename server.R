library(plotly)
library(Cairo)
options(shiny.usecairo=T)

if(FALSE){
  #To run this app
  library(shiny)
  runApp(".")
}



server <- function(input, output, session) {

  observeEvent(input$grstats_pool,{
    current_pool <- input$grstats_pool
    grstats <- all_grstats[[current_pool]]
    
    updateSelectizeInput(session, 'grstats_volcano', choices = names(grstats$volcano), server = TRUE)
    updateSelectizeInput(session, 'grstats_scatter', choices = names(grstats$scatterplot), server = TRUE)
    
    grstats <- all_timecourses[[current_pool]]
    updateSelectizeInput(session, 'grstats_gene', choices = grstats$genes, server = TRUE)
    
  })
  


  ################################################################################
  ########### Sample metadata ####################################################
  ################################################################################

  output$plotSamplemetaUmap <- renderPlot(height=700, {

    current_pool <- input$samplemeta_pool
    print(current_pool)
    
    samplemeta <- all_samplemeta[[current_pool]]
    
    ########## Interactive plot
    # fig <- plot_ly(type = 'scatter', mode = 'markers') 
    # fig <- fig %>%
    #   add_trace(
    #     x = samplemeta$umap1, 
    #     y = samplemeta$umap2, 
    #     text= samplemeta$barseqid,
    #     hoverinfo = 'text',
    #     marker = list(color='green'),
    #     showlegend = F
    #   )
    # fig
    
    if(nrow(samplemeta)>0){
    }
    
    p1 <- ggplot(samplemeta, aes(umap1,umap2,color=mouse_ref))+geom_point()
    p2 <- ggplot(samplemeta, aes(umap1,umap2,color=day))+geom_point()
    p3 <- ggplot(samplemeta, aes(umap1,umap2,color=mouse))+geom_point()
    p4 <- ggplot(samplemeta, aes(umap1,umap2,color=is_input))+geom_point()
    p5 <- ggplot(samplemeta, aes(umap1,umap2,color=genotype))+geom_point()
    p6 <- ggplot(samplemeta, aes(umap1,umap2,color=primed))+geom_point()
    #p7 <- ggplot(samplemeta, aes(umap1,umap2,color=pool))+geom_point()
    ptot <- p1/p2|p3/p4|p5/p6 #|p7
    ptot    
    

  })

  ################################################################################
  ########### GRstats - volcano ##################################################
  ################################################################################
  
  
  output$plot_grstats_volcano <- renderPlotly({
    
    current_pool <- input$grstats_pool
    print(current_pool)
    grstats <- all_grstats[[current_pool]]
    thecond <- input$grstats_volcano
    
    if(thecond %in% names(grstats$volcano)){
      toplot <- grstats$volcano[[thecond]]
      theplot <- ggplot(toplot, aes(fc1, p1, label=gene, color=genedesc)) + 
        geom_point(color="gray") + 
        geom_text() +
        xlab(paste("FC",thecond)) + 
        ylab(paste("-log10 pval",thecond))
    } else {
      print("missing volcando cond")
      theplot <- ggplot() + theme_void()
    }
    theplot  %>% ggplotly(source="plot_grstats_volcano") %>% event_register("plotly_click")
  })
  
  
  observeEvent(
    eventExpr = event_data("plotly_click", source = "plot_grstats_volcano"),
    handlerExpr = {
      event_data <- event_data("plotly_click", source = "plot_grstats_volcano")
      #print(event_data)
      clicked_gene <- get_current_scatter()$gene[event_data$pointNumber+1]  #plotly seems to do 0-indexing
      updateSelectInput(session, "grstats_gene", selected = clicked_gene)
    }
  )  
  
  
  ################################################################################
  ########### GRstats - scatter ##################################################
  ################################################################################

  get_current_scatter <- function(){
    current_pool <- input$grstats_pool
    grstats <- all_grstats[[current_pool]]
    thecond <- input$grstats_scatter
    
    if(thecond %in% names(grstats$scatterplot)){
      grstats$scatterplot[[thecond]]
    } else {
      data.frame()
    }
  }
  
  
  
  output$plot_grstats_scatterplot <- renderPlotly({
    
    current_pool <- input$grstats_pool
    grstats <- all_grstats[[current_pool]]
    thecond <- input$grstats_scatter
    
    if(thecond %in% names(grstats$scatterplot)){
      toplot <- grstats$scatterplot[[thecond]]
      thecond2 <- str_split_fixed(thecond," ",2)
      cond1 <- thecond2[1]
      cond2 <- thecond2[2]
      
      fc_range <- range(c(toplot$fc1, toplot$fc2))
      p_range <- range(c(toplot$p1, toplot$p))
      
      theplot <- ggplot(toplot, aes(fc1,fc2, label=gene, color=genedesc)) + geom_point(color="gray") + geom_text()+#size=1) +
            xlab(paste("FC",cond1)) + ylab(paste("FC",cond2)) +
            xlim(fc_range[1], fc_range[2]) + ylim(fc_range[1], fc_range[2])
    } else {
      print("missing scatter cond")
      theplot <- ggplot() + theme_void()
    }
    theplot %>% ggplotly(source="plot_grstats_scatterplot") %>% event_register("plotly_click")
  })
  
  
  
  observeEvent(
    eventExpr = event_data("plotly_click", source = "plot_grstats_scatterplot"),
    handlerExpr = {
      event_data <- event_data("plotly_click", source = "plot_grstats_scatterplot")
      #print(event_data)
      clicked_gene <- get_current_scatter()$gene[event_data$pointNumber+1]  #plotly seems to do 0-indexing
      updateSelectInput(session, "grstats_gene", selected = clicked_gene)
      }
  )  
  
  

  ################################################################################
  ########### GRstats - timecourse ###############################################
  ################################################################################
  
  output$plot_grstats_tcplot <- renderPlotly({
    
    current_pool <- input$grstats_pool
    grstats <- all_timecourses[[current_pool]]

    current_gene <- input$grstats_gene
    tctype <- input$grstats_genesummary  
    
    if(current_gene %in% grstats$genes){
      
      if(tctype=="Per condition"){
        #per cond
        all_gr_meta_percond <- grstats$all_gr_meta_percond
        ggplotly(ggplot(all_gr_meta_percond[all_gr_meta_percond$gene==current_gene,],
                        aes(x=day_int,y=rgr, group=condition_ref, color=condition_ref)) + 
                   geom_line()+
                   xlab("Day")+
                   ggtitle(current_gene))
      } else { #"Per mouse"
        #per mouse      
        all_gr_meta <- grstats$all_gr_meta
        ggplotly(ggplot(all_gr_meta[all_gr_meta$gene==current_gene,],aes(x=day_int,y=rgr, group=mouse_ref, color=mouse_ref)) +
          geom_line()+
            xlab("Day")+
            ggtitle(current_gene))
      }

    } else {
      print(paste("gene not in tc",current_gene))
      ggplotly(ggplot() + theme_void())
    }

  })
  
    
}

