library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)
library(magrittr)
source("qrmr_functions.R")
qrmrf_merge <- data.table::fread("data/qrmrf_merge.tsv")
linear <- data.table::fread("data/linear_merge_sub.tsv")

if(FALSE){
    input <- list()
    input$b0p_range <- c(0,306)
    input$b1p_range <- c(0,307)
    input$b2p_range <- c(0,303)
    input$b0_range <- c(-17,8)
    input$b1_range <- c(-18,26)
    input$b2_range <- c(-80,73)
    input$phenotype = 'a1c'
    
}

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    output$formula <- renderUI({
        withMathJax(
            paste0("$$y=\\beta_0+\\beta_1\\tau+\\beta_2(0.5-\\tau)^2+e$$"))
    })
    
    flag_qrmrf <- reactive({
        q1 <- qrmrf_merge %>% dplyr::mutate(
            match_filter = ifelse(
                MR_b0_log10p > {input$b0_range[1]} &
                    MR_b0_log10p < {input$b0p_range[2]} &
                    MR_b1_log10p > {input$b1p_range[1]} &
                    MR_b1_log10p < {input$b1p_range[2]} &
                    MR_b2_log10p > {input$b2p_range[1]} & 
                    MR_b2_log10p < {input$b2p_range[2]} &
                    MR_b0 > {input$b0_range[1]} &
                    MR_b0 < {input$b0_range[2]} &
                    MR_b1 > {input$b1_range[1]} &
                    MR_b1 < {input$b1_range[2]} &
                    MR_b2 > {input$b2_range[1]} &
                    MR_b2 < {input$b2_range[2]} &
                    phenotype == input$phenotype,
                             TRUE, FALSE) %>% 
                na.omit
        )
    })
    
    react_qrmrf <- reactive({
        q2 <- flag_qrmrf()[flag_qrmrf()$match_filter == TRUE,]
        q2 <- q2[order(q2$min_pval_12), ]
    })

    
    react_rsid <- reactive({
        r1 <- react_qrmrf()$rsid[react_qrmrf()$phenotype == input$phenotype]
    })
    output$choose_rsid <- renderUI({
        selectInput("rsid", "Choose RSID:", 
                    as.list({react_rsid()}), selected = {react_rsid()}[1])
    })

    output$table <- renderTable({react_qrmrf()}[
        ,c("chr", "rsid", "pos", "MR_b0", "MR_b1", "MR_b2", 
           "MR_b0_log10p", "MR_b1_log10p", "MR_b2_log10p", "log10_min_pval")])
    output$length <- renderText(length({react_rsid()}))

    output$beta_scatter <- renderPlotly({
        g1 <- gg_mr_beta_scatter({flag_qrmrf()}) +
            facet_grid(.~phenotype) +
            ggtitle(paste0("QRMR coefficient scatterplot"))
        g1  %>% ggplotly
    })

    output$g_manhattan <- renderPlotly({
        gg_manhattan({react_qrmrf()})
    })

    output$g_qrmr <- renderPlot({
        this_linear <- linear %>% dplyr::filter(phenotype == input$phenotype)
        gg_qrmr({react_qrmrf()}, this_linear, input$rsid)
    })


})
