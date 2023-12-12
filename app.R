library(shiny)
library(shinyBS)
library(shinydashboard)
library(shinyjs)
source("plots.R")
# source("pseudo.R")
# source("DEGs.R")

library(crayon)
message(red("loading data ... "))
# load("~/Downloads/multiome_data/scelist_azi_sr.RData")
# scelist <- readRDS("~/Downloads/multiome_data/scelistwt_shiny_14_June.RDS")
# save(scelist, mnnsce, file="./data/data_shiny_June23.RData")
# load("/Users/inzirio/Downloads/multiome_data/SleepMultiome_MNN.RData")
# mnnsce <- mnn.out
# saveRDS(mnnsce, file="/Users/inzirio/Downloads/multiome_data/SleepMultiome_MNN.Rds")
# mnnsce <- readRDS(file="/Users/inzirio/Downloads/multiome_data/SleepMultiome_MNN.Rds")
load("~/Downloads/multiome_data/data_shiny_Dec23.RData")
# names(assays(mnnsce)) <- "logcounts"
# mnnscelist <- list(WT=mnnsce[, mnnsce$condition=="WT"], SD=mnnsce[, mnnsce$condition=="SD"])
# pseudo <- readRDS("~/Downloads/multiome_data/GEX_pseudo_filtered_edger_5_0.RDS")

# pseudo <- pseudo[rowData(pseudo)$edgeR_5_0,]
# edgeRres <- readRDS("~/Downloads/multiome_data/edgeres_ruvs_1_6.RDS")
message(red("... Done!"))

ui <- navbarPage("Sleep Multiome",
    tabPanel("GEX investigation\nBy Sample",
        sidebarLayout(
            sidebarPanel(
                useShinyjs(),
                selectInput(inputId="sample_id", label="Sample", 
                    choices=names(scelist)),
                selectizeInput(inputId="gene_id", label="Gene",
                    choices=NULL),
                selectizeInput(inputId="ct_id", label="Cell Type",
                               choices=NULL),
                # tipify(checkboxInput(inputId="mnn_id", label="MNN"),
                #        "Check to show a unified plot for all samples")
            ),
            mainPanel(
               plotOutput(outputId="plot_id_out")
            )
        )
    ),
    tabPanel("GEX investigation\nBy Condition",
        sidebarLayout(
            sidebarPanel(
                useShinyjs(),
                selectInput(inputId="cond_id", label="Condition", 
                    choices=names(mnnscelist)),
                selectizeInput(inputId="gene_id1", label="Gene",
                    choices=NULL),
                selectizeInput(inputId="ct_id1", label="Cell Type",
                    choices=NULL)
                ),
                mainPanel(
                    plotOutput(outputId="plot_id_out1")
                )
        )
    )
    # tabPanel("GEX Pseudo-Bulk Norm",
    #     sidebarLayout(
    #         sidebarPanel(
    #             selectizeInput(inputId="ct_id", label="Cell Type",
    #                 choices=NULL),
    #             numericInput(inputId="ruvk_id", label="RUV k", value=1, step=1),
    #             checkboxInput(inputId="zoom_k_id", label="RUV k zoom")
    #         ),
    #         mainPanel(
    #             plotOutput(outputId="plot_pseudo_out")
    #         )
    #     )
    # ), 
    # tabPanel("GEX Pseudo-Bulk DEGs",
    #     sidebarLayout(
    #         sidebarPanel(
    #             selectizeInput(inputId="ct_id2", label="Cell Type", choices=NULL),
    #             radioButtons(inputId="ruvsg_id", label="RUVs/RUVg", choices=c("RUVs", "RUVg")),
    #             numericInput(inputId="ruvk_id2", label="RUV k", value=1, step=1),
    #             selectizeInput(inputId="cntrs_id", label="Contrasts", 
    #                 choices=NULL)
    #         ),
    #         mainPanel(
    #             tabsetPanel(
    #                 id = 'Contrasts',
    #                 tabPanel(title="S3 - WT", 
    #                     DT::dataTableOutput("tab_id1")
    #                 )
    #             )
    #         )
    #     )
    # )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    observe({
        updateSelectizeInput(session=session, inputId="gene_id",
            choices = rownames(scelist[[1]]), 
            selected="Arc", server = TRUE)
        
        updateSelectizeInput(session=session, inputId="ct_id",
                             choices = unique(scelist[[1]]$labels), 
                             selected="CA1-do",
                             server = TRUE)
        
        updateSelectizeInput(session=session, inputId="gene_id1",
                             choices = rownames(mnnscelist[[1]]), 
                             selected="Arc", server = TRUE)
        
        updateSelectizeInput(session=session, inputId="ct_id1",
                             choices = unique(mnnscelist[[1]]$labels), 
                             selected="CA1-do",
                             server = TRUE)
        # updateSelectizeInput(session=session, inputId="ct_id",
        #     choices = unique(pseudo$label), selected="CA1-do",
        #     server = TRUE)
        # updateSelectizeInput(session=session, inputId="cntrs_id",
        #     choices = c("All", 
        #                 "S3 - WT", 
        #                 "S3SD - WT",
        #                 "SD - WT",
        #                 "S3 - S3SD",
        #                 "S3 - SD",
        #                 "SD - S3SD",
        #                 "Custom"
        #     ), 
        #     selected="All",
        #     server = TRUE)
        # updateSelectizeInput(session=session, inputId="ct_id2",
        #                      choices = unique(pseudo$label), selected="CA1-do",
        #                      server = TRUE)
    })
    # observeEvent(input$mnn_id, {
    #     if(input$mnn_id){
    #         # shinyjs::hide(id = "sample_id")
    #         updateSelectizeInput(session=session, inputId="sample_id",
    #                              choices = names(mnnscelist), 
    #                              selected=names(mnnscelist)[1], server = TRUE)
    #         ctcho <- unique(mnnscelist[[input$sample_id]]$SingleR)
    #         gncho <- rownames(mnnscelist[[input$sample_id]])
    #     }else{
    #         # shinyjs::show(id = "sample_id")
    #         ctcho <- unique(scelist[[input$sample_id]]$SingleR)
    #         gncho <- rownames(scelist[[input$sample_id]])
    #     }
    # })
    
    
    
    output$plot_id_out <- renderPlot({
        # if(!input$mnn_id)
        # {
            idx <- which(names(scelist) == input$sample_id)
            sce <- scelist[[idx]]
            plotMultiomeSample(sce, input$gene_id, input$sample_id, input$ct_id)
        # } else {
            
        # }
        
    })
    
    output$plot_id_out1 <- renderPlot({
        # if(!input$mnn_id)
        # {
        idx <- which(names(mnnscelist) == input$cond_id)
        mnnsce <- mnnscelist[[idx]]
        plotMultiomeSample(mnnsce, input$gene_id1, input$cond_id, input$ct_id1)
        # } else {
        
        # }
        
    })
    
    
    # output$plot_pseudo_out <- renderPlot({
    #      pseudoRUVPlot(pseudo, input$ct_id, input$ruvk_id, input$zoom_k_id)
    # })
    # output$ui_contrasts <- renderUI({
    #     
    #     tabsetPanel(
    #         id = 'Contrasts',
    #         for(i in seq_along(degsgglist))
    #         tabPanel(names(degsgglist)[i], 
    #             DT::dataTableOutput(paste0("tab_id", i))#names(degsgglist)[i])
    #         )
    #         # tabPanel("mtcars", DT::dataTableOutput("mytable2")),
    #         # tabPanel("iris", DT::dataTableOutput("mytable3"))
    #         
    #     )
    # })
        # if(!is.null(input$cntrs_id))
        # {
        #     cntrs <- processInputContrasts(input$cntrs_id)
        #     degsgglist <<- useEdgeR(pseudo, input$ct_id2, input$ruvk_id2, cntrs,
        #         input$ruvsg_id)
        # }
    # degsgglist <- plotedgeres(edgeRres)
    # output$tab_id1 <- DT::renderDataTable({
    #     DT::datatable(degsgglist[[1]]$DEGs)
    #     })
    # output$tab_id1 <- DT::renderDataTable(degsgglist[[1]]$DEGs)
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
