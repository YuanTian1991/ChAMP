if(getRversion() >= "3.1.0") utils::globalVariables(c("myLoad","probe. features.epic","probe.features"))

CpG.GUI <- function(CpG=rownames(myLoad$beta),
                    arraytype="450K")
{
    if(arraytype=="EPIC") data(probe.features.epic) else data(probe.features)

    cgi.info <- table(probe.features[CpG,"cgi"])
    chromsome.info <- table(probe.features[CpG,"CHR"])
    feature.info <- table(probe.features[CpG,"feature"])
    type.info <- table(probe.features[CpG,"Type"])

    app <- shinyApp(
        ui = fluidPage(
                tags$head(tags$style("#cgi{height:40vh !important;}")),
                tags$head(tags$style("#chromsome{height:40vh !important;}")),
                tags$head(tags$style("#type{height:40vh !important;}")),
                tags$head(tags$style("#feature{height:40vh !important;}")),
                theme = shinytheme("readable"),
                titlePanel("CpG Overview"),
                tabsetPanel(
                            tabPanel("CpG",
                                     fluidRow(
                                              #column(width = 3,
                                              #       br(),
                                              #       sidebarPanel(  width=12,
                                              #                      numericInput('n', 'Number of obs', 100)
                                              #                   )
                                              #      ),
                                              column(width = 12,
                                                     fluidRow(
                                                              column(width = 6,
                                                                     plotlyOutput("chromsome"),
                                                                     plotlyOutput("feature")
                                                                    ),
                                                              column(width = 6,
                                                                     plotlyOutput("cgi"),
                                                                     plotlyOutput("type")
                                                                    )
                                                             )
                                                    )#column 
                                              )#fluidRow 
                                     )#tabPanel
                            )#tabsetPanel
                    ),#ui
    server = function(input, output){
        output$cgi <- renderPlotly({
                                     x <- list(
                                                title = ""
                                              )
                                      p <- plot_ly(x = names(cgi.info),
                                                   y = cgi.info,
                                                   name = "SF Zoo",
                                                   color = "orange",
                                                   type = "bar") %>%
                                           layout(xaxis = x)
                                  })
        output$chromsome <- renderPlotly({
                                     x <- list(
                                                title = ""
                                              )
                                      p <- plot_ly(x = paste("chr",names(chromsome.info),sep=""),
                                                   y = chromsome.info,
                                                   name = "SF Zoo",
                                                   type = "bar") %>%
                                           layout(xaxis = x)
                                  })
        output$feature <- renderPlotly({
                                     x <- list(
                                                title = ""
                                              )
                                      p <- plot_ly(x = names(feature.info),
                                                   y = feature.info,
                                                   name = "SF Zoo",
                                                   type = "bar") %>%
                                           layout(xaxis = x)
                                  })
        output$type <- renderPlotly({
                                     x <- list(
                                                title = ""
                                              )
                                      p <- plot_ly(x = names(type.info),
                                                   y = type.info,
                                                   name = "SF Zoo",
                                                   color = "orange",
                                                   type = "bar") %>%
                                           layout(xaxis = x)
                                  })
        }
    )
    runApp(app)
}
