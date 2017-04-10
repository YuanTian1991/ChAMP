if(getRversion() >= "3.1.0") utils::globalVariables(c("myLoad","probe.features.epic","probe.features","cmdscale","x","y","Sample_Name","economics","unemploy"))

QC.GUI <- function(beta=myLoad$beta,
                   pheno=myLoad$pd$Sample_Group,
                   arraytype="450K")
{
    innermdsplot <- function(beta,pheno)
    {
        o <- order(-matrixStats::rowVars(beta))[1:1000]
        d <- dist(t(beta[o, ]))
        fit <- cmdscale(d)
        col <- brewer.pal(8, "Dark2")[as.factor(pheno)]
        color_number <- nlevels(pheno)
        if(color_number < 3) color_number <- 3
        pal <- RColorBrewer::brewer.pal(color_number, "Set1")
        data <- data.frame(x=fit[,1],y=fit[,2],pheno=pheno,Sample_Name=colnames(beta))
        p <- plot_ly(data = data, x = ~x, y = ~y, text=~Sample_Name, color = ~pheno,colors = ~pal, type="scatter", mode = "markers", marker=list(size=15))
        m = list(l = 100,r = 50,b = 50,t = 100,pad = 10)
        p <- layout(p, title = 'MDS 1000 most variable positions',margin=m)
    }
    innerdensityPlot <- function(beta,pheno)
    {
        dens <- apply(beta, 2, density, na.rm = TRUE)
        df <- data.frame(x = unlist(lapply(dens, "[[", "x")),
                         y = unlist(lapply(dens, "[[", "y")),
                         cut = rep(names(dens), each = length(dens[[1]]$x))
                        )

        newdf <- do.call(rbind,lapply(dens,function(x) c(x$x,x$y)))

        color_number <- length(unique(pheno))
        if(color_number < 3) color_number <- 3
        mycol <- brewer.pal(color_number, "Set2")
        names(mycol) <- unique(pheno)
        newdf <- data.frame(newdf[,seq(1,1024,by=2)],names(dens),mycol[pheno])
        ###
        p <- plot_ly()
        for(i in 1:nrow(newdf))
        {
            p <- add_trace(p,
                           x = as.numeric(newdf[i,1:256]),  # x0, x1
                           y = as.numeric(newdf[i,257:512]),  # y0, y1
                           mode = "lines",
                           line = list(shape = "spline",color = newdf[i,514], width = 2),
                           showlegend = F,
                           hoverinfo = "text",
                           # Create custom hover text
                           text = newdf[i,513],
                           type="scatter")
        }

        ###
       m = list(l = 100,r = 50,b = 50,t = 100,pad = 10)
       p <- layout(p, title = 'Density Plot for each Sample',margin=m)
    }
    innertypePlot <- function(beta,arraytype)
    {
       if(arraytype=="EPIC") data(probe.features.epic) else data(probe.features)
       d1 <- density(beta[which(probe.features[rownames(beta),"Type"]=="I"),])
       d2 <- density(beta[which(probe.features[rownames(beta),"Type"]=="II"),])
       twolines <- data.frame(d1x = d1$x, d1y = d1$y, d2x = d2$x, d2y = d2$y)

       p <- plot_ly(data = twolines, x = ~d1x, y = ~d1y, mode = "lines+markers",type="scatter", name="Type-I probes",marker=list(size=1),line=list(width=4)) %>% 
       add_trace( x= ~d2x , y = ~d2y, mode="lines+markers", type="scatter",name="Type-II probes",marker=list(size=1),line=list(width=4))

       #p <- add_trace(p=p,x= ~d2$x,y= ~d2$y,mode="lines",name="Type-II probes")
       m = list(l = 100,r = 50,b = 50,t = 100,pad = 10)
       p <- layout(p, title = 'Type-I and Type-II density plot',margin=m)
    }
    innerdendrogram <- function(beta,pheno)
    {
        if(ncol(beta) <= 10)
        {
            hc <- hclust(dist(t(beta)))
        }else
        {
            SVD <- svd(beta)
            rmt.o <- EstDimRMT(beta - rowMeans(beta))
            k <- rmt.o$dim
            if(k < 2) k <- 2
            M <- SVD$v[,1:k]
            rownames(M) <- colnames(beta)
            colnames(M) <- paste("Component",c(1:k))
            hc <- hclust(dist(M))
        }
        dend <- as.dendrogram(hc)
        MyColor <- rainbow(length(table(pheno)))
        names(MyColor) <- names(table(pheno))
        labels_colors(dend) <- MyColor[pheno[order.dendrogram(dend)]]
        dend <- dend %>% set("labels_cex",0.8)
        dend <- dend %>% set("leaves_pch", 19) %>% set("leaves_cex",0.6) %>% set("leaves_col",MyColor[pheno[order.dendrogram(dend)]])

        plot(dend,center=TRUE,main=paste("Dendrogram for ",nrow(beta), " probes",sep=""))
        legend("topright",fill=MyColor,legend=names(MyColor)) 
    }
    innerheatmap <- function(beta)
    {
        if(nrow(beta) < 1000) ncpg <- nrow(beta) else ncpg <- 1000
        o <- order(-matrixStats::rowVars(beta))[1:ncpg]
        data <- beta[o,]
        a <- hclust(dist(data))$order
        b <- hclust(dist(t(data)))$order
        data <- data[a,b]
        m = list(l = 200,r = 50,b = 50,t = 100,pad = 10)
        p <- plot_ly(z = data,x = colnames(data), y = rownames(data),type = "heatmap")
        p <- layout(p, title = paste("Heatmap for top ",ncpg," variable CpGs",sep=""),margin=m)
    }
    innertest <- function()
    {
        m <- loess(unemploy / pop ~ as.numeric(date), data = economics)
        p <- economics %>%
        plot_ly(x = date, y = unemploy / pop) %>%
        add_trace(y = fitted(m)) %>%
        layout(showlegend = F)
    }

    app <- shinyApp(
        ui = fluidPage(
                tags$head(tags$style("#mdsPlot{height:80vh !important;}")),
                tags$head(tags$style("#densityPlot{height:80vh !important;}")),
                tags$head(tags$style("#dendrograme{height:80vh !important;}")),
                tags$head(tags$style("#typedensityPlot{height:80vh !important;}")),
                tags$head(tags$style("#heatmap{height:80vh !important;}")),

                theme = shinytheme("readable"),
                titlePanel("QC Overview"),
                tabsetPanel(
                            tabPanel("mdsPlot",
                                     align = "center",
                                     plotlyOutput("mdsPlot")
                                     ),
                            tabPanel("TypeDensity",
                                     align = "center",
                                     plotlyOutput("typedensityPlot")
                                     ),
                            tabPanel("QCPlot",
                                     align = "center",
                                     plotlyOutput("densityPlot")
                                     ),
                            tabPanel("Dendrogram",
                                     align = "center",
                                     plotOutput("dendrograme")
                                     ),
                            tabPanel("Heatmap",
                                     align = "center",
                                     plotlyOutput("heatmap")
                                     )
                            )#tabsetPanel
                    ),#ui
    server = function(input, output){
        output$mdsPlot <- renderPlotly(
                                   {
                                       innermdsplot(beta,pheno)
                                   }
                                  )
        output$typedensityPlot <- renderPlotly(
                                   {
                                       innertypePlot(beta,arraytype)
                                   }
                                  )
        output$densityPlot <- renderPlotly(
                                   {
                                       innerdensityPlot(beta,pheno)
                                   }
                                  )
        output$dendrograme <- renderPlot(
                                   {
                                       innerdendrogram(beta,pheno) 
                                   }
                                  )
        output$heatmap <- renderPlotly(
                                   {
                                       innerheatmap(beta)
                                   }
                                  )
        }
    )
    runApp(app)
}
