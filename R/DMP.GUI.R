if(getRversion() >= "3.1.0") utils::globalVariables(c("myLoad","myDMP","myNorm","probe","proportion","cgi","gene","Number","Variance","Value","ID","Sample","diamonds","logFC","negtive_log10P.adj","info","size","Beta_Value","Pheno","economics","unemploy","pop","MatchGeneName"))

DMP.GUI <- function(DMP=myDMP[[1]],
                    beta=myNorm,
                    pheno=myLoad$pd$Sample_Group,
                    cutgroupnumber=4)
{
    data(MatchGeneName)
    innerheatmap <- function(data)
    {
        if(nrow(data) < 5000) ncpg <- nrow(data) else ncpg <- 5000
        o <- order(-matrixStats::rowVars(data))[1:ncpg]
        data <- data[o,]

        a <- hclust(dist(data))$order
        b <- hclust(dist(t(data)))$order
        data <- data[a,b]
        m = list(l = 200,r = 50,b = 50,t = 100,pad = 10)
        SampleColfunc <- colorRampPalette(c("brown1","deepskyblue"))
        p <- plot_ly(z = data,x = colnames(data), y = rownames(data),type = "heatmap")
        p <- layout(p, title = paste("Heatmap for ",nrow(data)," 0.05 significant CpGs",sep=""),margin=m)
    }
    innercgibarplot <- function(h)
    {
        h <- apply(h,1,function(x) if(sum(x>0)!=0) prop.table(x) else x)
        h <- data.frame(cgi=rep(colnames(h),each=nrow(h)),probe=rep(rownames(h),ncol(h)),proportion=as.vector(h))
        p <- plot_ly(data=h,x = ~probe, y = ~proportion, type = "bar", color =~cgi)
        m = list(l = 70,r = 30,b = 50,t = 50,pad = 10)
        p <- layout(p, title = "cgi proportion barplot",margin=m)
    }
    innerfeaturebarplot <- function(h)
    {
        h <- apply(h,1,function(x) if(sum(x>0)!=0) prop.table(x) else x)
        h <- data.frame(cgi=rep(colnames(h),each=nrow(h)),probe=rep(rownames(h),ncol(h)),proportion=as.vector(h))
        p <- plot_ly(data=h,x = ~probe, y = ~proportion, type = "bar", color = ~cgi)
        m = list(l = 70,r = 30,b = 50,t = 50,pad = 10)
        p <- layout(p, title = "feature proportion barplot",margin=m)
    }
    innerfeaturecgibarplot <- function(h)
    {
        h <- apply(h,1,function(x) if(sum(x>0)!=0) prop.table(x) else x)
        h <- data.frame(cgi=rep(colnames(h),each=nrow(h)),probe=rep(rownames(h),ncol(h)),proportion=as.vector(h))
        p <- plot_ly(data=h,x = ~probe, y = ~proportion, type = "bar", color = ~cgi)
        m = list(l = 70,r = 30,b = 100,t = 50,pad = 10)
        p <- layout(p, title = "feature-cgi proportion barplot",margin=m)
    }
    innergeneenrich <- function(select)
    {
        rank <- sort(table(select[,"gene"]),decreasing=TRUE)
        if(length(rank) > 70) rank <- rank[1:70]
        select <- select[select$gene %in% names(rank),]
        select$logFC <- select$logFC>0
        h <- table(as.character(select$gene),select$logFC)
        h <- h[order(rowSums(h),decreasing=T),]
        colnames(h) <- c("Hypo","Hyper")
        h <- data.frame(Variance=rep(colnames(h),each=nrow(h)),gene=rep(rownames(h),ncol(h)),Number=as.vector(h))
        p <- plot_ly(data=h,x = ~gene, y = ~Number, type = "bar", color = ~Variance)
        m = list(l = 70,r = 30,b = 150,t = 50,pad = 10)
        p <- layout(p, title = paste("Top ",length(rank)," Gene Mostly Enriched by Significant CpGs",sep=""),margin=m,barmode = "stack")
    }

    innerGenePlot <- function(select,Group)
    {
        message("<< Generating dmrplot >>")
        #select <- select[order(select$MAPINFO),]
        G <- lapply(Group,function(x) t(x))
        G <- lapply(G,function(h) data.frame(Sample=rep(colnames(h),each=nrow(h)),ID=rep(rownames(h),ncol(h)),pos=rep(as.numeric(as.factor(select$MAPINFO)),ncol(h)),Value=as.vector(h)))
        for(i in 1:length(G)) G[[i]] <- data.frame(G[[i]],pheno=names(G)[i])
        X <- do.call(rbind,G)
        X <- cbind(X,showtext=paste("cpgID:",X$ID," Sample:",X$Sample,sep=""))

        p <- plot_ly()
        p <- add_trace(p,data=X, x=~pos, y=~Value,type="scatter", mode ="markers",text=~showtext,color=~pheno,marker=list(opacity = 0.6,size = 4))

        message("<< Dots Plotted >>")

        Fit <- lapply(G,function(h) data.frame(pos=unique(h$pos),ID=unique(h$ID),Mean=aggregate(h$Value,by=list(h$pos),mean)[,2]))
        for(i in 1:length(Fit)) Fit[[i]] <- data.frame(Fit[[i]],pheno=paste(names(Fit)[i],"Mean"))
        df <- data.frame(x = unlist(lapply(Fit, "[[", "pos")),
                         y = unlist(lapply(Fit, "[[", "Mean")),
                         ID = unlist(lapply(Fit, "[[", "ID")),
                         cut = unlist(lapply(Fit, "[[", "pheno")))
        p <- add_trace(p,data=df, x = ~x, y = ~y, text = ~ID, type="scatter", color = ~cut,mode = 'lines+markers',line = list(shape = "spline",width = 3,opacity = 0.3),marker=list(size = 1))

        message("<< Mean line Plotted >>")

        Fit2 <- lapply(G,function(h) data.frame(pos=unique(h$pos),ID=unique(h$ID),Fitted=unique(fitted(loess(h$Value ~ h$pos)))))
        for(i in 1:length(Fit2)) Fit2[[i]] <- data.frame(Fit2[[i]],pheno=paste(names(Fit2)[i],"Loess"))
        df2 <- data.frame(x = unlist(lapply(Fit2, "[[", "pos")),
                         y = unlist(lapply(Fit2, "[[", "Fitted")),
                         ID = unlist(lapply(Fit2, "[[", "ID")),
                         cut = unlist(lapply(Fit2, "[[", "pheno")))
        p <- add_trace(p,data=df2, x = ~x, y = ~y, text = ~ID, color = ~cut, type="scatter" ,mode = 'lines+markers',line =list(shape ="spline",width=3,opacity=0.3,dash = "dash"),marker=list(size = 1))

        message("<< Loess line Plotted >>")

        regioncol <- c("1stExon"="#00eeee","3'UTR"="#8b008b","5'UTR"="#00ee00","Body"="#ffd700","IGR"="#9e9e9e","TSS1500"="#ff3030","TSS200"="#00688b")
        featureregion <- list()
        for(i in names(regioncol))
        {
            index <- which(select$feature==i)
            if(length(index)==0) next

            oldw <- getOption("warn")
            options(warn = -1)

            featureregion[[i]]<- data.frame(x1 = index-0.5,
                                            x2 = index+0.5,
                                            y = -0.05,
                                            IDfeature = i,
                                            color = regioncol[i])

            options(warn = oldw)
        }
        cgicol <- c("island"="#ff83fa","opensea"="#ffdab9","shelf"="#bbffff","shore"="#98fb98")
        cgiregion <- list()
        for(i in names(cgicol))
        {
            index <- which(select$cgi==i)
            if(length(index)==0) next

            oldw <- getOption("warn")
            options(warn = -1)

            cgiregion[[i]]<- data.frame(x1 = index-0.5,
                                        x2 = index+0.5,
                                        y = -0.1,
                                        IDfeature = i,
                                        color = cgicol[i])

            options(warn = oldw)

        }
        #############################################
        df <- rbind(do.call(rbind,featureregion),do.call(rbind,cgiregion))
        colorrecord <- c()
        for(i in 1:nrow(df))
        {
            legend = F
            if(!df$color[i] %in% colorrecord)
            {
                colorrecord <- c(colorrecord,as.character(df$color[i]))
                legend = T
            }
            p <- add_trace(p,
                           x = c(df$x1[i], df$x2[i]),  # x0, x
                           y = c(df$y[i],df$y[i]),  # y0, y1
                           name = df$IDfeature[i],
                           mode = "lines",
                           line = list(color = df$color[i], width = 10 , opacity=0.5),
                           showlegend = legend,
                           hoverinfo = "text",
                           # Create custom hover text
                           text = df$IDfeature[i],
                           type="scatter")
        }
        #############################################

        message("<< Cgi Bar Plotted >>")


        p <- layout(p, title = select$gene[1])
    }

    innervolcanoplot <- function(select)
    {
        h <- diamonds[sample(nrow(diamonds),nrow(select),replace=T),c("carat","price","carat","clarity")]
        colnames(h) <- c("logFC","negtive_log10P.adj","size","info")
        h[[1]] <- select$logFC
        h[[2]] <- -(1)*log(select$adj.P.Val,10)
        h[[3]] <- (abs(h[[1]])-min(abs(h[[1]])))/(max(abs(h[[1]])-min(abs(h[[1]]))))+(h[[2]]-min(h[[2]]))/(max(h[[2]]-min(h[[2]])))
        h[[4]] <- paste(rownames(select),select$gene,select$feature,select$cgi)
        if(nrow(h) > 1500) h2 <- h[sample(1:nrow(h),1500),]
        h2 <- h2[order(h2$size),]
        SampleColfunc <- colorRampPalette(c("brown1","deepskyblue"))
        p <- plot_ly(h2, x=logFC, y=negtive_log10P.adj, text=info,mode = "markers",showlegend=FALSE,marker=list(opacity=0.5,size=3+(size*4),color=SampleColfunc(nrow(h2))))
        m = list(l = 70,r = 70,b = 50,t = 50,pad = 10)
        if(nrow(h) <= 1500)
        {
            p <- layout(p, title = paste("Volcano plot for",nrow(h2),"significant CpGs"),margin=m)
        }else
        {
            p <- layout(p, title = paste("Volcano plot for 1500 out of total",nrow(h)," significant CpGs (Random)"),margin=m)
        }
    }


    innercpgplot <- function(c,cpgname)
    {
        if(class(c$Pheno)=="numeric")
        {
            p <- plot_ly(c, x = ~Pheno, y = ~Beta_Value, text = ~Sample_Name, color = ~Beta_Value, mode = "markers",marker=list(opacity=0.5,size=5))
            m = list(l = 70,r = 70,b = 50,t = 50,pad = 10)
            p <- layout(p, title = paste("Scatter plot for",cpgname),margin=m)
        } else {
            p <- plot_ly(c, y = ~Beta_Value, color = ~Pheno, type = "box", boxpoints = "all", jitter = 0.3,pointpos = 0)
            m = list(l = 70,r = 70,b = 50,t = 50,pad = 10)
            p <- layout(p, title = paste("Boxplot for",cpgname),margin=m)
        }
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
                tags$head(tags$style("#test{height:40vh !important;}")),
                tags$head(tags$style("#heatmap{height:80vh !important;}")),
                tags$head(tags$style("#cgibarplot{height:40vh !important;}")),
                tags$head(tags$style("#featurebarplot{height:40vh !important;}")),
                tags$head(tags$style("#featurecgibarplot{height:40vh !important;}")),
                tags$head(tags$style("#geneenrich{height:40vh !important;}")),
                tags$head(tags$style("#geneplot{height:40vh !important;}")),
                tags$head(tags$style("#volcanoplot{height:40vh !important;}")),
                tags$head(tags$style("#cpgplot{height:40vh !important;}")),
                tags$head(tags$style("#inc{height:40vh !important;}")),

                theme = shinytheme("readable"),
                titlePanel("DMP Overview"),
                
		        sidebarLayout(
                    sidebarPanel(
                        br(),
                        width=3,
                        sliderInput("pvaluebin","P value Cutoff:",min = 0,max = 1,step = 0.025,value = 0.05),
                        sliderInput("logFCbin","abs(logFC) Cutoff:",min = 0,max = 1,step = 0.025,value = 0),
                        actionButton("go", "Submit"),

                        br(),
                        br(),
                        br(),
                        textInput("genebin", "Gene Symbol:", value = "NFIX"),
                        actionButton("genesearch", "Submit"),

                        br(),
                        br(),
                        br(),
                        textInput("cpgbin", "CpG ID:", value = "cg06822689"),
                        actionButton("cpgsearch", "Submit")
                                ),
                        mainPanel( 
                             tabsetPanel(
                                 tabPanel("DMPtable",
                                          align = "left",
                                          div(dataTableOutput("table"), style = "font-size:75%")
                                         ),
                                 tabPanel("Heatmap",
                                          align = "center",
                                          plotlyOutput("heatmap")
                                         ),
                                 tabPanel("Feature&Cgi",
                                          align = "center",
                                          br(),
                                          fluidRow(
                                              column(width = 6,
                                                  plotlyOutput("cgibarplot")
                                                     ),
                                              column(width = 6,
                                                  plotlyOutput("featurebarplot")
                                                     ),#column
                                              column(width = 12,
                                                  plotlyOutput("featurecgibarplot")
                                                     )#column
                                                  )#fluidRow
                                         ),
                                 tabPanel("Gene",
                                          align = "center",
                                          fluidRow(
                                              align = "center",
                                              br(),
                                              column(width = 12,
                                                  plotlyOutput("geneplot")
                                                     ),
                                              column(width = 12,
                                                  align = "center",
                                                  br(),
                                                  htmlOutput("inc")
                                                     )#column
                                                  )#fluidRow
                                         ),#tabPanel
                                 tabPanel("CpG",
                                          br(),
                                          align = "center",
                                          fluidRow(
                                           #   column(width = 6,
                                           #       plotlyOutput("volcanoplot")
                                           #          ),
                                              column(width = 12,
                                                  plotlyOutput("geneenrich")
                                                     ),
                                              column(width = 6,
                                                  plotlyOutput("cpgplot")
                                                     )
                                                  )#fluidRow
                                         )#tabPanel
                                    )#tabsetPanel
                                 )
                             )
                    ),#ui
    server = function(input, output){

        CpG_repalce <- eventReactive(input$cpgsearch,
                    {
                        pvalueCutoff <- as.numeric(input$pvaluebin)
                        abslogFCCutoff <- as.numeric(input$logFCbin)
                        cpgid <- input$cpgbin

                        ### Generate Data for CpGplot
                        mycpg <- data.frame(Beta_Value=beta[cpgid,],Pheno=pheno,Sample_Name=colnames(beta))
                        Parameter <- list(mycpg=mycpg,cpgid=cpgid)
                    })
        Gene_repalce <- eventReactive(input$genesearch,
                    {
                        pvalueCutoff <- as.numeric(input$pvaluebin)
                        abslogFCCutoff <- as.numeric(input$logFCbin)
                        genename <- input$genebin

                        ### Generate Data for Geneplot
                        mygeneselect=DMP[which(DMP$gene==genename & DMP$adj.P.Val <= pvalueCutoff & abs(DMP$logFC) >= abslogFCCutoff),]
                        mygeneselect <- mygeneselect[order(mygeneselect$MAPINFO),]
                        if(class(pheno)=="numeric")
                        {
                            cut_group <- as.character(cut(pheno,4))
                            mygroup <- split(as.data.frame(t(beta[rownames(mygeneselect),])),cut_group)
                        } else {
                            mygroup <- split(as.data.frame(t(beta[rownames(mygeneselect),])),pheno)
                        }
                        Parameter <- list(mygeneselect=mygeneselect,mygroup=mygroup,genename=genename)
                    })
        Cutoff_repalce <- eventReactive(input$go,
                    {
                        gc()
                        pvalueCutoff <- as.numeric(input$pvaluebin)
                        abslogFCCutoff <- as.numeric(input$logFCbin)

                        mydmp <- DMP[which(DMP$adj.P.Val <= pvalueCutoff & abs(DMP$logFC) >= abslogFCCutoff),]
                        message(dim(mydmp))
                        ### Generate Data for Heatmap
                        myheatmapdata=beta[rownames(mydmp),] 

                        ### Generate Data for cgi Barplot
                        h.cgi <- rbind(AllProbe=table(as.factor(DMP$cgi)),
                                   HyperProbe=table(as.factor(mydmp$cgi)[mydmp$logFC >= abslogFCCutoff]),
                                   HypoProbe=table(as.factor(mydmp$cgi)[mydmp$logFC <= (-1)*abslogFCCutoff]))

                        ### Generate Data for feature Barplot
                        h.feature <- rbind(AllProbe=table(as.factor(DMP$feature)),
                                   HyperProbe=table(as.factor(mydmp$feature)[mydmp$logFC >= abslogFCCutoff]),
                                   HypoProbe=table(as.factor(mydmp$feature)[mydmp$logFC <= (-1)*abslogFCCutoff]))

                        ### Generate Data for feature.cgi Barplot
                        h.feature.cgi <- rbind(AllProbe=table(as.factor(DMP$feat.cgi)),
                                               HyperProbe=table(as.factor(mydmp$feat.cgi)[mydmp$logFC >= abslogFCCutoff]),
                                               HypoProbe=table(as.factor(mydmp$feat.cgi)[mydmp$logFC <= (-1)*abslogFCCutoff]))

                        ### Generate Data for geneenrich Barplot
                        myallgeneselect=mydmp[mydmp$gene!="",c("gene","adj.P.Val","logFC")]



                        Parameter <- list(pvalueCutoff=pvalueCutoff,
                                       abslogFCCutoff=abslogFCCutoff,
                                       mydmp=mydmp,
                                       myheatmapdata=myheatmapdata,
                                       h.cgi=h.cgi,
                                       h.feature=h.feature,
                                       h.feature.cgi=h.feature.cgi,
                                       myallgeneselect=myallgeneselect)
                    })

        output$table <- renderDataTable({
                                         if(class(pheno)=="numeric") {
                                            datatable <- Cutoff_repalce()$mydmp[,c("CHR","MAPINFO","gene","feature","cgi",names(DMP)[1:6])]
                                         } else {
                                            datatable <- Cutoff_repalce()$mydmp[,c("CHR","MAPINFO","gene","feature","cgi",names(DMP)[1:9])]
                                         }
                                            datatable <- data.frame(ID=rownames(datatable),datatable)
                                        },options = list(pageLength=20))
        output$heatmap <- renderPlotly({   
                                           innerheatmap(Cutoff_repalce()$myheatmapdata)
                                       })
        output$cgibarplot <- renderPlotly({   
                                           innercgibarplot(Cutoff_repalce()$h.cgi)
                                       })
        output$featurebarplot <- renderPlotly({   
                                           innerfeaturebarplot(Cutoff_repalce()$h.feature)
                                       })
        output$featurecgibarplot <- renderPlotly({   
                                           innerfeaturecgibarplot(Cutoff_repalce()$h.feature.cgi)
                                       })
        output$geneenrich <- renderPlotly({   
                                           innergeneenrich(Cutoff_repalce()$myallgeneselect)
                                       })
        output$geneplot <- renderPlotly({   
                                           innerGenePlot(Gene_repalce()$mygeneselect,Gene_repalce()$mygroup)
                                       })
#        output$volcanoplot <- renderPlotly({   
#                                           innervolcanoplot(Cutoff_repalce()$mydmp)
#                                       })
        output$cpgplot <- renderPlotly({   
                                           innercpgplot(CpG_repalce()$mycpg,CpG_repalce()$cpgid)
                                       })
        output$inc <- renderUI({
                 webpage <- tags$iframe(src=paste("https://www.wikigenes.org/e/gene/e/",MatchGeneName[which(MatchGeneName[,1]==Gene_repalce()$genename),2],".html",sep=""),width="60%",height="100%")
                   })
        }
    )
    runApp(app)
}
