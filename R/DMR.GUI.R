if(getRversion() >= "3.1.0") utils::globalVariables(c("myDMR","myLoad","myNorm","probe.features.epic","probe.features","Value","ID","Sample","aggregate","x","y","fitted","loess","gene","Number","Variance"))

DMR.GUI <- function(DMR=myDMR,
                    beta=myNorm,
                    pheno=myLoad$pd$Sample_Group,
                    runDMP=TRUE,
                    compare.group=NULL,
                    arraytype="450K")
{
    message("!!! important !!! Since we just upgrated champ.DMP() function, which is now can support multiple phenotypes. Here in DMR.GUI() function, if you want to use \"runDMP\" parameter, and your pheno contains more than two groups of phenotypes, you MUST specify compare.group parameter as compare.group=c(\"A\",\"B\") to get DMP value between group A and group B.")

    message("\n[ Section 1: Calculate DMP Start  ]\n")
    if(runDMP)
    {
        tmpbeta <- beta
        tmppheno <- pheno
        if(class(pheno)=="numeric") {
            message("  Your pheno parameter is numeric, champ.DMP() function would calculate linear regression for your CpGs.")
        } else {
            message("  You pheno is ",class(pheno)," type.")
            message("    Your pheno information contains following groups. >>")
            sapply(unique(pheno),function(x) message("    <",x,">:",sum(pheno==x)," samples."))
            
            if(length(unique(pheno))==2) {
                message("  Your pheno contains EXACTLY two phenotypes, which is good, compare.group is not needed.")
            } else {
                message("  Your pheno contains more than 2 phenotypes, please use compare.group to specify only two of them here.")
                if(is.null(compare.group)){
                    stop("  compare.group is needed here, please specify compare.group.")
                } else if (sum(compare.group %in% unique(pheno))==2) {
                    message("  Your compare.group is in accord with your pheno, which is good, now we are about to extract information for your compare.group.")
                    tmpbeta <- beta[,which(pheno %in% compare.group)]
                    tmppheno <- pheno[which(pheno %in% compare.group)]
                } else {
                    stop("  Seems your compare.group is not in accord with your pheno, please recheck your pheno and your compare.group.")
                }
            }
        }
        message("Calculating DMP")
        DMP <- champ.DMP(beta=tmpbeta,
                         pheno=tmppheno,
                         adjPVal=1,
                         adjust.method="BH",
                         compare.group=compare.group,
                         arraytype=arraytype)
        DMP <- DMP[[1]]
    }
    message("\n[ Section 1: Calculate DMP Done  ]\n")

    if(arraytype == "EPIC") data(probe.features.epic) else data(probe.features)
    probe.features <- probe.features[rownames(beta),]

    message("\n[ Section 2: Mapping DMR to annotation Start  ]\n")

    message("  Generating Annotation File")
    DMR[[1]]$seqnames <- as.factor(substr(DMR[[1]]$seqnames,4,100))
    index <- apply(DMR[[1]],1,function(x) which(probe.features$CHR==x[1] & probe.features$MAPINFO >= as.numeric(x[2]) & probe.features$MAPINFO <= as.numeric(x[3])))
    Anno <- data.frame(DMRindex=unname(unlist(sapply(names(index),function(x) rep(x,length(index[[x]]))))),probe.features[do.call(c,index),1:8])
    message("  Generating Annotation File Success")

    if(identical(names(DMR),"DMRcateDMR"))
    {
        sig <- data.frame(DMR.pvalue=rep(0,nrow(DMR$DMRcateDMR)),
                          DMR.probes=unlist(lapply(index,length)))
    }else if(identical(names(DMR),"BumphunterDMR"))
    {
        sig <- data.frame(DMR.pvalue=DMR$BumphunterDMR$p.valueArea,
                          DMR.probes=unlist(lapply(index,length)))
    }else if(identical(names(DMR),"ProbeLassoDMR"))
    {
        sig <- data.frame(DMR.pvalue=DMR$ProbeLassoDMR$dmrP,
                          DMR.probes=unlist(lapply(index,length)))
    }
    message("\n[ Section 2: Mapping DMR to annotation Done  ]\n")

    innerdmrplot <- function(select,Group,dmr.idx)
    {
        message("<< Generating dmrplot >>")
        #select <- select[order(select$MAPINFO),]
        G <- lapply(Group,function(x) t(x))
        G <- lapply(G,function(h) data.frame(Sample=rep(colnames(h),each=nrow(h)),ID=rep(rownames(h),ncol(h)),pos=rep(as.numeric(as.factor(select$MAPINFO)),ncol(h)),Value=as.vector(h)))
        for(i in 1:length(G)) G[[i]] <- data.frame(G[[i]],pheno=names(G)[i])
        X <- do.call(rbind,G)
        X <- cbind(X,showtext=paste("cpgID:",X$ID," Sample:",X$Sample,sep=""))
        p <- plot_ly()
        p <- add_trace(p,data=X, x=~pos, y=~Value,type="scatter", mode = "markers", text =~showtext,color=~pheno,marker=list(opacity = 0.6,size = 4))

        message("<< Dots Plotted >>")

        Fit <- lapply(G,function(h) data.frame(pos=unique(h$pos),ID=unique(h$ID),Mean=aggregate(h$Value,by=list(h$pos),mean)[,2]))
        for(i in 1:length(Fit)) Fit[[i]] <- data.frame(Fit[[i]],pheno=paste(names(Fit)[i],"Mean"))
        df <- data.frame(x = unlist(lapply(Fit, "[[", "pos")),
                         y = unlist(lapply(Fit, "[[", "Mean")),
                         ID = unlist(lapply(Fit, "[[", "ID")),
                         cut = unlist(lapply(Fit, "[[", "pheno")))
        p <- add_trace(p,data=df, x =~x, y =~y, text =~ID,type="scatter",mode="lines+markers", color =~cut,line = list(shape = "spline",width = 3,opacity = 0.3),marker=list(size = 1))

        message("<< Mean line Plotted >>")

        Fit2 <- lapply(G,function(h) data.frame(pos=unique(h$pos),ID=unique(h$ID),Fitted=unique(fitted(loess(h$Value ~ h$pos)))))
        for(i in 1:length(Fit2)) Fit2[[i]] <- data.frame(Fit2[[i]],pheno=paste(names(Fit2)[i],"Loess"))
        df2 <- data.frame(x = unlist(lapply(Fit2, "[[", "pos")),
                         y = unlist(lapply(Fit2, "[[", "Fitted")),
                         ID = unlist(lapply(Fit2, "[[", "ID")),
                         cut = unlist(lapply(Fit2, "[[", "pheno")))
        p <- add_trace(p,data=df2, x =~x, y =~y, text =~ID,type="scatter",mode="lines+markers", color =~cut,line =list(shape ="spline",width=3,opacity=0.3,dash = "dash"),marker=list(size = 1))

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

        p <- layout(p, title = paste("DMR",dmr.idx,sep="_"))
    }
    innergeneenrichplot <- function(select)
    {
        print(dim(select))
        message("<< Calculating Gene Enrich Plot >>")
        select <- select[which(select$gene != ""),]
        if(c("adj.P.Val") %in% colnames(select))
        {
            select$t <- select$t>0
            select$adj.P.Val[which(select$adj.P.Val <= 0.05)] <- 1
            select$adj.P.Val[which(select$adj.P.Val != 1)] <- 0
            h <- table(as.character(select$gene),paste(select$t,select$adj.P.Val))
            h <- h[order(rowSums(h),decreasing=T),]
            h <- cbind(h[,c(4,2)],rowSums(h[,c(1,3)]))
            colnames(h) <- c("Hyper","Hypo","Not Significance")
            h <- data.frame(Variance=rep(colnames(h),each=nrow(h)),gene=rep(rownames(h),ncol(h)),Number=as.vector(h))
            p <- plot_ly(data=h,x =~gene, y =~Number, type = "bar", color =~Variance)
        }else
        {
            h <- data.frame(table(as.character(select$gene)))
            colnames(h) <- c("gene","Number")
            h <- h[order(h$Number,decreasing=T),]
            p <- plot_ly(data=~h,x=~gene,y=~Number,type="bar")
        }
        m = list(l = 70,r = 30,b = 150,t = 50,pad = 10)
        p <- layout(p, title = paste("All ",length(unique(h$gene))," Gene Enriched by DMR-related CpGs",sep=""),margin=m,barmode = "stack")
    }
    innerheatmap <- function(data)
    {
        m = list(l = 170,r = 70,b = 60,t = 50,pad = 10)
        p <- plot_ly(z = data,x = colnames(data), y = rownames(data),type = "heatmap")
        p <- layout(p, title = paste("Heatmap for ",nrow(data)," CpGs in all DMRs",sep=""),margin=m)
    }
    innertestplot <- function()
    {
        plot_ly(x=1:10,y=1:10)
    }

    app <- shinyApp(
        ui = fluidPage(

                tags$head(tags$style("#dmrplot{height:40vh !important;}")),
                tags$head(tags$style("#dmrcateplot{height:80vh !important;}")),
                tags$head(tags$style("#geneenrichplot{height:40vh !important;}")),
                tags$head(tags$style("#heatmap{height:40vh !important;}")),

                theme = shinytheme("readable"),
                titlePanel(paste(names(DMR)[1],"Overview")),
                
		        sidebarLayout(
                    sidebarPanel(
                        br(),
                        width=3,
                        sliderInput("pvaluebin","P value Cutoff:",min = 0,max = 1,step = 0.025,value = 0.05),
                        sliderInput("minprobebin","Min Number Probes:",min = 1,max = 100,step = 1,value = 7),
                        actionButton("go", "Submit"),

                        br(),
                        br(),
                        br(),
                        sliderInput("dmrbin","DMR Index:",min = 1,max = nrow(sig),step = 1,value = 1),
                        actionButton("dmrsearch", "Submit")
                                ),#sidebarPanel

                        mainPanel( 
                             tabsetPanel(
                                 tabPanel("DMRtable",
                                          align = "left",
                                          div(dataTableOutput("table"), style = "font-size:75%")
                                         ),
                                 tabPanel("CpGtable",
                                          align = "left",
                                          div(dataTableOutput("cpgtable"), style = "font-size:75%")
                                         ),
                                 tabPanel("DMRPlot",
                                          align = "center",
                                          fluidRow(
                                                   align = "center",
                                                   br(),
                                                   column(width = 12,
                                                          align = "left",
                                                          div(dataTableOutput("probetable"), style = "font-size:75%")
                                                         ),#column
                                                   column(width = 12,
                                                          align = "center",
                                                          plotlyOutput("dmrplot")
                                                         )#column
                                                  )#fluidRow
                                         ),
                                 tabPanel("Summary",
                                          align = "center",
                                          fluidRow(
                                                   align = "center",
                                                   br(),
                                                   column(width = 12,
                                                          align = "center",
                                                          plotlyOutput("geneenrichplot")
                                                         ),#column
                                                   column(width = 12,
                                                          align = "center",
                                                          plotlyOutput("heatmap")
                                                         )#column
                                                  )#fluidRow
                                         )
                                    )#tabsetPanel
                                 )#mainPanel
                             )#sidebarLayout
                    ),#ui
    server = function(input, output){

        DMR_repalce <- eventReactive(input$dmrsearch,
                    {
                        gc()
                        pvalueCutoff <- as.numeric(input$pvaluebin)
                        minprobeCutoff <- as.numeric(input$minprobebin)
                        dmr.idx <- as.numeric(input$dmrbin)

                        ### Generate Data for dmrplot
                        mydmrselect <- Anno[Anno$DMRindex==paste("DMR",dmr.idx,sep="_"),]
                        mydmrselect <- mydmrselect[order(mydmrselect$MAPINFO),]
                        mygroup <- split(as.data.frame(t(beta[rownames(mydmrselect),])),pheno)
                        list(mydmrselect=mydmrselect,
                             mygroup=mygroup,
                             dmr.idx=dmr.idx)
                    })
        Cutoff_repalce <- eventReactive(input$go,
                    {
                        gc()
                        pvalueCutoff <- as.numeric(input$pvaluebin)
                        minprobeCutoff <- as.numeric(input$minprobebin)
                        mydmr <- DMR[[1]][which(sig$DMR.pvalue <= pvalueCutoff & sig$DMR.probes >= minprobeCutoff),]

                        ## Generate data for geneenrichplot
                        mygeneselect <- Anno[which(Anno$DMRindex %in% rownames(mydmr)),]

                        if(runDMP) mygeneselect <- data.frame(mygeneselect,DMP[rownames(mygeneselect),c("t","adj.P.Val")])
                        
                        ## Generate data for heatmap
                        mydata <- beta[rownames(mygeneselect),]
                        rownames(mydata) <- paste(mygeneselect$DMRindex,rownames(mydata),sep="-")

                        list(pvalueCutoff=pvalueCutoff,
                             minprobeCutoff=minprobeCutoff,
                             mydmr=mydmr,
                             mygeneselect=mygeneselect,
                             mydata=mydata)
                    })

        output$probetable <- renderDataTable({
                                            datatable <- DMR_repalce()$mydmrselect
                                            datatable <- data.frame(ID=rownames(datatable),datatable)
                                        },options = list(pageLength=8))
        output$table <- renderDataTable({
                                            datatable <- Cutoff_repalce()$mydmr
                                            datatable <- data.frame(ID=rownames(datatable),datatable)
                                        },options = list(pageLength=20))
        output$cpgtable <- renderDataTable({
                                            datatable <- Cutoff_repalce()$mygeneselect
                                            datatable <- data.frame(ID=rownames(datatable),datatable)
                                           },options = list(pageLength=20))
        output$dmrplot <- renderPlotly({
                                            innerdmrplot(DMR_repalce()$mydmrselect,DMR_repalce()$mygroup,DMR_repalce()$dmr.idx)
                                       })
        output$geneenrichplot <- renderPlotly({
                                            innergeneenrichplot(Cutoff_repalce()$mygeneselect)
                                       })
        output$heatmap <- renderPlotly({
                                            innerheatmap(Cutoff_repalce()$mydata)
                                       })

        }
    )
    runApp(app)
}
