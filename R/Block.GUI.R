if(getRversion() >= "3.1.0") utils::globalVariables(c("myBlock","myLoad","myNorm","probe.features.epic","probe.features","Value","ID","Sample","aggregate","x","y","fitted","loess"))

Block.GUI <- function(Block=myBlock,
                      beta=myNorm,
                      pheno=myLoad$pd$Sample_Group,
                      runMVP=TRUE,
                      compare.group=NULL,
                      arraytype="450K")
{
    if(arraytype == "EPIC") data(probe.features.epic) else data(probe.features)
    probe.features <- probe.features[names(Block$allCLID.v),]

    design <- model.matrix( ~ 0 + pheno)
    fit <- lmFit(Block$avbetaCL.m, design)
    fit2 <- eBayes(fit)
    clusterMVP <- topTable(fit2,coef=1,number=nrow(Block$avbetaCL.m),adjust.method="BH")

    B2 <- makeGRangesFromDataFrame(Block$Block[,1:4])
    A2 <- makeGRangesFromDataFrame(data.frame(seqnames=Block$posCL.m[,2],start=Block$posCL.m[,1],end=Block$posCL.m[,1]))
    C <- as.data.frame(findOverlaps(B2,A2))
    clustertmp <- data.frame(Blockindex=paste("Block",C$queryHits,sep="_"),Block$posCL.m[C$subjectHits,2:1])
    clustertmp <- data.frame(clustertmp,clusterMVP[rownames(clustertmp),])

    message("Generation CpG information")
    cpgtmp <- probe.features[names(Block$allCLID.v[which(Block$allCLID.v %in% rownames(clustertmp))]),c(1,2,5,6,7)]
    if(runMVP)
    {
        message("Calculating MVP")
        cpgMVP <- champ.MVP(beta=beta,
                            pheno=pheno,
                            adjPVal=1,
                            adjust.method="BH",
                            compare.group=compare.group,
                            arraytype=arraytype)
        cpgtmp <- data.frame(cpgtmp,cpgMVP[rownames(cpgtmp),c(1,3,4,5)])
    }
    cpgtmp <- data.frame(Clusterindex=Block$allCLID.v[rownames(cpgtmp)],cpgtmp)
    cpgtmp <- data.frame(Blockindex=clustertmp[as.character(cpgtmp$Clusterindex),"Blockindex"],cpgtmp)
    cpgtmp <- cpgtmp[order(cpgtmp$Blockindex,cpgtmp$MAPINFO),]

    innerblockplot <- function(G,cpgidf,block.idx)
    {
        message("<< Generating dmrplot >>")

        for(i in 1:length(G)) G[[i]] <- data.frame(G[[i]],pheno=names(G)[i])
        X <- do.call(rbind,G)
        p <- plot_ly(data=X, x=pos, y=Value, mode = "markers", text = paste("cpgID:",ID," Sample:",Sample,sep=""),color=pheno,marker=list(opacity = 0.6,size = 4))

        message("<< Dots Plotted >>")

        Fit <- lapply(G,function(h) data.frame(pos=unique(h$pos),ID=unique(h$ID),Mean=aggregate(h$Value,by=list(h$pos),mean)[,2]))
        for(i in 1:length(Fit)) Fit[[i]] <- data.frame(Fit[[i]],pheno=paste(names(Fit)[i],"Mean"))
        df <- data.frame(x = unlist(lapply(Fit, "[[", "pos")),
                         y = unlist(lapply(Fit, "[[", "Mean")),
                         ID = unlist(lapply(Fit, "[[", "ID")),
                         cut = unlist(lapply(Fit, "[[", "pheno")))
        p <- add_trace(data=df, x = x, y = y, text = ID, color = cut,line = list(shape = "linear",width = 3,opacity = 0.3),marker=list(size = 1))

        message("<< Mean line Plotted >>")

        Fit2 <- lapply(G,function(h) data.frame(pos=unique(h$pos),ID=unique(h$ID),Fitted=unique(fitted(loess(h$Value ~ h$pos)))))
        for(i in 1:length(Fit2)) Fit2[[i]] <- data.frame(Fit2[[i]],pheno=paste(names(Fit2)[i],"Loess"))
        df2 <- data.frame(x = unlist(lapply(Fit2, "[[", "pos")),
                         y = unlist(lapply(Fit2, "[[", "Fitted")),
                         ID = unlist(lapply(Fit2, "[[", "ID")),
                         cut = unlist(lapply(Fit2, "[[", "pheno")))
        p <- add_trace(data=df2, x = x, y = y, text = ID, color = cut,line =list(shape ="spline",width=3,opacity=0.3,dash = "dash"),marker=list(size = 1))

        message("<< Loess line Plotted >>")

        if(!is.null(cpgidf))
            p <- add_trace(data=cpgidf,x=pos,y=y,text=ID,mode="markers",name="CpG Islands",marker=list(opacity = 0.6,size = 10,color="#8b3a62"))
        p <- layout(p, title = paste("Block",block.idx,sep="_"))
    }
    innertestplot <- function()
    {
        plot_ly(x=1:10,y=1:10)
    }

    app <- shinyApp(
        ui = fluidPage(
                tags$head(tags$style("#blockplot{height:40vh !important;}")),

                theme = shinytheme("readable"),
                titlePanel("Block Overview"),
                
		        sidebarLayout(
                    sidebarPanel(
                        br(),
                        width=3,
                        sliderInput("pvaluebin","P value Cutoff:",min = 0,max = 1,step = 0.025,value = 0.05),
                        sliderInput("mincluserbin","Min Number Clusters:",min = 1,max = 350,step = 1,value = 100),
                        actionButton("go", "Submit"),

                        br(),
                        br(),
                        br(),
                        sliderInput("blockbin","Block Index:",min = 1,max = nrow(Block$Block),step = 1,value = 1),
                        actionButton("blocksearch", "Submit")
                                ),#sidebarPanel

                        mainPanel( 
                             tabsetPanel(
                                 tabPanel("Blocktable",
                                          align = "left",
                                          div(dataTableOutput("blocktable"), style = "font-size:75%")
                                         ),
                                 tabPanel("CpGtable",
                                          align = "left",
                                          div(dataTableOutput("cpgtable"), style = "font-size:75%")
                                         ),
                                 tabPanel("Blockplot",
                                          align = "center",
                                          br(),
                                          fluidRow(
                                                   align = "center",
                                                   br(),
                                                   column(width = 12,
                                                          align = "left",
                                                          div(dataTableOutput("clustertable"), style = "font-size:75%")
                                                         ),#column
                                                   column(width = 12,
                                                          align = "center",
                                                          plotlyOutput("blockplot")
                                                         )#column
                                                  )#fluidRow
                                         )
                                    )#tabsetPanel
                                 )#mainPanel
                             )#sidebarLayout
                    ),#ui
    server = function(input, output){

        Block_repalce <- eventReactive(input$blocksearch,
                    {
                        gc()
                        pvalueCutoff <- as.numeric(input$pvaluebin)
                        mincluserCutoff <- as.numeric(input$mincluserbin)
                        block.idx <- as.numeric(input$blockbin)
                        
                        index <- which(Block$posCL.m[,"Chr"] == Block$Block[block.idx,"chr"] & 
                                       Block$posCL.m[,"Pos"] >= Block$Block[block.idx,"start"] & 
                                       Block$posCL.m[,"Pos"] <= Block$Block[block.idx,"end"])

                       clustersinblock <- clustertmp[which(clustertmp$Blockindex==paste("Block",block.idx,sep="_")),] 
                       cpg <- myBlock$allCLID.v[which(myBlock$allCLID.v %in% rownames(clustersinblock))]

                       cpg.pos <- aggregate(probe.features[names(cpg),"MAPINFO"],by=list(cpg),function(x) c(min(x),max(x),max(x)-min(x)))
                       cpg.pos.final <- cpg.pos[,2]
                       rownames(cpg.pos.final) <- cpg.pos[,1]
                       colnames(cpg.pos.final) <- c("Start","End","Width")

                       cpg.gene <- aggregate(probe.features[names(cpg),"gene"],by=list(cpg),function(x) paste(unique(x),collapse=","))
                       cpg.gene.final <- cpg.gene[,2]
                       names(cpg.gene.final) <- cpg.gene[,1]

                       clustersinblock <- data.frame(clustersinblock[,1:3],
                                                     cpg.pos.final[rownames(clustersinblock),],
                                                     Number=table(cpg)[rownames(clustersinblock)],
                                                     gene=cpg.gene.final[rownames(clustersinblock)],
                                                     clustersinblock[4:9])

                        OSindex <- index[which(substr(names(index),1,2)=="OS")]
                        CPGIindex <- index[which(substr(names(index),1,2)=="CP")]
                        SHindex <- index[which(substr(names(index),1,2)=="SH")]

                        Group <- split(as.data.frame(t(Block$avbetaCL.m[OSindex,])),pheno)
                        CPGIpos <- Block$posCL.m[CPGIindex,1]
                        G <- lapply(Group,function(x) t(x))
                        G <- lapply(G,function(h) data.frame(Sample=rep(colnames(h),each=nrow(h)),ID=rep(rownames(h),ncol(h)),pos=rep(Block$posCL.m[OSindex,1],ncol(h)),Value=as.vector(h)))

                        if(length(CPGIindex)>0)
                            cpgidf <- data.frame(ID=rownames(Block$posCL.m)[CPGIindex],pos=Block$posCL.m[CPGIindex,1],y=1)
                        else
                            cpgidf=NULL

                        list(pvalueCutoff=pvalueCutoff,
                             mincluserCutoff=mincluserCutoff,
                             G=G,
                             cpgidf=cpgidf,
                             block.idx=block.idx,
                             clustersinblock=clustersinblock)
                    })
        Cutoff_repalce <- eventReactive(input$go,
                    {
                        gc()
                        pvalueCutoff <- as.numeric(input$pvaluebin)
                        mincluserCutoff <- as.numeric(input$mincluserbin)
                        myblock <- Block$Block[which(Block$Block$p.valueArea <= pvalueCutoff & Block$Block$L >= mincluserCutoff),]
                        mycpgtable <- cpgtmp[which(cpgtmp$Blockindex %in% rownames(myblock)),]

                        list(pvalueCutoff=pvalueCutoff,
                             mincluserCutoff=mincluserCutoff,
                             myblock=myblock,
                             mycpgtable=mycpgtable)
                    })

        output$blocktable <- renderDataTable({
                                            datatable <- Cutoff_repalce()$myblock
                                            datatable <- data.frame(ID=rownames(datatable),datatable)
                                        },options = list(pageLength=20))
        output$cpgtable <- renderDataTable({
                                            datatable <- Cutoff_repalce()$mycpgtable
                                            datatable <- data.frame(ID=rownames(datatable),datatable)
                                        },options = list(pageLength=20))
        output$blockplot <- renderPlotly({
                                            innerblockplot(Block_repalce()$G,Block_repalce()$cpgidf,Block_repalce()$block.idx)
                                        })
        output$clustertable <- renderDataTable({
                                            datatable <- Block_repalce()$clustersinblock
                                            datatable <- data.frame(ID=rownames(datatable),datatable)
                                        },options = list(pageLength=8))

        }
    )
    runApp(app)
}
