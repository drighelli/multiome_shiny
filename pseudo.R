pseudoRUVPlot <- function(pseudo, cellType, K, zoom)
{
    require(RUVSeq)
    require(darioscripts)
    
    pseudo <- logNormCounts(pseudo)
    pseudo <- runPCA(pseudo)
    
    gg1 <- plotPCA(pseudo, colour_by="label", shape_by="condition") + 
        ggtitle("Pseudo-Bulk PCA", subtitle="logcounts")
    
    pseudo.ct <- pseudo[,pseudo$label==cellType]
    pseudo.ct <- logNormCounts(pseudo.ct)
    pseudo.ct <- runPCA(pseudo.ct)
    gg2 <- plotPCA(pseudo.ct, colour_by="condition", shape_by="replicate") + 
        ggtitle(paste0("Pseudo-Bulk PCA ", cellType), subtitle="logcounts")
    
    nc <- readRDS("~/Downloads/multiome_data/SD_neg_ctrl.RDS")
    nc <- nc[nc[,2] %in% rownames(pseudo.ct),]
    # dim(nc)
    
    groups <- makeGroups(pseudo.ct$condition)
    
    ggl <- list()
    ggl <- .myappend(ggl, list(gg1, gg2))
    if(!zoom)
    {
        for (k in 1:K)
        {
            normExprData <- .RUVsnorm(counts(pseudo.ct), nc, groups, k)
            gg1 <- PlotPCAPlotlyFunction(log1p(normExprData),
                design.matrix=colData(pseudo.ct),
                shapeColname="replicate", colorColname="condition",
                prefix.plot=paste0("RUVs k: ",k), show.plot.flag=FALSE, 
                cowplot=TRUE)
            
            normExprData <- .RUVgnorm(counts(pseudo.ct), nc, k)
            gg2 <- PlotPCAPlotlyFunction(log1p(normExprData),
                design.matrix=colData(pseudo.ct),
                shapeColname="replicate", colorColname="condition",
                prefix.plot=paste0("RUVg k: ",k), show.plot.flag=FALSE,
                cowplot=TRUE)
            ggl <- .myappend(ggl, list(gg1, gg2))
        }
    }else{
        normExprData <- .RUVsnorm(counts(pseudo.ct), nc, groups, K)
        ggl[[3]] <- PlotPCAPlotlyFunction(log1p(normExprData),
            design.matrix=colData(pseudo.ct),
            shapeColname="replicate", colorColname="condition",
            prefix.plot=paste0("RUVs k: ", K), show.plot.flag=FALSE,
            cowplot=TRUE)
        normExprData <- .RUVgnorm(counts(pseudo.ct), nc, K)
        ggl[[4]] <- PlotPCAPlotlyFunction(log1p(normExprData),
            design.matrix=colData(pseudo.ct),
            shapeColname="replicate", colorColname="condition",
            prefix.plot=paste0("RUVg k: ", K), show.plot.flag=FALSE,
            cowplot=TRUE)
    }
    
    
    gg <- ggarrange(plotlist=ggl)#, ncol=2,nrow=(length(ggl)/2))
    return(gg)
}

.myappend <- function(x, values)
{
    idx <- length(x)
    if(length(values) == 1)
    {
        x[[idx+1]] <- values
    }
    
    for(i in seq_along(values))
    {
        # print(i)
        x[[idx+i]] <- values[[i]]
    }
    
    return(x)
}

.RUVsnorm <- function(counts, nc, groups, k)
{
    ruvedSExprData <- RUVs(as.matrix(round(counts)), 
                           cIdx=nc[,2], scIdx=groups, k=k)
    return(ruvedSExprData$normalizedCounts)
}

.RUVgnorm <- function(counts, nc, k)
{
    ruvedSExprData <- RUVg(as.matrix(round(counts)), 
                           cIdx=nc[,2], k=k)
    return(ruvedSExprData$normalizedCounts)
}
# .processNegativeControls <- function(file.path="~/Downloads/multiome_data/SD_Negative_Controls.txt")
# {
#     require(biomaRt)
#     neg.ctrls <- read.table(file.path)
#     colnames(neg.ctrls) <- c("ensembl_id")
#     
#     ensembl.biomaRt <- useEnsembl(biomart = "genes", 
#                                   dataset = "mmusculus_gene_ensembl")
#     
#     ensembl.row <- neg.ctrls$ensembl_id
#     ensembl.id <- getBM(attributes=c('external_gene_name','ensembl_gene_id'), 
#                         filters = 'ensembl_gene_id', values = ensembl.row, mart = ensembl.biomaRt)
#     
#     ensembl.id <- ensembl.id[order(ensembl.id$ensembl_gene_id),]
#     
#     neg.ctrls <- as.data.frame(neg.ctrls[order(neg.ctrls$ensembl_id),])
#     colnames(neg.ctrls) <- c("ensembl_id")
#     neg.ctrls$symbol <- ensembl.id$external_gene_name
#     saveRDS(neg.ctrls, file="~/Downloads/multiome_data/SD_neg_ctrl.RDS")
#     return(neg.ctrls)
# }
