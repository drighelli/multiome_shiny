processInputContrasts <- function(cntrs_input)
{
    if(cntrs_input == "All")
    {
        return(c("S3 - WT", "S3SD - WT", "SD - WT", "S3 - S3SD", "S3 - SD", 
                 "SD - S3SD"))
    } else if(cntrs_input == "Custom") {
        return(NULL)
    } else {
        return(cntrs_input)
    }
}

useEdgeR <- function(pseudo, cellType, ruvk, contrasts, norm)
{
    require(darioscripts)
    require(RUVSeq)
    pseudo.ct <- pseudo[,pseudo$label==cellType]
    groups <- makeGroups(pseudo.ct$condition)
    nc <- readRDS("~/Downloads/multiome_data/SD_neg_ctrl.RDS")
    nc <- nc[nc[,2] %in% rownames(pseudo.ct),]
    
    if(norm=="RUVs")
    {
        ruvedSExprData <- RUVs(as.matrix(round(counts(pseudo.ct))), 
                               cIdx=nc[,2], scIdx=groups, k=ruvk)
    } else {
        ruvedSExprData <- RUVg(as.matrix(round(counts(pseudo.ct))), 
                               cIdx=nc[,2], k=ruvk)
    }
    
    
    colData(pseudo.ct) <- cbind.DataFrame(colData(pseudo.ct), ruvedSExprData$W)
    
    edgeRres <- applyEdgeR(counts(pseudo.ct), design.matrix=colData(pseudo.ct), 
       factors.column="condition", weight.columns=colnames(ruvedSExprData$W), 
       contrasts=contrasts, is.normalized=FALSE)
    pc <- readRDS(file="~/Downloads/multiome_data/SD_pos_ctrls.RDS")
    
    degggl <- list()
    for( i in seq_along(edgeRres) )
    {
        edgeRres[[i]]$gene <- rownames(edgeRres[[i]])
        edgeRres[[i]]$DEG <- FALSE
        edgeRres[[i]]$DEG[edgeRres[[i]]$FDR<0.05] <- TRUE
        print(sum(edgeRres[[i]]$DEG))
        # ggl[[i]] <- PlotVolcanoPlot(edgeRres[[i]], counts.dataframe=counts(pseudo.ct),
        #    design.matrix=colData(pseudo.ct), positive.ctrls.list=pc[,1],
        #    prefix.plot=names(edgeRres)[i], threshold=0.5, show.plot.flag=FALSE)
        gg <- luciaVolcanoPlot1(edgeRres[[i]], prefix=names(edgeRres)[i],
            positive.controls.df=pc, threshold=0.05)
        degggl[[i]] <- list(DEGs=edgeRres, ggp=gg)
    }
    if(length(edgeRres)!=1) names(degggl) <- names(edgeRres)
    # ggarrange(plotlist=ggl)
    # 
    return(degggl)
}

plotedgeres <- function(edgeRres)
{
    pc <- readRDS(file="~/Downloads/multiome_data/SD_pos_ctrls.RDS")
    
    degggl <- list()
    for( i in seq_along(edgeRres) )
    {
        edgeRres[[i]]$gene <- rownames(edgeRres[[i]])
        edgeRres[[i]]$DEG <- FALSE
        edgeRres[[i]]$DEG[edgeRres[[i]]$FDR<0.05] <- TRUE
        print(sum(edgeRres[[i]]$DEG))
        # ggl[[i]] <- PlotVolcanoPlot(edgeRres[[i]], counts.dataframe=counts(pseudo.ct),
        #    design.matrix=colData(pseudo.ct), positive.ctrls.list=pc[,1],
        #    prefix.plot=names(edgeRres)[i], threshold=0.5, show.plot.flag=FALSE)
        gg <- luciaVolcanoPlot1(edgeRres[[i]], prefix=names(edgeRres)[i],
                                positive.controls.df=pc, threshold=0.05)
        degggl[[i]] <- list(DEGs=edgeRres, ggp=gg)
    }
    if(length(edgeRres)!=1) names(degggl) <- names(edgeRres)
    # ggarrange(plotlist=ggl)
    # 
    return(degggl)
}