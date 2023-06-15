
# gene <- "Arc"

library(ggpubr)
library(scater)
# lapply(seq_along(scelist), function(i)
# {
    # sce <- scelist[[i]]
plotMultiomeSample <- function(sce, gene, sampleid, ctid)
{
    
    # colname <- "Azimuth"
    # gga <- plotTSNE(sce, colour_by=colname) + ggtitle(paste0(sampleid, " ", colname))
    colname <- "SingleR"
    ggs <- plotTSNE(sce, colour_by=colname) + ggtitle(paste0(sampleid, " ", colname))
    gexg <- as.vector(logcounts(sce[which(rownames(sce) == gene),]))
    colData(sce)[paste0(`gene`,"_col")] <- gex
    sce  <- sce[,order(colData(sce)[paste0(`gene`,"_col")])]
    ggg1 <- plotTSNE(sce, colour_by=gene) + ggtitle(paste0(sampleid, " ", gene))
    
    # gex <- as.vector(logcounts(sce[which(rownames(sce) == gene),]))
    idxhigh <- which(gexg>0)
    idx0 <- which(gexg==0)
    idxlow <- which(gexg<0)
    gex <- gexg
    
    if ( length(idxhigh) > 0 ) {gex[idxhigh] <- "> 0"}
    if ( length(idx0) > 0 ) {gex[idx0] <- "= 0"}
    if ( length(idxlow) > 0 ) {gex[idxlow] <- "< 0"}
    gex <- as.factor(gex)
    colors <- NULL
    alpha <- NULL
    for (c in levels(gex)) {
        if(c=="= 0") {colors = c(colors, "antiquewhite"); alpha=c(alpha,"0.1")} 
        if(c=="< 0") {colors = c(colors, "red"); alpha=c(alpha,"1")} 
        if(c=="> 0") {colors = c(colors, "green"); alpha=c(alpha,"1")} 
    }
    colData(sce)[paste0(`gene`,"_col")] <- gex
    sce  <- sce[,order(colData(sce)[paste0(`gene`,"_col")])]
    ggg2 <- plotTSNE(sce, colour_by=paste0(gene,"_col")) +
        ggtitle(paste0(sampleid, " low/0/high ", gene)) +
        scale_colour_manual(values=colors) +
        scale_alpha_manual(values=alpha)
    
    ggg3 <- plotTSNE(sce[,sce$SingleR==ctid], colour_by=gene) + 
        ggtitle(paste(sampleid, ctid, gene))
    
    gg <- ggarrange(ggs, ggg3, ggg1, ggg2)
    return(gg)
    # ggsave(filename=paste0("plots/tSNE_", names(scelist)[i], ".pdf"), plot=gg, device="pdf")
}
# })


 
# gene <- "Arc"
# 
# library(ggpubr)
# lapply(seq_along(scelist), function(i)
# {
#     sce <- scelist[[i]]
#     colname <- "Azimuth"
#     gga <- plotTSNE(sce, colour_by=colname) + ggtitle(paste0(names(scelist)[[i]], " ", colname))
#     colname <- "SingleR"
#     ggs <- plotTSNE(sce, colour_by=colname) + ggtitle(paste0(names(scelist)[[i]], " ", colname))
# 
#     ggg <- plotTSNE(sce, colour_by=gene) + ggtitle(paste0(names(scelist)[[i]], " ", gene))
# 
#     gex <- as.vector(logcounts(sce[which(rownames(sce) == gene),]))
# 
#     gex[which(gex!=0)] <- "> 0"
#     gex[which(gex==0)] <- "= 0"
#     gex <- as.factor(gex)
#     colData(sce)[paste0(`gene`,"_col")] <- gex
#     sce  <- sce[,order(colData(sce)[paste0(`gene`,"_col")])]
#     ggg2 <-
# 
#         plotTSNE(sce, colour_by=paste0(gene,"_col")) +
#         ggtitle(paste0(names(scelist)[[i]], " 0/1 ", gene)) +
#         scale_colour_manual(values=c("antiquewhite", "red")) +
#         scale_alpha_manual(values=c("0.1", "1"))
# 
#     gg <- ggarrange(gga, ggs, ggg, ggg2)
# 
#     ggsave(filename=paste0("plots/tSNE_", names(scelist)[i], ".pdf"), plot=gg, device="pdf")
# })
