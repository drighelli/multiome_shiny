library(scran)
library(scater)
library(scDblFinder)
scelistdbl <- lapply(seq_along(scelist), function(i)
{
    sce <- scelist[[i]]
    topgs <- getTopHVGs(sce, prop=0.1)
    dbl.scores <- computeDoubletDensity(sce, subset.row=topgs)
    dbl.calls <- doubletThresholding(data.frame(score=dbl.scores),
                                     method="griffiths", returnType="call")
    print(paste0("------------- ", names(scelist)[i], " -------------"))
    print(dim(sce))
    print(summary(dbl.calls))
    # sce <- sce[,dbl.calls=="singlet"]
    colData(sce) <- cbind.DataFrame(colData(sce), dbl.calls, dbl.scores)
    sce
})
names(scelistdbl) <- names(scelist)


gglist <- lapply(seq_along(scelistdbl), function(i)
{
    sce <- scelistdbl[[i]]
    gdb1 <- plotTSNE(sce[,sce$labels=="CA1-do"], colour_by="dbl.calls") + ggtitle(names(scelistdbl)[i])
    gdb2 <- plotTSNE(sce[,sce$labels=="CA1-do"], colour_by="dbl.scores") + ggtitle(names(scelistdbl)[i])
    ggarrange(gdb1, gdb2)
})


# it is not possible to compute doubles on MNN sce