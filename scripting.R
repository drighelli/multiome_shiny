scelist <- readRDS("~/Downloads/multiome_data/scelist_azi_sr_altexp_filtered_upd.RDS")
BiocManager::install("scran")
library(scran)
scelist1 <- lapply(scelist, function(sce)
{
    altExp(sce) <- NULL
    sce <- logNormCounts(sce)
    hvg.sce.var <- getTopHVGs(sce, n=1000)
    sce <- runPCA(sce, subset_row=hvg.sce.var)
    sce <- runTSNE(sce, dimred="PCA")
    sce
})
saveRDS(scelist1, file="~/Downloads/multiome_data/scelist_shiny_14_June.RDS")
scelist <- readRDS("~/Downloads/multiome_data/scelist_shiny_14_June.RDS")
names(scelist)
scelistwt <- scelist[grep("_SD|_WT",names(scelist))]
saveRDS(scelistwt, file="~/Downloads/multiome_data/scelistwt_shiny_14_June.RDS")


# library(scater)
# plotTSNE(scelist1[[1]])
