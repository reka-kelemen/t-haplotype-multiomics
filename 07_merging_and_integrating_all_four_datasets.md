# Merging and integrating all four datasets

We merged all four datasets in Seurat and clustered the nuclei:

    mc <- merge(m1,list(m2,m3,m4,murat),add.cell.ids=1:5,collapse=F,merge.dr=F,merge.data=F)
    mc <- SCTransform(mc,vars.to.regress = c("percent.mt","percent.ribo"),variable.features.n = 10000)
    mc[["xchr.autosome.ratio"]] <- colMeans(mc[["SCT"]]$data[rownames(mc)%in%XChrGenes,])/
    colMeans(mc[["SCT"]]$data[rownames(mc)%in%AutosomalGenes,])
    mc <- RunPCA(mc)
    mc <- FindNeighbors(mc, dims = 1:30)
    mc <- FindClusters(mc, resolution = 2)
    mc <- RunUMAP(mc, dims = 1:30, reduction = "pca")
    mci <- IntegrateLayers(object = mc, method = CCAIntegration, orig.reduction = "pca",
     new.reduction = "integrated.cca",normalization.method = "SCT",verbose = FALSE)
    mci <- FindNeighbors(mci, dims = 1:20,reduction="integrated.cca")
    mci <- FindClusters(mci, resolution = 1.2,reduction="integrated.cca")
    mci <- RunUMAP(mci, dims = 1:20, reduction="integrated.cca",
    reduction.name="umap.integrated",min.dist=0.1)
