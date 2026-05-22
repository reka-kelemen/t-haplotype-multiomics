# Removing doublets

To remove barcodes that likely represent two nuclei, we used the DoubletFinder package in R. First we performed dimensionality reduction and clustered nuclei with Seurat in R:

    m1 <- FindVariableFeatures(m1, selection.method = "vst", nfeatures = 10000)
    ## get rid of t complex and sex chr genes
    VariableFeatures(m1) <- VariableFeatures(m1)[(VariableFeatures(m1)%in%AutosomalGenes)] 
    m1 <- ScaleData(m1)
    m1 <- RunPCA(m1, features = VariableFeatures(object = m1))
    m1 <- FindNeighbors(m1, dims = 1:20)
    m1 <- FindClusters(m1, resolution = 0.5)
    m1 <- RunUMAP(m1, dims = 1:20)

    ### find doublets, and subset object to remove them ####
    sweep.list <- paramSweep(m1,PCs = 1:20,sct = F)
    sweep.stats <- summarizeSweep(sweep.list, GT = F)
    my.pK <- find.pK(sweep.stats)
    pKMax <- as.numeric(paste(my.pK$pK[which(my.pK$BCmetric==max(my.pK$BCmetric))])) 
    nExp_poi <- round(0.075*nrow(m1@meta.data))
    m1 <- doubletFinder(m1,PCs = 1:20,pN = 0.25,pK = pKMax,nExp = nExp_poi,reuse.pANN = FALSE, sct = FALSE)
    m1 <- m1[,m1@meta.data[,ncol(m1@meta.data)]=="Singlet"]
