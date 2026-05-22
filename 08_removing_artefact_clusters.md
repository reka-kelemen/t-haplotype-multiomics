# Removing artefact clusters

We assessed the following aspects of cell clusters in order to find potential artefact clusters:

- Clear expression of the markers of one specific cell type

- UMAP location and cell density of cluster in the Kelemen and Murat datasets

- Low mean X chromosome to mean autosome expression ratios if the cluster shows meiotic marker gene expression

- Bimodal distribution of mean X chromosome to mean autosome expression ratios if the cluster expressed haploid spermatid marker genes

- Fraction of intronic reads compared to other clusters of similar cell types

- Fraction of nuclei present in cluster compared to the fraction observed in the Murat dataset

We removed clusters that showed deviations in all of these categories and integrated the filtered dataset again. Then we searched for artefact clusters again using the above criteria, removed them, and inegrated all four of our datasets without the Murat dataset for further analyses.

    mc <- merge(m1,list(m2,m3,m4),add.cell.ids=1:4,collapse=F,merge.dr=F,merge.data=F)
    mc <- SCTransform(mc,vars.to.regress = c("percent.mt","percent.ribo","intronic.read.fraction"),
    variable.features.n = nVarFeatures)
    mc <- RunPCA(mc,verbose=F)
    mci <- IntegrateLayers(object = mc, method = CCAIntegration, orig.reduction = "pca",
     new.reduction = "integrated.cca",
                          normalization.method = "SCT",verbose = FALSE)
