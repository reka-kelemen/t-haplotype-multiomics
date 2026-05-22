# Filtering for high-quality nuclei

We used the R package Seurat to filter for high-quality nuclei. First, to be able to compare to an external guide dataset when looking for artefact clusters, we filtered for genes present in the single nucleus RNA-seq dataset from ([Murat et al. 2023](references.md#ref-murat2023molecular)). To remove the effect of the *t*-haplotype in identifying cell types, we removed genes between 3-42 Mbs on chromosome 17. We filtered for nuclei that have their percent of mitochondrial UMIs below 5 percent, and their intronic read fractions above 40 percent.

    m1.data <-  Read10X_h5("223616_output_3_filtered.h5", use.names = TRUE)
    m1.data <- m1.data[rownames(m1.data)%in%muratgenes,]
    m1 <- CreateSeuratObject(counts = m1.data, project = paste(id))
    m1[["percent.mt"]] <- PercentageFeatureSet(m1, pattern = "^mt-")
    m1 <- PercentageFeatureSet(m1,pattern = "^(Rps|Rpl)",col.name = "percent.ribo")
    m1$percent.mt[is.na(m1$percent.mt)] <- 100
    m1$percent.ribo[is.na(m1$percent.ribo)] <- 100
    nf <- read.table(paste("DataTables/NuclearFractions_",id,sep = ""))
    nf2 <- nf[,2]
    names(nf2) <- nf[,1]
    m1[["intronic.read.fraction"]] <- nf2
    m1$intronic.read.fraction[is.na(m1$intronic.read.fraction)] <- 0
    m1[["genotype"]] <- ifelse(id%in%c("223616","223617","243692","269132"),"t-carrier","non-carrier")
    if(id%in%c("243691","269131")){m1[["genotype"]] <- "t-hom"}
    m1 <- subset(m1, subset = percent.mt < 5 & intronic.read.fraction > 0.4)

    m1 <- SCTransform(m1,vars.to.regress = c("percent.mt","percent.ribo","intronic.read.fraction"),
      variable.features.n = 10000)

Further, after [Germain, Sonrel, and Robinson 2020], we filtered out nuclei that were at least two median absolute deviations (MADs) away from the median of at least two of the following nuclei statitistics:

- log10 of total UMI count

- log10 of total gene count

- ratio of log10 total gene count and log10 total UMI count

- percent of counts coming from the top 20 most abundant genes (only the right tail was removed)

- percent mitochondrial UMIs

- percent ribosomal UMIs

```
    m1[["nFeaturePernCount"]] <- m1$nFeature_RNA/m1$nCount_RNA
    m1$nFeaturePernCount[is.na(m1$nFeaturePernCount)] <- 1
    top20Expr <- names(head(sort(rowSums(m1[["SCT"]]$data),decreasing = T),n=20))
```
    m1$percent.top20 <- PercentageFeatureSet(m1,features = top20Expr)
    m1$percent.top20[is.na(m1$percent.top20)] <- 100
`
