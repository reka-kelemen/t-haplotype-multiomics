# Differential expression analysis

We conducted tests of differential expression with the MAST method in Seurat, which was designed to take into account the properties of data from single cell experiments. For each diploid cell type we compared +/*t* vs. +/+ samples, and passed the replicate and batch information to MAST as the latent variables.

    m.combined <- PrepSCTFindMarkers(m.combined)
    markers <- FindMarkers(m.combined, subset.ident = "Spermatocytes 1",
    recorrect_umi=F,group.by = "genotype", ident.1 = "t-carrier",ident.2 = "non-carrier",
    min.pct=0.25,logfc.threshold=0.1,test.use = "MAST",latent.vars = c("replicate","batch"))

For haploid spermatids we grouped nuclei based on the *t*-score-based inferred cell genotypes and conducted different comparisons:

    markers <- FindMarkers(m.combined, subset.ident = "Early round spermatids",
    recorrect_umi=F,group.by = "cell.genotype.t.score.based", ident.1 = "t",ident.2 = "nt",
    min.pct=0.25,logfc.threshold=0.1,test.use = "MAST",latent.vars = c("replicate","batch"))
