# Inferring genotype in the *t* complex for each haploid spermatid nucleus

We calculated a *t*-haplotype score for each nucleus in the haploid spermatid clusters of the +/*t* samples in the following way. Considering only those genes in the *t* complex, which showed expression either from their reference or their pseudo-*t*-alleles, we calculated the difference between the mean gene expression from the pseudo-*t*-alleles and the reference alleles:

    m.combined[["t.score"]] <- colMeans(m.combined[["RNA"]]$data[c(tgenes,"t-Ppp1cb","t-Smok-Tcr",
    "t-Rnpepl1"),])-colMeans(m.combined[["RNA"]]$data[homgenes,])

We fitted two normal distributions to the bimodal t-score distributions using the *mix* function of the *mixdist* R package, and identified the 90th percentile of the lower distribution as the highest limit for +-spermatids and the 10th percentile of the top distribution as the lowest limit for *t*-spermatids (see Supplementary figure 5).

    for (stage in levels(m.combined)[9:14]) {
    t.score.data <- m.combined$t.score[((m.combined$orig.ident%in%c("223616","223617"))+
    (m.combined$seurat_clusters==stage))==2]
    his <- hist(t.score.data,breaks=breaks,plot = F)  
    df <- data.frame(mid=his$mids, cou=his$counts)  
    fitpro <- mix(as.mixdata(df), mixparam(mu=c(-0.02,0.02), sigma=c(0.01,0.01)), dist="norm")
    nt.upper.limit <- qnorm(0.9,mean = fitpro$parameters$mu[1],sd = fitpro$parameters$sigma[1])
    t.lower.limit <- qnorm(0.1,mean = fitpro$parameters$mu[2],sd = fitpro$parameters$sigma[2])
    }

We computed differential expression statistics between the putative *t*-spermatids and +-spermatids in each cluster, and we retained *t* complex genes, whose reference or *t*-alleles were significantly differentially expressed with an adjusted p-value below 0.05:

    FindMarkers(m.combined,ident.1 = Cells(m.combined)[((m.combined$genotype=="t-carrier")+
                                                        (m.combined$seurat_clusters==stage)+
                                                        (m.combined$t.score>t.lower.limit))==3],
                             ident.2 = Cells(m.combined)[((m.combined$genotype=="t-carrier")+
                                                         (m.combined$seurat_clusters==stage)+
                                                         (m.combined$t.score<nt.upper.limit))==3],
                                                         assay = "RNA",test.use = "MAST")

Based only on the genes, whose alleles were significantly differentially expressed, we computed the *t*-score again:

    m.combined[["t.score.DE"]] <- colMeans(m.combined[["RNA"]]$data[unique(tgenes.DE),])-
    colMeans(m.combined[["RNA"]]$data[unique(homgenes.DE),])

We fitted again two normal distributions to the bimodal distributions of the new *t*-score within each cluster using the *mixdist* package in R, using the same code as above. Since the *t*-score distributions are less bimodal in the round spermatid stages than in the elongating and condensing spermatid stages, we used the 75th percentiles of the lower *t*-score distributions as the upper limits for being classified as +-spermatids, and the 25th percentiles of the upper *t*-score distributions as the lower limits for being classified as *t*-spermatids. For elongating and condensing spermatid clusters we used the 90th and 10th percentiles, respectively.
