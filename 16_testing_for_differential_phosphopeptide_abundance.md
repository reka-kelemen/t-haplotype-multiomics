# Testing for differential phosphopeptide abundance

We compared phosphopeptide abundances using the PhosR package in R. First we created a table of normalized protein abundances:

    library(msqrob2)
    pe <- readQFeatures(assayData=peptideTableLog,colData=columnInfo,name="pepLog")
    rowData(pe[["pepLog"]])$nNonZero <- rowSums(assay(pe[["pepLog"]]) > 0)
    Protein_filter <- rowData(pe[["pepLog"]])$Protein.Group %in%
     smallestUniqueGroups(rowData(pe[["pepLog"]])$Protein.Group)
    pe <- pe[Protein_filter,,]
    pe <- filterFeatures(pe, ~ nNonZero >= 2)
    pe <- normalize(pe,i = "pepLog",name = "pepNormQuantilesRobust",method = "quantiles.robust")
    pe <- aggregateFeatures(pe,i = "pepNormQuantilesRobust", fcol = "Protein.Group",name = "proteinQuantNormRob")
    protabund <- assay(pe[["proteinQuantNormRob"]])

Then we use this table to normalize phosphopeptide abundances:

    pp <- read.delim("report.phosphosites_90.tsv",sep = "\t",quote = "", comment.char = "",header = T)
    quant <- log2(apply(as.matrix(pp[,(ncol(pp)-9):ncol(pp)]), 2, as.numeric)) ## transform it to log2 base values
    ppe <- PhosphoExperiment(assays = list(Quantification = quant),UniprotID=pp$Protein,GeneSymbol=pp$Gene.Names,
                             Sequence=pp$Sequence,Residue = pp$Residue,Site = pp$Site)
    SummarizedExperiment::colData(ppe) <- sample.info[colnames(quant),]
    gntype <- SummarizedExperiment::colData(ppe)$genotype
    ppe <- selectGrps(ppe, gntype, 1,n=1,assay = "Quantification")
    ppe <- medianScaling(ppe, scale = F, assay = "Quantification")
    ppe.scaled <- ppe@assays@data$scaled
    ppe.scaled.protabund <- merge(ppe.scaled,protabund,by.x = "DP.uniprots",by.y = 0,all.x = T,
    suffixes = c("phos","prot"))
    ppe.scaled.protabund[is.na(ppe.scaled.protabund)] <- 0
    ppe.scaled.protnorm <- matrix(as.numeric(unlist(ppe.scaled.protabund[,2:11])),
    nrow=nrow(ppe.scaled.protabund[,2:11]))-matrix(as.numeric(unlist(ppe.scaled.protabund[,
    15:24])),nrow=nrow(ppe.scaled.protabund[,15:24]))
    for (i in 1:nrow(ppe.scaled.protnorm)) {
      temprow <- ppe.scaled.protnorm[i,]
      temprow[temprow==0] <- min(temprow)-3
      ppe.scaled.protnorm[i,] <- temprow
    }

We use the limma paackage to test for a significant genotype effect on each phosphopeptide abundance:

    design <- model.matrix(~ gntype)
    fit <- lmFit(ppe.scaled.protnorm, design)
    contrast.matrix <- makeContrasts(gntypet, levels=design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2,robust = T)
    DP <- topTable(fit2, coef="gntypet", number = Inf,adjust.method = "BH")

Only phosphopeptides with non-overlapping intensity values between genotype groups and with p-value below 0.05 were selected:

    tsmaller <- rownames(ppe.scaled.protnorm)
    [rowMaxs(ppe.scaled.protnorm[,gntype=="t"])<rowMins(ppe.scaled.protnorm[,gntype=="nt"])]
    tbigger <- rownames(ppe.scaled.protnorm)
    [rowMaxs(ppe.scaled.protnorm[,gntype=="nt"])<rowMins(ppe.scaled.protnorm[,gntype=="t"])]
    nonoverlapping <- c(tsmaller,tbigger)
    DP.nonoverlapping <- DP[nonoverlapping,]
    DP.nonoverlapping.sign <- DP.nonoverlapping[DP.nonoverlapping$P.Value<0.05,]
