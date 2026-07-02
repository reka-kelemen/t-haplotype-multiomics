# Creating gene-aggregated count matrices

library(Seurat)
library(hdf5r)
The count matrix produced by CellBender for each sample was read into R, and the counts of contigs mapping to the same gene (sense and antisense mappings were treated as two different genes) were added up. Only the previously identified high quality nuclei were retained:

    id <- 223616
    mouse.data <-  Read10X_h5(paste("../",id,"_cellbender_output_filtered.h5",sep = ""), use.names = TRUE)
      colnames(mouse.data) <- paste(i,colnames(mouse.data),sep = "_")
      mouse <- mouse.data[,colnames(mouse.data)%in%Cells(m.combined.ref)[m.combined.ref$orig.ident==id]]
      gene.aggregated.matrix <- c()
      for (gene in unique(gene.contig[,1])) {
        contigsToTake <- gene.contig[gene.contig[,1]==gene,2]
        gene.aggregated.matrix <- rbind(gene.aggregated.matrix,colSums(rbind(rep(0,ncol(mouse)),
        mouse[contigsToTake,])))
      }
      rownames(gene.aggregated.matrix) <- unique(gene.contig[,1])
      mouse <- CreateSeuratObject(counts = gene.aggregated.matrix, project = paste(id))
      mouse[["genotype"]] <- gntype

Samples were merged into a single Seurat object, and UMI counts were normalized by the method SCTransform:

```
library(sctransform)

m1 <- readRDS("mouse_223616_geneAggregated.rds")
m2 <- readRDS("mouse_223617_geneAggregated.rds")
m3 <- readRDS("mouse_223618_geneAggregated.rds")
m4 <- readRDS("mouse_272100_geneAggregated.rds")
mc <- merge(m1,list(m2,m3,m4),collapse=F,merge.dr=F,merge.data=F)

nVarFeatures <- 5000
m.combined.ref <- readRDS("Reference_based_object.rds")

mc[["percent.mt"]] <- PercentageFeatureSet(mc, pattern = "^mt-")
mc[["percent.ribo"]] <- PercentageFeatureSet(mc,pattern = "^(Rps|Rpl)")
mc[["intronic.read.fraction"]] <- m.combined.ref$intronic.read.fraction
mc$percent.mt[is.na(mc$percent.mt)] <- 100
mc$percent.ribo[is.na(mc$percent.ribo)] <- 100
mc$intronic.read.fraction[is.na(mc$intronic.read.fraction)] <- 0

mc$batch <- ifelse(mc$orig.ident==4,2,1)
mc <- subset(mc, subset = nCount_RNA > 0) ## there are nuclei with 0 counts that cause problems
mc <- SCTransform(mc,vars.to.regress = c("percent.mt","percent.ribo","intronic.read.fraction","replicate","batch"),variable.features.n = nVarFeatures)
saveRDS(mc, file = "m.merged.assemblyMapped.geneAggregated.rds")

```

