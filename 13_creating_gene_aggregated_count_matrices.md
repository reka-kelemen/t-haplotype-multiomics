# Creating gene-aggregated count matrices

The count matrix produced by CellBender for each sample was read into R, and the counts of contigs mapping to the same gene (sense and antisense mappings were treated as two different genes) were added up:

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
