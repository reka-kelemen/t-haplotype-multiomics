# Mapping and counting of UMIs (unique molecular identifiers)

CellRanger was used to map and count the reads for each gene and barcode.

    cellranger mkref --genome GRCm39_plusPseudo-t 
    --fasta Mus_musculus.GRCm39.dna.primary_assembly_Psuedo-t_GainedGenes.fa 
    --genes Mus_musculus.GRCm39.110_andPseudo-t_andTcrPpp1cbRnpepl1.gtf

    cellranger count --include-introns 
    --fastqs HTK5FDSX5_4_R15092_20230319/demultiplexed/223616 --id 223616 
    --transcriptome GRCm39_plusPseudo-t
