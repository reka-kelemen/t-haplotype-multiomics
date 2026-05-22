# Mapping to and counting UMIs for each contig in the combined assembly

We used the whole set of assembled contigs as the reference to map the single nucleus RNA-seq reads against using the software CellRanger (version 7.1.0).

    cellranger mkref --genome Assembly --fasta t_nt_assembly.fa.long.fa --genes Assembly.gtf
    cellranger count --include-introns true --fastqs HTK5FDSX5_4_R15092_20230319/demultiplexed/223616 
    --id 223616 --transcriptome Assembly
