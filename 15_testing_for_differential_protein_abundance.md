# Testing for differential protein abundance

We predicted protein sequences from a published *t*-haplotype sequence ([Swanepoel et al. 2025](references.md#ref-swanepoel2025acquisition)):

    samtools faidx tw5.fasta chr_tw5:1-45000000 > tw5_proximal45Mbs.fa
    miniprot -t100 -d tw5.mpi tw5_proximal45Mbs.fa #index
    miniprot -t100 --gff tw5.mpi UP000000589_10090.plusSmokTcrt-Ppp1cb.fasta > tw5.gff ##detect mouse proteins
    ## output translated amino acid sequence
    gffread tw5.gff -g tw5_proximal45Mbs.fa -M -K -E -V -H -C -y tw5_translated_cds.fa 

And kept only those proteins that cover at least 90 percent of a mouse protein:

    blat -prot -out=psl UP000000589_10090.plusSmokTcrt-Ppp1cb.fasta tw5_translated_cds.fa Hits.blat
    sort -k10 Hits.blat > Hits_sorted.blat
    perl besthitblat.pl Hits_sorted.blat

    tail +2 Hits_sorted.blat.besthit | awk '(($1+$2)/$15)>0.9 {print $10,$14}' | tr "|" " " | 
    cut -f1,3 -d " " | tr " " "\t" > tw5protein_besthit_UniProtProtein_90percAlnCoverage

We used DiaNN to predict peptide abundance from the mouse proteome and the predicted *t*-haplotype proteins:

    diann --f ../../RawData/1-A14463-wildtype_Slot1-1_1_7274.d \
          --f ../../RawData/2-A14464-mutant_Slot1-2_1_7288.d \
          --f ../../RawData/3-A14465-mutant_Slot1-3_1_7280.d \
          --f ../../RawData/4-A14466-wildtype_Slot1-4_1_7286.d \
          --f ../../RawData/5-A14467-mutant_Slot1-5_1_7276.d \
          --f ../../RawData/6-A14468-wildtype_Slot1-6_1_7290.d \
          --f ../../RawData/7-A14469-mutant_Slot1-7_1_7284.d \
          --f ../../RawData/8-A14470-wildtype_Slot1-8_1_7282.d \
          --f ../../RawData/9-A14471-wildtype_Slot1-9_1_7278.d \
          --f ../../RawData/10-A14472-mutant_Slot1-10_1_7272.d \
          --lib  --threads 64 --verbose 1 --out report.parquet --qvalue 0.01 --matrices --out-lib lib.parquet \
          --gen-spec-lib --predictor \
          --fasta uniprotkb_proteome_UP000000589_2025_11_14.fasta --fasta tw5_predicted_proteins.fa \
          --fasta-search --min-fr-mz 200 --max-fr-mz 1800 --met-excision --min-pep-len 7 --max-pep-len 30 \
          --min-pr-mz 400 --max-pr-mz 1350 --min-pr-charge 2 --max-pr-charge 4 --cut K*,R* \
          --missed-cleavages 1 --unimod4 --var-mods 3 --var-mod UniMod:35,15.994915,M \
          --var-mod UniMod:1,42.010565,*n --var-mod UniMod:21,79.966331,STY \
          --mass-acc 15 --mass-acc-ms1 15 --peptidoforms --reanalyse --rt-profiling

We compared protein abundances using the msqrob2 package in R:

    library(msqrob2)
    pe <- readQFeatures(assayData=peptideTableLog,colData=columnInfo,name="pepLog")
    rowData(pe[["pepLog"]])$nNonZero <- rowSums(assay(pe[["pepLog"]]) > 0)
    Protein_filter <- rowData(pe[["pepLog"]])$Protein.Group %in% 
    smallestUniqueGroups(rowData(pe[["pepLog"]])$Protein.Group)
    pe <- pe[Protein_filter,,]
    pe <- filterFeatures(pe, ~ nNonZero >= 2)

    pe <- normalize(pe,i = "pepLog",name = "pepNorm",method = "center.median")
    pe <- aggregateFeatures(pe,i = "pepNorm", fcol = "Protein.Group",name = "protein")

    pe <- msqrob(object = pe, i = "protein", formula = ~genotype)
    getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])
    L <- makeContrast("genotypet=0", parameterNames = c("genotypet"))
    pe <- hypothesisTest(object = pe, i = "protein", contrast = L)
    res <- rowData(pe[["protein"]])$genotypet
    res[res$adjPval<0.05,]
