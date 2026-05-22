# Assembling a *t*-spermatid transcriptome and a +/+ testis transcriptome

<img width="704" height="929" alt="image" src="https://github.com/user-attachments/assets/b100201f-00cf-427d-a5cb-62041a48c308" />


*Figure 1. Pipeline of assembling a t-spermatid and a +/+ testis transcriptome*

## Extracting reads of *t*-spermatids

To create a *t*-spermatid assembly, we extracted reads with barcodes that we inferred to be *t*-spermatids.

    umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN \
    --stdin 223616_S1_L004_R1_001.fastq.gz \
    --stdout 223616_t-spermatids_1.fastq \
    --read2-in 223616_S1_L004_R2_001.fastq.gz \
    --read2-out=223616_t-spermatids_2.fastq \
    --whitelist=barcodes_t-spermatids_223616

## Trimming reads

We trimmed the 5’ [template switch oligo](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/algorithms/overview#read-trimming) sequence from the second reads, the polyT tails from the 5’ ends of the first reads, and performed quality trimming for both reads with FQTrim (version 0.9.7).

    fqtrim -o trimmed.fq -5 AAGCAGTGGTATCAACGCAGAGTACATGGG -a 10 -q 10 223616_t-spermatids_1.fastq 
    223616_t-spermatids_2.fastq

We combined the trimmed second reads from the two +/*t* samples, and truncated the read names so that the assembler, *Trinity*, would accept the second reads as single-end reads.

    cat 22361*_t-spermatids_trimmed_2.fq > t-spermatids_TSOtrimmed_2.fq
    cut -f1 -d " " t-spermatids_TSOtrimmed_2.fq > t-spermatids_TSOtrimmed_nameChanged_2.fq

## Assembling transcriptomes

We used the software *Trinity* (version 2.15.1), to assemble transcripts at least 500 basepairs long:

    Trinity --seqType fq --single t-spermatids_TSOtrimmed_nameChanged_2.fq --CPU 15 --max_memory 400G 
    --min_contig_length 500 --output t-spermatids_read2_trinity_output

We filtered for contigs above 0.5 TPM in both samples of the focal genotype (+/*t* or +/+). We estimated expression levels for each contig in each sample using the software *Kallisto* (version 0.50.1):

    kallisto index -i TSpermatidAssembly tSpermatids_stranded_trinity_output.Trinity.fasta
    kallisto quant -i TSpermatidAssembly -o kallisto_outdir_223616 --single --fr-stranded -l 50 -s 10
     FASTQ_Files/223616_t-spermatids.trimmed_nameChanged_2.fq
    kallisto quant -i TSpermatidAssembly -o kallisto_outdir_223617 --single --fr-stranded -l 50 -s 10
     FASTQ_Files/223617_S2_L004_R2_001_t-spermatids.trimmed_nameChanged.fq
     
    ###### create matrix of TPM estimates per contig
    ~/Software/trinityrnaseq-v2.15.1/util/abundance_estimates_to_matrix.pl --est_method kallisto 
    --out_prefix kallisto --gene_trans_map tSpermatids_stranded_trinity_output.Trinity.fasta.gene-trans-map --name_sample_by_basedir kallisto_outdir_223616/abundance.tsv kallisto_outdir_223617/abundance.tsv

    #### Filter for contigs with TPMs above 0.5 in both samples #############
    samtools faidx tSpermatids_stranded_trinity_output.Trinity.fasta 
    -r <(awk '$2>0.5&&$3>0.5 {print $1}' kallisto.isoform.TPM.not_cross_norm) > tSpermatidsAssembly_min05TPMAllSamples.fa

Next, we reduced the redundancy of our assembly by removing the shorter contig in pairs that have less than a 1% divergence over an alignment of at least a 100 basepairs. This was achieved using our own custom script, called SpliceFinder, which uses the software BLAT to produce alignments between transcripts:

    perl ./SpliceFinder_2_pblat.pl tSpermatidsAssembly_min05TPMAllSamples.fa

## Combining the *t*-spermatid assembly and the +/+ testis assembly

One of the +/+ testis samples (sample ID 223618) was used to assemble a +/+ transcriptome using the same steps and code as used for the *t*-spermatid transcriptome, above (without the first read-extraction step). We merged the two assemblies, and removed redundant contigs by discarding the shorter transcripts in pairs that have less than a 1% divergence from each other over alignments of at least a 100 basepairs, using our custom SpliceFinder script.

    sed 's/>/>t_/g' tSpermatidsAssembly_min05TPMAllSamples.fa.long > t_assembly.fa
    sed 's/>/>nt_/g' ntTestisAssembly_min05TPMAllSamples.fa.long > nt_assembly.fa
    cat t_assembly.fa nt_assembly.fa > t_nt_combined.fa 

    perl ./SpliceFinder_2_pblat.pl t_nt_combined.fa

## Annotating the contigs

We mapped the contigs against the mouse transcriptome (GRCm39) and the *t*-haplotype-specific sequences, *Ppp1cb<sup>t</sup>, Rnpepl1<sup>t</sup>* and *Smok<sup>Tcr</sup>*, using the tool BLAT (version 36x2).

    blat -out=blast8 Mus_musculus.GRCm39.cdna_and_ncrna_tSpecificTranscripts.fa 
    t_nt_assembly.fa.long TranscriptomeHits.blast

Then we assigned contigs to the gene whose transcript they had their best BLAST score against, using our custom script *besthitblat.pl*, and a minimum alignment length of 200 basepairs. Contigs with an antisense alignment were annotated as the gene’s antisense version (e.g. GeneA-as).

    perl besthitblat_blastScore.pl TranscriptomeHits_sorted.blast
    join -1 2 <(awk '$4>199' TranscriptomeHits_sorted.blast.besthit | sort -V -k2,2) 
    <(sort -V -k1,1 TranscriptID_GeneName_202406) | awk '{ if($9>$10){print $13"-as",$2}
    else {print $13,$2} }' | tr ' ' '\t' | sort -k1,1 > Gene_Contig
