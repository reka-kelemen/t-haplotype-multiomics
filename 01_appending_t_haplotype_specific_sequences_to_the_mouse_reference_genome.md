# Appending *t*-haplotype-specific sequences to the mouse reference genome

We created a pseudo-reference sequence for the *t*-haplotype by substituting *t*-haplotype-specific single nucleotide and INDEL variants, identified in [Kelemen and Vicoso 2018], into the *Mus musculus* reference genome (GRCm39).

    bcftools index ERR899402_Tcomplex_tHaplSubset_ConvertedToGRCm39.vcf.gz
    samtools faidx Mus_musculus.GRCm39.dna.primary_assembly.fa 17:5000000-40000000 >
     Mus_musculus.GRCm39.dna.primary_assembly_tcomplex.fa
    bcftools consensus -H 2 -f Mus_musculus.GRCm39.dna.primary_assembly_tcomplex.fa -o Pseudo-t-haplotype.fa
     ERR899402_Tcomplex_tHaplSubset_ConvertedToGRCm39.vcf.gz

We concatenated the GRCm39 reference genome with the above generated pseudo-*t*-haplotype sequence, and the trascript sequences of two *t*-haplotype-specific gene duplicates, *Ppp1cb<sup>t</sup>* and *Rnpepl1<sup>t</sup>*, identified in [Kelemen et al. 2022]. We also appended *t*-haplotype-specific *Smok<sup>Tcr</sup>* sequence identified by [Herrmann et al. 1999].
