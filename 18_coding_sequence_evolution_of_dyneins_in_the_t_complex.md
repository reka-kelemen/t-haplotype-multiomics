# Coding sequence evolution of dyneins in the *t* complex

We predicted dynein coding sequences on the *t*-haplotype sequence of [Swanepoel et al. 2025]:

    miniprot -t100 -d tw5.mpi tw5_proximal45Mbs.fa #index
    miniprot -t100 --gff tw5.mpi Mus_musculus.GRCm39.pep.all.plusSmokTcrt-Ppp1cb.fa > tw5.gff
    gffread tw5.gff -g tw5_proximal45Mbs.fa -M -K -E -V -H -C -x tw5_cds.fa ## output spliced CDSs

We estimated a phylogenetic tree from the alignments with IQTree with 1000 bootstrap iterations:

    iqtree2 -s alignment.fa -B 1000

    seqfile = alignment.fa * sequence data file name
    outfile = PAML_Results * main result file name
    treefile = tree_base.newick * tree structure file name
    noisy = 3  * 0,1,2,3,9: how much rubbish on the screen
    verbose = 1  * 1: detailed output, 0: concise output
    runmode = 0  * 0: user tree
    seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table *

    ndata = 1
    clock = 0   * 0:no clock, 1:clock; 2:local clock
    aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a, 7:AAClasses
    aaRatefile = wag.dat * only used for aa seqs with model=empirical(_F)
    model = 0   * models for codons: 0:one, 1:b, 2:2 or more dN/dS ratios for branches
    NSsites = 0  * 0:one w;1:neutral;2:selection; 
    icode = 0  * 0:universal code; 1:mammalian mt; 2-11:see below
    Mgene = 0  * 0:rates, 1:separate;
    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate
    omega = 0.5  * initial omega value
    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
    alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)
    Malpha = 0  * different alphas for genes
    fix_rho = 1  * 0: estimate rho; 1: fix it at rho
    rho = 0. * initial or fixed rho,   0:no correlation
    getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
    RateAncestor = 0  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
    Small_Diff = .5e-6 *
    * cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)? *
    * fix_blength = 0   * 0: ignore, -1: random, 1: initial, 2: fixed
    method = 0   * 0: simultaneous; 1: one branch at a time

For the branch models we changed the parameter model=2, and for the branch-site models the parameter NSsites=2.
