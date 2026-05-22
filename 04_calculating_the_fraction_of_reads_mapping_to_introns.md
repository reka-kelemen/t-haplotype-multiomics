# Calculating the fraction of reads mapping to introns

We used the software DropletQC to output the intronic read fraction, a measure of nuclear, unspliced RNA for each nucleus:

    nuclear_fraction_tags(outs = "223616/outs",tiles = 1, cores = 1, verbose = TRUE)
