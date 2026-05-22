# Testing expression evolution in spermatocytes along the *Mus* phylogeny

We used the evemodel package (version 0.0.0.9008) in R to test for lineage-specific shifts in expression for dyneins in the house mouse genome. We provided the species tree, as estimated by [Chevret, Veyrunes, and Britton-Davidian 2005]; [López-Antoñanzas et al. 2024], and used the expression estimates in spermatocytes from the study of [Kopania et al. 2022].

    phy.string <- "(((Mm:2,Ms:2):5,Mp:7):2.9,Rn:9.9);"
    tre <- read.tree(text = phy.string)
    thetaShiftBool_full <- ifelse(1:Nedge(tre)%in%c(3),T,F) ## Mm
    test.thetaShift_full <- twoThetaTest(tree = tre, gene.data = table.of.expr.values,
    isTheta2edge = thetaShiftBool_full, colSpecies = species.names.of.columns)
