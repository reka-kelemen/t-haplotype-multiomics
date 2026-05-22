# Removing ambient RNA contamination

We used the software CellBender (version 0.3.2) to remove putative ambient RNA from each nucleus, and to remove barcodes that likely represent droplets that do not contain a nucleus.

    cellbender remove-background \
         --input  223616/outs/raw_feature_bc_matrix/ \
         --output ./223616_cellbender_output \
         --expected-cells 7000 \
         --epochs 150 --cuda --total-droplets-included 30000 --fpr 0.05 --learning-rate 0.00005 
         --empty-drop-training-fraction 0.3 --low-count-threshold 500 \
         --z-dim 50 --z-layers 200
