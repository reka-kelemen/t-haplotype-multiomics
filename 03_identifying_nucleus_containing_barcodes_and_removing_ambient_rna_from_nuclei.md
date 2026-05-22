# Identifying nucleus-containing barcodes and removing ambient RNA from nuclei

We used Cellbender to remove ambient RNA and identify usable barcodes:

    cellbender remove-background \
         --input  223616/outs/raw_feature_bc_matrix/ \
         --output ./223616_output_3 \
         --expected-cells 7000 \
         --epochs 150 --cuda --total-droplets-included 30000 --fpr 0.05 --learning-rate 0.00005 
         --empty-drop-training-fraction 0.3 --low-count-threshold 500
         --z-dim 50 --z-layers 200
