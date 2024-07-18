
bam_file="./data/cellranger_out/human_brain_3k_atac_possorted_bam.bam"
barcode_list="./data/cellranger_out/filtered_feature_bc_matrix/barcodes.tsv"
chrM_ref="./data/annotation/chrM.fasta"

cellsnp-lite -s $bam_file \
    -O ./tmp/cellsnp_out \
    -b $barcode_list \
    -p 7 \
    -f $chrM_ref \
    --UMItag None \
    --countORPHAN \
    --exclFLAG UNMAP \
    --chrom chrM
    # --minCOUNT 0 \
    # --minLEN  0 \
    # --minMAPQ 0 \
    # --minMAF 0 \
    
