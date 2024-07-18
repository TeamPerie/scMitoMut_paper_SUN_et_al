
bam_file="./data/cellranger_out/bj_mkn45_1pct_possorted_bam.bam"
barcode_list="./data/cellranger_out/bj_mkn45_1pct_per_cell_barcode.tsv"
chrM_ref="./data/annotation/chrM.fasta"

cellsnp-lite -s $bam_file \
    -O ./tmp/cellsnp_out \
    -b $barcode_list \
    -p 1 \
    -f $chrM_ref \
    --UMItag None \
    --countORPHAN \
    --exclFLAG UNMAP \
    --chrom MT
    # --minCOUNT 0 \
    # --minLEN  0 \
    # --minMAPQ 0 \
    # --minMAF 0 \
    
