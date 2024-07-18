#!/bin/bash
#PBS -l walltime=10:0:00
#PBS -l mem=50gb
#PBS -l nodes=1:ppn=15
#PBS -q batch
#PBS -N mquad
#PBS -o mquad.log
#PBS -e mquad.err


# cd /data/users/wsun/projects/nature_barcode/analysis_3/run_example_10xcnv

# source ~/.bashrc
#eval `modulecmd bash load igor`
# export PS1=

#mamba activate R
#Rscript ./run_mgatk2mtx.R

# mamba activate py3
# mkdir -p ./tmp/mquad_output
# mquad -m ./data/atac_out_bj_mkn45_1pct/atac_out_bj_mkn45_1pct_filtered.alt.mtx,./data/atac_out_bj_mkn45_1pct/atac_out_bj_mkn45_1pct_filtered.dp.mtx  -o ./tmp/mquad_output -p 10
mquad -m ./tmp/cellsnp_out/cellSNP.tag.AD.mtx,./tmp/cellsnp_out/cellSNP.tag.DP.mtx  -o ./tmp/mquad_output_cellsnp-lite -p 10 --minDP 5

# echo "Subject: [Curie Bioinfo Cluster] mquad is done" | sendmail sunwjie@gmail.com
