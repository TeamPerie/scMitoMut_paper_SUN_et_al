#!/bin/bash
#PBS -l walltime=30:0:00
#PBS -l mem=40gb
#PBS -l nodes=1:ppn=5
#PBS -q batch
#PBS -N cell_line
#PBS -o cell_line.log
#PBS -e cell_line.err

cd /data/users/wsun/projects/nature_barcode/analysis_3/run_example_10x_brain
source ~/.bashrc
#eval `modulecmd bash load igor`
export PS1=
mamba activate py3

lib_id=mgatk_atac
output=./tmp/atac
mkdir -p $output

bam_file=./data/cellranger_out/human_brain_3k_atac_possorted_bam.bam
cell_list=./data/cellranger_out/filtered_feature_bc_matrix/barcodes.tsv

mgatk tenx -i $bam_file \
  -g rCRS \
  -n $lib_id \
  -o $output -c 7 \
  -bt CB -b $cell_list \
  --keep-duplicates \
  --alignment-quality -1

echo "Subject: [Curie Bioinfo Cluster] mgatk is done" | sendmail sunwjie@gmail.com
