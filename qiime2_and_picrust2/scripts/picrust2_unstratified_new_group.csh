#!/bin/sh
#
####SBATCH
#SBATCH --job-name=picrust2
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --partition=defq

#load picrust2 environment
ml anaconda
conda activate picrust2

picrust2_pipeline.py -s ../PICRUSt2_new_group/combined_reps.fasta -i ../PICRUSt2_new_group/unrarefied_asv_table_new_group.txt -o ../PICRUSt2_new_group/unstratified_new_group -p 48
