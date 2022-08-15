#!/bin/sh
#
###SBATCH
#SBATCH --job-name=merge
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --partition=lrgmem
#import qiime environment
module load qiime2/2018.8

#make sure to edit the merge_config.txt and put into the combined directory
source ../config_files/qiime2_merge_sequence.config

#echo the time for each
echo "Starting qiime2 two file merge"
date

#Check TMPDIR path is working
echo $TMPDIR
qiime feature-table merge --i-tables ${RUN1_DADA2} --i-tables ${RUN2_DADA2} --o-merged-table ${PREFIX}_dada2.qza
qiime feature-table merge-seqs --i-data ${RUN1_REPS} --i-data ${RUN2_REPS} --o-merged-data ${PREFIX}_reps.qza

qiime feature-table summarize --i-table ${PREFIX}_dada2.qza --o-visualization ${PREFIX}_dada2_summarize.qzv --m-sample-metadata-file ${METADATA}

echo "End of Script"
date

