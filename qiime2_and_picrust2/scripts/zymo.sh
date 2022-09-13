#import zymo 16s reference database (homemade) 
qiime tools import --type 'FeatureData[Sequence]' --input-path zymo_16S_contig_ID.fasta --output-path zymo_16S_sequence.qza
qiime tools import --type 'FeatureData[Taxonomy]'  --input-path zymo_16S_contigs.txt --output-path zymo_16S_taxonomy.qza

#subset feature table and rep seqs
qiime feature-table filter-samples --i-table ../combined_dada2.qza --m-metadata-file ../metadata_rockcreek_16s_combined.txt --p-where "[type] = 'zymo'" --o-filtered-table ./zymo_table.qza
qiime feature-table filter-seqs --i-data ../combined_reps.qza --i-table zymo_table.qza --o-filtered-data zymo_reps.qza

#blast against the zymo database
qiime feature-classifier classify-consensus-blast --i-query zymo_reps.qza  --i-reference-reads ../../db/zymo_16S_sequence.qza --i-reference-taxonomy ../../db/zymo_16S_taxonomy.qza --o-classification zymo_taxonomy.qza

#blast with only saving the top hit
qiime feature-classifier classify-consensus-blast --i-query zymo_reps.qza  --i-reference-reads ../../db/zymo_16S_sequence.qza --i-reference-taxonomy ../../db/zymo_16S_taxonomy.qza --o-classification zymo_taxonomy.qza --p-maxaccepts 1
