
ml qiime2

#collapse feature table to taxonomy level-6
qiime taxa collapse --i-table combined_dada2.qza --p-level 6 --i-taxonomy depth-94840-_taxonomy.qza --o-collapsed-table combined_dada2_collapsed.qza
#calculate relative frequency for the collapsed table
qiime feature-table relative-frequency --i-table combined_dada2_collapsed.qza --o-relative-frequency-table combined_dada2_collapsed_frequency.qza

#convert .qza to .biom format
qiime tools export --input-path combined_dada2_collapsed_frequency.qza --output-path collapsed_feature_table/

#create a new folder and mv the output file to the folder. change names as well
#convert .biom to desired .txt format
biom convert -i combined_dada2_collapsed_frequency.biom -o combined_dada2_collapsed_frequency.txt --to-tsv --header-key taxonomy

