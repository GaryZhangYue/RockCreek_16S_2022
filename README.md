# RockCreek_16S_2022
Folder data_analysis_with_R has two files: rarefied_table_115_samples.R and unrarefied_asv_table_new_group_manually_saved.R
The first file records the R codes used for data analysis and visualization for all samples by grouping samples into three categories ('on', 'altered' and 'off').
The second file records the R codes used for data analysis and visualization for a subset of samples, which was grouped into two aeration phases ('on and high' and 'off or low' aeration).

Folder qiime2_and_picrust2 contains the bash scripts and the associated configuration files used for producing the ASV tables, taxonomy classification, phylogenetic tree, and PICRUSt2 functional gene abundance inference.
In the sub-directory script/:
qiime2_merge_sequence_2_runs.csh was used to merge the qiime2 feature table and representative sequences;
moving_picture_analysis.csh was used to call QIIME2 function to perform data analysis, including rarefication, creating a rooted phylogenetic tree, taxonomy classification,alpha and beta diversity analysis; some of the output files were used as input to the R data analysis files
picrust2_unstratified_new_group.csh was used to call picrust2 software to infer functional gene abundance

The configuration files associated with qiime2_merge_sequence_2_runs.csh and moving_picture_analysis.csh was included in the directory config_files
