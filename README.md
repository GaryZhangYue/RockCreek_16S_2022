# RockCreek_16S_2022

This repository is associated with the manuscript of "Zhang, Y. and Preheim, S.P., 2022. Influence of artificial aeration on the composition, diversity and potential function of microbial communities in a hypoxic estuarine system. Applied and Environmental Microbiology, in preparation". The related data will be published soon.

Directory data_analysis_with_R has two files: 1) rarefied_table_115_samples.R and 2) unrarefied_asv_table_new_group_manually_saved.R
The first file records the R codes used for data analysis and visualization for all samples by grouping them into three categories ('on', 'altered' and 'off').
The second file records the R codes used for data analysis and visualization for a subset of samples, which was grouped into two aeration phases ('on and high' and 'off or low' aeration).

Directory qiime2_and_picrust2 contains the bash scripts and the associated configuration files for producing the ASV tables, taxonomy classification, phylogenetic tree, and PICRUSt2 functional gene abundance inference.

In the sub-directory script/:
qiime2_merge_sequence_2_runs.csh was used to merge the qiime2 feature tables and representative sequences;
moving_picture_analysis.csh was used to call QIIME2 function to perform data analysis, including rarefaction, generation of a rooted phylogenetic tree and taxonomy classification, as well as alpha and beta diversity analysis; some of the output files were used as input to the R data analysis files
picrust2_unstratified_new_group.csh was used to call picrust2 to infer functional gene abundance

The configuration files associated with qiime2_merge_sequence_2_runs.csh and moving_picture_analysis.csh were included in the directory config_files/
