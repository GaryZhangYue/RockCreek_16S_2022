#Creator: Yue Zhang
#The function of this script is to:
#1) parse the results of the best BLAST of the positive controls (output of zymo.sh).
#2) document the zymo community theoretical composition

#load required packages
packages_to_load = c('stringr','phyloseq','vegan','tidyr','ggplot2','qiime2R',
                     'dplyr','ANCOMBC','usedist','rioja','ggpubr','spaa','ape',
                     'cowplot','ggforestplot','GGally','scales','metagMisc','ggsci','rstatix','wesanderson')
lapply(packages_to_load,require, character.only = T)

col21 <- rev(c("tomato1","darkblue","turquoise1","lightblue","darkred","mediumblue","purple","bisque","greenyellow","yellow","violetred2","darkgreen","darkgoldenrod1","deeppink3","cadetblue4","orchid2","seagreen3","purple4","dodgerblue2","red","gray27"))

#read the best BLAST results
pr<-qza_to_phyloseq(
  features="R_analysis_input/zymo_table.qza",
  taxonomy="R_analysis_input/zymo_taxonomy_maxxacepts1.qza")
taxa = as.data.frame(tax_table(pr))
asv = as.data.frame(otu_table(pr))

df = merge(x = subset(taxa, select =  c(Kingdom), drop = F),
           y = asv,
           by = 0,
           all.y = T)
df = df[,-1]
df = aggregate(.~Kingdom, df, sum)
#rename
names(df)[1] = 'strain'
#input the theoretical composition of 16s copy number
df.theoretical = data.frame(strain = c('Pseudomonas_aeruginosa_16S_170923.fasta',
                                       'Escherichia_coli_16S_170923.fasta',
                                       'Salmonella_enterica_16S_170923.fasta',
                                       'Lactobacillus_fermentum_16S_170924.fasta',
                                       'Enterococcus_faecalis_16S_170923.fasta',
                                       'Staphylococcus_aureus_16S_170923.fasta',
                                       'Listeria_monocytogenes_16S_170923.fasta',
                                       'Bacillus_subtilis_16S_170923.fasta'),
                            theoretical_percentage = c(4.2,10.1,10.4,18.4,9.9,15.5,14.1,17.4))

sum(df.theoretical$theoretical_percentage) #sanity check
df = merge(x = df, y = df.theoretical, by = 'strain')
df[,c(2:ncol(df))] = decostand(df[,c(2:ncol(df))],2,method = 'total')
colSums(df[,c(2:ncol(df))])
#stack columns
df.s = gather(df, "sample",'frequency',-strain)
df.s$strain = sub('_16S.*','',df.s$strain)

write.csv(df,'R_analysis_input/zymo.csv')
