#Creator: Yue Zhang
#This scripts is used to analyze the subset of 85 samples by grouping them into 2 aeration phases -- aerators "on and high" and "off or low"

##########Preparation###############
#load required packages
packages_to_load = c('stringr','phyloseq','vegan','tidyr','ggplot2','qiime2R',
                     'dplyr','ANCOMBC','usedist','rioja','ggpubr','spaa','ape',
                     'cowplot','ggforestplot','GGally','scales','metagMisc','ggsci','rstatix',
                     'wesanderson','iCAMP')
lapply(packages_to_load,require, character.only = T)

#generating color palette
col1 = wes_palette("Rushmore1")
col2 = wes_palette("GrandBudapest2")
col3 = wes_palette('Darjeeling2')
col21 <- rev(c("tomato1","darkblue","turquoise1","lightblue","darkred","mediumblue","purple","bisque","greenyellow","yellow","violetred2","darkgreen","darkgoldenrod1","deeppink3","cadetblue4","orchid2","seagreen3","purple4","dodgerblue2","red","gray27"))

################load variables################
PFIG = FAS = 'R_analysis_output/85_samples'
#dir.create(PFIG)
#set random seed
set.seed(666)

############metadata pre-processing###################
metadata = read.csv('R_analysis_input/metadata_simplified_115.csv')
#keep the column type row
metadata_typerow = metadata[1,]
metadata = metadata[-1,]
#create a new identifier column in metadata
metadata$date = as.character(metadata$date)
metadata$date[metadata$date == '8'] = '08' 
metadata$date[metadata$date == '9'] = '09'
metadata$`depth.station.date.time` = paste(metadata$depth.watercolumn,metadata$station,
                                           metadata$date,metadata$time,sep = '-')

metadata$date_time = as.character(metadata$date_time)
metadata$date_time[metadata$date_time == '8.5'] = '08.5'
metadata$date_time[metadata$date_time == '9'] = '09'
metadata$date_time[metadata$date_time == '9.5'] = '09.5'

#change "t" to "top" and "b" to "bottom"
metadata$depth.watercolumn = as.character(metadata$depth.watercolumn)
metadata$depth.watercolumn[metadata$depth.watercolumn == 'b'] = 'bottom'
metadata$depth.watercolumn[metadata$depth.watercolumn == 't'] = 'top'
metadata$depth.watercolumn = factor(metadata$depth.watercolumn,levels = c("top", "bottom"))

#assign groups
metadata$date_numeric = as.numeric(metadata$date)
summary(metadata$date_numeric)
metadata %>% mutate(group = ifelse(date_numeric<=9 | date_numeric >= 26, 'on and high',
                                   ifelse(date_numeric >= 16 & date_numeric <= 23, 'off or low', 'transition'))) -> metadata
#remove the transitioning phase
metadata = metadata[metadata$group != "transition",]

#combine with the column type row
metadata = bind_rows(metadata,metadata_typerow)
metadata$group[nrow(metadata)] = 'categorical'
metadata$date_numeric[nrow(metadata)] = 'numeric'
metadata$date[nrow(metadata)] = 'categorical'
metadata$date_time[nrow(metadata)] = 'categorical'
metadata$depth.station.date.time[nrow(metadata)] = 'categorical'
metadata = metadata[c(nrow(metadata),1:nrow(metadata)-1),]

write.table(metadata,paste(FAS,'/metadata_85.txt',sep = ''), sep = '\t', quote = F,row.names = F)
metadata = NA

#####################read data using phyloseq#####################
#Creating a Phyloseq Object
pr<-qza_to_phyloseq(
  features="R_analysis_input/rarefied_table.qza",
  tree="R_analysis_input/rooted-tree.qza",
  taxonomy="R_analysis_input/taxonomy.qza",
  metadata=paste(FAS,'/metadata_85.txt',sep = ''))
pr
asv = data.frame(otu_table(pr)) 
metadata= data.frame(sample_data(pr))
tax = data.frame(tax_table(pr))

#set the depth column to factor
metadata$depth.watercolumn = factor(metadata$depth.watercolumn,levels = c('top','bottom'))

####################ASV table pre-processing##############################
#order the asv table based on the number of total counts across all samples
asv$sum = apply(asv,1,sum)
asv.o = asv[order(asv$sum,decreasing = T),]
head(asv.o)

#remove the ghost ASVs
asv.o %>% dplyr::filter(sum >0) %>% as.data.frame -> asv.of 

#simplify the ASV ID
nums <- sprintf('%0.5d', 1:nrow(asv.of)) 
ASVnames <- paste("ASV",nums, sep="")
ASVconversion <- data.frame(ASV_ID = rownames(asv.of),ASVnames, stringsAsFactors = F) 
rownames(asv.of) = ASVconversion$ASVnames
write.csv(x = ASVconversion, file = paste(FAS,"/ASVconversion.csv", sep = ''),
          row.names = T)

#replace the column names in asv.of with a meaningful identifier "depth.station.date.time"
asv.of = asv.of[,order(names(asv.of))]
metadata.o = metadata[order(rownames(metadata)),]
metadata.or = metadata.o
rownames(metadata.or) <- metadata.or$depth.station.date.time
SIDconversion = data.frame(SIDasv = names(asv.of)[c(1:ncol(asv.of)-1)], SIDmetadata = rownames(metadata.o),SD = metadata.o$depth.station.date.time)
SIDconversion$SIDmetadata.conv = str_replace_all(SIDconversion$SIDmetadata,'-','.')
summary(SIDconversion$SIDasv == SIDconversion$SIDmetadata.conv) #sanity check -- all TRUE
asv.of = asv.of[,c(1:ncol(asv.of)-1)]#remove the sum column
(names(asv.of) <- SIDconversion$SD[match(names(asv.of),SIDconversion$SIDmetadata.conv)])
asv.of = asv.of[,order(names(asv.of))]

#feature count and total number of reads
sum(asv.of)

#calculate frequency -- normalize by sample
asv.off = decostand(asv.of,2,method = 'total')
colSums(asv.off) #col sum should be equal to 1

#############taxonomic table pre-processing################
#remove ghost ASVs
tax = subset(tax, rownames(tax) %in% ASVconversion$ASV_ID)
rownames(tax) = ASVconversion$ASVnames[match(rownames(tax),ASVconversion$ASV_ID)]
tax.o = tax[order(rownames(tax)),]

########################Taxonomic barplot########################
####################Phylum level bar plot##################
#aggregate ASV table to Phylum level
tax.o %>%
  dplyr::select(`Phylum`) %>%
  merge(x = ., y = asv.off, all.y = T, by = 0) -> asv.offp
rownames(asv.offp) <- asv.offp$Row.names
asv.offp = asv.offp[,-c(1)]

#change NA in Phylum column to unknown
asv.offp$Phylum[is.na(asv.offp$Phylum)] = 'Unknown'
asv.offp = aggregate(. ~ `Phylum`,asv.offp, sum)
summary(colSums(asv.offp[,-c(1)])) #sanity check
rownames(asv.offp) = asv.offp$Phylum
asv.offp = subset(asv.offp, select = -c(Phylum))

#retain the 9 most abundant Phyla
ft = 9

asv.offp = asv.offp[order(rowSums(asv.offp), decreasing = T),]  
asv.offps = bind_rows(asv.offp[c(1:ft),],colSums(asv.offp[c(c(ft+1):nrow(asv.offp)),]))
rownames(asv.offps)[c(ft+1)] = 'All others'
summary(colSums(asv.offps))

#stack the ordered asv table
rownames(asv.offps) -> asv.offps$Phylum
asv.offps.s = gather(asv.offps,"SD","frequency",-Phylum)

#merge the stacked table with a subset of metadata to include some categorization columns
metadata.o %>%
  dplyr::select(`depth.station.date.time`,`depth.watercolumn`,`group`,`date_time`,`station`) %>%
  merge(x = ., y = asv.offps.s, all.x = T, by.x = "depth.station.date.time", by.y = "SD") -> asv.offps.s

#change they phylum column to ordered factors for plotting
asv.offps.s$Phylum = factor(asv.offps.s$Phylum, levels = c(rownames(asv.offps)))
#make bar plots
#with legend
ggplot(data = asv.offps.s, mapping = aes(x = `date_time`, y = `frequency`,fill = Phylum)) +
  geom_bar(position = 'fill', stat="identity",color = 'black') +
  scale_fill_manual(values = col21) +
  facet_grid(depth.watercolumn~station) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.3,'cm'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10, face = 'bold', color = 'black'),
        axis.text.y = element_text(size=15, face = 'bold', color = 'black'),
        axis.title=element_text(size=18,face="bold"),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"),
        strip.text.y = element_text(size = 15, color = "black", face = "bold"),
        strip.background = element_rect(color="black", fill="grey95", size=1.5, linetype="solid"))+
  xlab("Date") + ylab('Relative Abundance') +
  guides(fill=guide_legend(nrow=2))
ggsave(paste(PFIG,"/",'Taxonomy_top_9_Phyla_barplot.png',sep = ""),width = 10, height = 5, device = 'png')

#without legend
ggplot(data = asv.offps.s, mapping = aes(x = `date_time`, y = `frequency`,fill = Phylum)) +
  geom_bar(position = 'fill', stat="identity",color = 'black') +
  scale_fill_manual(values = col21) +
  facet_grid(depth.watercolumn~station) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'none',
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.3,'cm'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10, face = 'bold', color = 'black'),
        axis.text.y = element_text(size=15, face = 'bold', color = 'black'),
        axis.title=element_text(size=18,face="bold"),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"),
        strip.text.y = element_text(size = 15, color = "black", face = "bold"),
        strip.background = element_rect(color="black", fill="grey95", size=1.5, linetype="solid"))+
  xlab("Date") + ylab('Relative Abundance') +
  guides(fill=guide_legend(nrow=2))
ggsave(paste(PFIG,"/",'Taxonomy_top_9_Phyla_barplot_no_legend.png',sep = ""),width = 10, height = 5, device = 'png')


########ASV bar plot for Cyanobacteria##########
tax.om <- tax.o
tax.om.cya = subset(tax.om, Phylum == 'Cyanobacteria')
tax.om.cya %>%
  dplyr::select(Class) %>%
  merge(x = ., y = asv.off, by = 0, all.x = T) -> tax.om.cya
tax.om.cya$Class = paste(tax.om.cya$Row.names, tax.om.cya$Class, sep = ' ')
rownames(tax.om.cya) <- tax.om.cya$Class
tax.om.cya = tax.om.cya[,-c(1:2)]
tax.om.cyaso = bind_rows(tax.om.cya[1:ft,],colSums(tax.om.cya[c(ft+1):nrow(tax.om.cya),]))
rownames(tax.om.cyaso)[c(ft+1)] <- 'Other Cyanobacteria'

#stack the ASV table
tax.om.cyaso$Class <- rownames(tax.om.cyaso)
tax.om.cyaso.s = gather(tax.om.cyaso,"SD","frequency",-Class)

#merge the stacked table with a subset of metadata to add some categorization columns
metadata.o %>%
  dplyr::select(`depth.station.date.time`,`depth.watercolumn`,`group`,`date_time`,`station`) %>%
  merge(x = ., y = tax.om.cyaso.s, all.x = T, by.x = "depth.station.date.time", by.y = "SD") -> tax.om.cyaso.s

#plot
#with legend
ggplot(data = tax.om.cyaso.s, mapping = aes(x = `date_time`, y = `frequency`,fill = Class)) +
  geom_bar(position = 'stack', stat="identity",color = 'black') +
  scale_fill_jco()+
  facet_grid(depth.watercolumn~station) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.3,'cm'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10, face = 'bold', color = 'black'),
        axis.text.y = element_text(size=15, face = 'bold', color = 'black'),
        axis.title=element_text(size=18,face="bold"),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"),
        strip.text.y = element_text(size = 15, color = "black", face = "bold"),
        strip.background = element_rect(color="black", fill="grey95", size=1.5, linetype="solid"))+
  xlab("Date") + ylab('Relative Abundance') +
  guides(fill=guide_legend(nrow=4))
ggsave(paste(PFIG,"/",'Cyanobacteria_barplot.png',sep = ""),width = 12, height = 7, device = 'png')

#no legend
ggplot(data = tax.om.cyaso.s, mapping = aes(x = `date_time`, y = `frequency`,fill = Class)) +
  geom_bar(position = 'stack', stat="identity",color = 'black') +
  scale_fill_jco()+
  facet_grid(depth.watercolumn~station) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'none',
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.3,'cm'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10, face = 'bold', color = 'black'),
        axis.text.y = element_text(size=15, face = 'bold', color = 'black'),
        axis.title=element_text(size=18,face="bold"),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"),
        strip.text.y = element_text(size = 15, color = "black", face = "bold"),
        strip.background = element_rect(color="black", fill="grey95", size=1.5, linetype="solid"))+
  xlab("Date") + ylab('Relative Abundance') +
  guides(fill=guide_legend(nrow=4))
ggsave(paste(PFIG,"/",'Cyanobacteria_barplot_no_legend.png',sep = ""),width = 10, height = 5, device = 'png')

##########################alpha diversity analysis#################################
#coerce the group column in metadata.o to factors
metadata.o$group = factor(metadata.o$group, levels = c('on and high', 'off or low'))

#shannon diversity matrix
alpha = read_qza("R_analysis_input/shannon_vector.qza")
alpha <- alpha$data %>% tibble::rownames_to_column("SampleID")
alpha = merge(x = alpha,
              y = metadata.o,
              by.x = 'SampleID',
              by.y = 0,
              all.Y = T)

(p1t=ggboxplot(subset(alpha, depth.watercolumn == 'top'), x = 'group', y = 'shannon',
               add = "jitter", fill = 'group', size = 1)+
    xlab("Phase") + ylab('Shannon index') + ggtitle('Top') +
    coord_cartesian(ylim = c(2,10))+
    theme(text = element_text(size = 16), legend.position = 'none')+
    scale_fill_lancet(alpha=0.8) +
    stat_compare_means(method = "wilcox.test",label.x = 0.7,label.y = 2,size = 6))

(p1b=ggboxplot(subset(alpha, depth.watercolumn == 'bottom'), x = 'group', y = 'shannon',
               add = "jitter", fill = 'group', size = 1)+
    xlab("Phase") + ylab('Shannon index') + ggtitle('Bottom')+
    coord_cartesian(ylim = c(2,10))+
    theme(text = element_text(size = 16), legend.position = 'none')+
    scale_fill_jama(alpha=0.8) +
    stat_compare_means(method = "wilcox.test",label.x = 0.7,label.y = 2,size = 6))

#Faith PD matrix
alpha = read_qza("R_analysis_input/faith_pd_vector.qza")
alpha <- alpha$data %>% tibble::rownames_to_column("SampleID")
alpha = merge(x = alpha,
              y = metadata.o,
              by.x = 'SampleID',
              by.y = 0,
              all.Y = T)

(p2t=ggboxplot(subset(alpha, depth.watercolumn == 'top'), x = 'group', y = 'faith_pd',
               add = "jitter", fill = 'group', size = 1)+
    xlab("Phase") + ylab('Faith PD') + ggtitle('Top') +
    coord_cartesian(ylim = c(50,250))+
    theme(text = element_text(size = 16), legend.position = 'none')+
    scale_fill_lancet(alpha=0.8) +
    stat_compare_means(method = "wilcox.test",label.x = 0.7,label.y = 50,size = 6))
(p2b=ggboxplot(subset(alpha, depth.watercolumn == 'bottom'), x = 'group', y = 'faith_pd',
               add = "jitter", fill = 'group', size = 1)+
    xlab("Phase") + ylab('Faith PD') + ggtitle('Bottom') +
    coord_cartesian(ylim = c(50,250))+
    theme(text = element_text(size = 16), legend.position = 'none')+
    scale_fill_jama(alpha=0.8) +
    stat_compare_means(method = "wilcox.test",label.x = 0.7,label.y = 50,size = 6))

#evenness matrix
alpha = read_qza("R_analysis_input/evenness_vector.qza")
alpha <- alpha$data %>% tibble::rownames_to_column("SampleID") 
alpha = merge(x = alpha,
              y = metadata.o,
              by.x = 'SampleID',
              by.y = 0,
              all.Y = T)

(p3t=ggboxplot(subset(alpha, depth.watercolumn == 'top'), x = 'group', y = 'pielou_e',
               add = "jitter", fill = 'group', size = 1)+
    xlab("Phase") + ylab('Pielou evenness') + ggtitle('Top') +
    coord_cartesian(ylim = c(0.2,0.8))+
    theme(text = element_text(size = 16), legend.position = 'none')+
    scale_fill_lancet(alpha=0.8) +
    stat_compare_means(method = "wilcox.test",label.x = 0.7,label.y = 0.2,size = 6))

(p3b=ggboxplot(subset(alpha, depth.watercolumn == 'bottom'), x = 'group', y = 'pielou_e',
               add = "jitter", fill = 'group', size = 1)+
    xlab("Phase") + ylab('Pielou evenness') + ggtitle('Bottom') + 
    coord_cartesian(ylim = c(0.2,0.8))+
    theme(text = element_text(size = 16), legend.position = 'none')+
    scale_fill_jama(alpha=0.8) +
    stat_compare_means(method = "wilcox.test",label.x = 0.7,label.y = 0.2,size = 6))

#observed otus
alpha = read_qza("R_analysis_input/observed_otus_vector.qza")
alpha <- alpha$data %>% tibble::rownames_to_column("SampleID") 
alpha = merge(x = alpha,
              y = metadata.o,
              by.x = 'SampleID',
              by.y = 0,
              all.Y = T)

(p4t=ggboxplot(subset(alpha, depth.watercolumn == 'top'), x = 'group', y = 'observed_otus',
               add = "jitter", fill = 'group', size = 1)+
    xlab("Phase") + ylab('Observed ASVs') + ggtitle('Top') + 
    coord_cartesian(ylim = c(0,4000))+
    theme(text = element_text(size = 16), legend.position = 'none')+
    scale_fill_lancet(alpha=0.8) +
    stat_compare_means(method = "wilcox.test",label.x = 0.7,label.y = 0.2,size = 6))

(p4b=ggboxplot(subset(alpha, depth.watercolumn == 'bottom'), x = 'group', y = 'observed_otus',
               add = "jitter", fill = 'group', size = 1)+
    xlab("Phase") + ylab('Observed ASVs') + ggtitle('Bottom') + 
    coord_cartesian(ylim = c(0,4000))+
    theme(text = element_text(size = 16), legend.position = 'none')+
    scale_fill_jama(alpha=0.8) +
    stat_compare_means(method = "wilcox.test",label.x = 0.7,label.y = 0.2,size = 6))

plot_grid(p4t,p2t,p3t,p1t,
          p4b,p2b,p3b,p1b,
          nrow = 2)
ggsave(paste(PFIG,"/",'alpha_diversity.png',sep = ""),width = 18, height = 9, device = 'png')

##############PCoA plot using the Weighted UniFrac distance###############
#read the phylogenetic distance matrix output by QIIME2
dm = read_qza("R_analysis_input/weighted_unifrac_distance_matrix.qza")
#select the rows we need
dm = dist_subset(dm$data,rownames(metadata.o))
# Some distance measures may result in negative eigenvalues. In that case, add a correction:
pcoa.wunifrac = pcoa(dm,correction = 'cailliez')
#transform the pcoa.wunifrac object to a df for making plottable with ggplot2
coordinates.pcoa.wunifrac = data.frame(pcoa.wunifrac$vectors,check.names = FALSE)

###############plot for all samples#############
###############by depths#############
#Adonis test between depths
padonis = adonis_pairwise(metadata.o,dm,'depth.watercolumn',permut = 999)
padonis$Adonis.tab
padonis$Betadisper.tab
#plot
color = metadata.o$depth.watercolumn
(pcoa_all_depth=ggplot()+
    geom_point(data=subset(coordinates.pcoa.wunifrac[,c(1,2)]),
               mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = color),
               alpha = 0.8, size = 3) + 
    labs(title = 'All samples',
         subtitle = paste('PERMANOVA = ', padonis$Adonis.tab[6][[1]], ', ',
                          'PERMDISP = ', padonis$Betadisper.tab[5][[1]], sep = ''))+
    stat_ellipse(data=subset(coordinates.pcoa.wunifrac[,c(1,2)]),
                 mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = color),
                 type = 'norm') +
    xlab(paste('PC1 (',label_percent()(pcoa.wunifrac$values$Rel_corr_eig[1]),')',sep = '')) +
    ylab(paste('PC2 (',label_percent()(pcoa.wunifrac$values$Rel_corr_eig[2]),')',sep = '')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = 'black'),
          legend.title=element_blank(),
          legend.position = 'bottom',
          legend.text = element_text(size=15))+
    scale_color_jama())
ggsave(paste(PFIG,"/",'pcoa_all_depth.png',sep = ""),width = 4, height = 4, device = 'png')

###############by phases#############
padonis = adonis_pairwise(metadata.o,dm,'group',permut = 999)
padonis$Adonis.tab
padonis$Betadisper.tab
#plot
color = metadata.o$group
(pcoa_all_phase=ggplot()+
    geom_point(data=subset(coordinates.pcoa.wunifrac[,c(1,2)]),
               mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = color),
               alpha = 0.8, size = 3) +
    labs(title = 'All samples',
         subtitle = paste('PERMANOVA = ', padonis$Adonis.tab[6][[1]], ', ',
                          'PERMDISP = ', padonis$Betadisper.tab[5][[1]], sep = ''))+
    stat_ellipse(data=subset(coordinates.pcoa.wunifrac[,c(1,2)]),
                 mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = color),
                 type = 'norm') +
    xlab(paste('PC1 (',label_percent()(pcoa.wunifrac$values$Rel_corr_eig[1]),')',sep = '')) +
    ylab(paste('PC2 (',label_percent()(pcoa.wunifrac$values$Rel_corr_eig[2]),')',sep = '')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = 'black'),
          legend.title=element_blank(),
          legend.position = 'bottom',
          legend.text = element_text(size=15))+
    scale_color_lancet())
ggsave(paste(PFIG,"/",'pcoa_all_phase.png',sep = ""),width = 4, height = 4, device = 'png')

###############plot for top samples by phase#############
#subset
metadata.o.s = subset(metadata.o, depth.watercolumn == 'top')
#repeat the pcoa step
dm = read_qza("R_analysis_input/weighted_unifrac_distance_matrix.qza")
dm = dist_subset(dm$data,rownames(metadata.o.s))
pcoa.wunifrac = pcoa(dm,correction = 'cailliez')
coordinates.pcoa.wunifrac = data.frame(pcoa.wunifrac$vectors,check.names = FALSE)

#differences between phases
padonis = adonis_pairwise(metadata.o.s,dm,'group',permut = 999)
padonis$Adonis.tab
padonis$Betadisper.tab
#plot
color = metadata.o.s$group
(pcoa_top=ggplot()+
    geom_point(data=subset(coordinates.pcoa.wunifrac[,c(1,2)]),
               mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = color),
               alpha = 0.8, size = 3) +
    labs(title = 'Top',
         subtitle = paste('PERMANOVA = ', padonis$Adonis.tab[6][[1]], ', ',
                          'PERMDISP = ', padonis$Betadisper.tab[5][[1]], sep = ''))+
    stat_ellipse(data=subset(coordinates.pcoa.wunifrac[,c(1,2)]),
                 mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = color),
                 type = 'norm') +
    xlab(paste('PC1 (',label_percent()(pcoa.wunifrac$values$Rel_corr_eig[1]),')',sep = '')) +
    ylab(paste('PC2 (',label_percent()(pcoa.wunifrac$values$Rel_corr_eig[2]),')',sep = '')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = 'black'),
          legend.title=element_blank(),
          legend.position = 'bottom',
          legend.text = element_text(size=15))+
    scale_color_lancet())
ggsave(paste(PFIG,"/",'pcoa_top_phase.png',sep = ""),width = 4, height = 4, device = 'png')

###############plot for bottom samples by phases#############
#subset
metadata.o.s = subset(metadata.o, depth.watercolumn == 'bottom')
#repeat the pcoa step
dm = read_qza("R_analysis_input/weighted_unifrac_distance_matrix.qza")
dm = dist_subset(dm$data,rownames(metadata.o.s))
pcoa.wunifrac = pcoa(dm,correction = 'cailliez')
coordinates.pcoa.wunifrac = data.frame(pcoa.wunifrac$vectors,check.names = FALSE)

#differences between phases
padonis = adonis_pairwise(metadata.o.s,dm,'group',permut = 999)
padonis$Adonis.tab
padonis$Betadisper.tab
#plot
color = metadata.o.s$group
(pcoa_bottom=ggplot()+
    geom_point(data=subset(coordinates.pcoa.wunifrac[,c(1,2)]),
               mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = color),
               alpha = 0.8, size = 3) +
    labs(title = 'Bottom',
         subtitle = paste('PERMANOVA = ', padonis$Adonis.tab[6][[1]], ', ',
                          'PERMDISP = ', padonis$Betadisper.tab[5][[1]], sep = ''))+
    stat_ellipse(data=subset(coordinates.pcoa.wunifrac[,c(1,2)]),
                 mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = color),
                 type = 'norm') +
    xlab(paste('PC1 (',label_percent()(pcoa.wunifrac$values$Rel_corr_eig[1]),')',sep = '')) +
    ylab(paste('PC2 (',label_percent()(pcoa.wunifrac$values$Rel_corr_eig[2]),')',sep = '')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = 'black'),
          legend.title=element_blank(),
          legend.position = 'bottom',
          legend.text = element_text(size=15))+
    scale_color_lancet())
ggsave(paste(PFIG,"/",'pcoa_bottom_phase.png',sep = ""),width = 4, height = 4, device = 'png')

##########################Differential abundance analysis###########################
#use the un-rarefied table for ANCOMBC analysis as suggested by the authors of ANCOM-BC
pu<-qza_to_phyloseq(
  features="R_analysis_input/unrarefied_table.qza",
  tree="R_analysis_input/rooted-tree.qza",
  taxonomy="R_analysis_input/taxonomy.qza")
pu.asv = data.frame(otu_table(pu),check.names = F)
pu.tax = data.frame(tax_table(pu),check.names = F)
pu.tre = phy_tree(pu)

#remove discarded samples
pu.asv = pu.asv[,names(pu.asv) %in% rownames(metadata.o)]
#remove the ghost ASVs
pu.asv = subset(pu.asv, rowSums(pu.asv)>0)

#check if sample id matches between metadata and the ASV table
pu.sampleid.check = match.name(cn.list=list(comm=pu.asv), rn.list = list(env=metadata.o))
#check if the ASV id matches across the ASV table, phylogenetic tree tips and taxonomic table
pu.sampleid.check = match.name(rn.list = list(asv=pu.asv,tax=pu.tax),tree.list=list(tree=pu.tre)) # ghost ASVs have been removed from the tax tabla and tree tips
pu.asv <- pu.sampleid.check$asv
pu.tax <- pu.sampleid.check$tax
pu.tre <- pu.sampleid.check$tree
dim(pu.asv) #22717,84
dim(pu.tax) #22717,7
dim(metadata.o) #84,13

#reconstruct the phyloseq object pu
pu <- NA
pu = phyloseq(otu_table(as.matrix(pu.asv),taxa_are_rows = T),
              tax_table(as.matrix(pu.tax)),
              phy_tree(pu.tre),
              sample_data(metadata.o))

#define threshold to subsample taxa, e.g. if threshold = 0.1 
#then only the taxa > 0.01% relative abundance (total counts across all samples) will be plotted
threshold = 0.0001

####Genus level comparison######
#agglomerating rank (all ASV with phylum = NA will be removed)
pu.ph = tax_glom(pu, taxrank = 'Genus', NArm = F)
summary(rownames(pu.ph@otu_table) ==rownames(pu.ph@tax_table)) #sanity check

#subset to bottom
pu.ph.s = subset_samples(pu.ph, depth.watercolumn %in% c('bottom'))
pu.ph.s.ac = ancombc(phyloseq = pu.ph.s, formula = 'group',
                     p_adj_method = 'holm', zero_cut = 0.90, lib_cut = 0,
                     group = 'group', struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                     max_iter = 100, conserve = TRUE, alpha = 0.05, global = F)

#coefficient (beta) is the log-fold change in comparison with the reference group
res = pu.ph.s.ac$res

#gather beta (log change)
res.beta = res$beta
res.beta$ASV = rownames(res.beta)
res.beta = gather(res.beta,"comparison","beta",-ASV)
#gather se
res.se = res$se
res.se$ASV = rownames(res.se)
res.se = gather(res.se, "comparison", "se", -ASV)
#gather adjusted p values
res.qval = res$q_val
res.qval$ASV = rownames(res.qval)
res.qval = gather(res.qval,"comparison","qval",-ASV)

#sanity check
summary(res.beta$ASV == res.se$ASV & res.beta$ASV == res.qval$ASV)
summary(res.beta$comparison == res.se$comparison & res.beta$comparison == res.qval$comparison)

pu.ph.b = data.frame(res.beta, 
                     se = res.se$se,
                     qval = res.qval$qval,
                     sample = 'bottom') 

#rename ASV to taxonomy
tax = as.data.frame(tax_table(pu.ph.s))
tax$concatenated_tax = paste(tax$Family,tax$Genus, sep = '_')
pu.ph.b$tax = tax$concatenated_tax[match(pu.ph.b$ASV,rownames(tax))]
anc.res.ge.b <- pu.ph.b

#subset to ASVs with relative abundance > 0.01% 
pu.ph.s = prune_taxa(taxa_sums(pu.ph.s)>sum(taxa_sums(pu.ph.s))*threshold, pu.ph.s)
(pu.ph.b = subset(pu.ph.b, ASV %in% rownames(otu_table(pu.ph.s))))

#remove the unclassified ASVs
pu.ph.b = pu.ph.b[complete.cases(pu.ph.b),]
pu.ph.b = subset(pu.ph.b, grepl('NA',pu.ph.b$tax) == F)

#select taxa with q <= 0.05
pu.ph.b = subset(pu.ph.b, qval <= 0.05)
#change adjusted q value to stars
pu.ph.b$qval = ifelse(pu.ph.b$qval > 0.01 & pu.ph.b$qval <= 0.05, "*", 
                      ifelse(pu.ph.b$qval > 0.001 & pu.ph.b$qval <= 0.01, '**',
                             ifelse(pu.ph.b$qval <= 0.001, '***', '')))

#start plotting
ggplot(data = pu.ph.b, 
       mapping = aes(x = tax,y = beta,fill = comparison))+
  geom_stripped_cols(odd = "ivory2", even = 'white')+
  geom_bar(stat = 'identity',
           position=position_dodge())+
  geom_errorbar(mapping = aes(ymin = beta - se,
                              ymax = beta + se),
                width = 0.5,
                size = 0.15,
                position = position_dodge(width = 0.9)) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5)+
  facet_grid(.~sample, scales = "free_x")+
  labs(x=NULL, y="log(off or low/on and high)")+
  coord_flip()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        text = element_text(size = 16, colour = 'black'),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"),
        strip.text.y = element_text(size = 15, color = "black", face = "bold")) +
  geom_text(mapping = aes(x = tax, y = beta + sign(beta), label = qval),
            position = position_dodge(width = 1),
            size = 3, color = 'red')+
  scale_fill_uchicago(alpha=0.8)+
  ylim(-4,4)
ggsave(paste(PFIG,'ancombc_genus_b.png',sep = '/'),width = 8,height = 3, device = 'png')

#subset "top" samples
pu.ph.s = subset_samples(pu.ph, depth.watercolumn %in% c('top')) #change here
pu.ph.s.ac = ancombc(phyloseq = pu.ph.s, formula = 'group',
                     p_adj_method = 'holm', zero_cut = 0.90, lib_cut = 0,
                     group = 'group', struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                     max_iter = 100, conserve = TRUE, alpha = 0.05, global = F)

#coefficient (beta) is the log-fold change in comparison with the reference group
res = pu.ph.s.ac$res

#gather beta (log change)
res.beta = res$beta
res.beta$ASV = rownames(res.beta)
res.beta = gather(res.beta,"comparison","beta",-ASV)
#gather se
res.se = res$se
res.se$ASV = rownames(res.se)
res.se = gather(res.se, "comparison", "se", -ASV)
#gather adjusted p value
res.qval = res$q_val
res.qval$ASV = rownames(res.qval)
res.qval = gather(res.qval,"comparison","qval",-ASV)

#sanity check
summary(res.beta$ASV == res.se$ASV & res.beta$ASV == res.qval$ASV)
summary(res.beta$comparison == res.se$comparison & res.beta$comparison == res.qval$comparison)

pu.ph.b = data.frame(res.beta,
                     se = res.se$se,
                     qval = res.qval$qval,
                     sample = 'top') 

#rename ASV to taxonomy
tax = as.data.frame(tax_table(pu.ph.s))
tax$concatenated_tax = paste(tax$Family,tax$Genus, sep = '_')
pu.ph.b$tax = tax$concatenated_tax[match(pu.ph.b$ASV,rownames(tax))]
anc.res.ge.t <- pu.ph.b
pu.ph.s = prune_taxa(taxa_sums(pu.ph.s)>sum(taxa_sums(pu.ph.s))*threshold, pu.ph.s) 
(pu.ph.b = subset(pu.ph.b, ASV %in% rownames(otu_table(pu.ph.s))))

#remove the unclassified ASVs
pu.ph.b = pu.ph.b[complete.cases(pu.ph.b),]
pu.ph.b = subset(pu.ph.b, grepl('NA',pu.ph.b$tax) == F)

#select taxa with q <= 0.05
pu.ph.b = subset(pu.ph.b, qval <= 0.05)
#change q value to stars
pu.ph.b$qval = ifelse(pu.ph.b$qval > 0.01 & pu.ph.b$qval <= 0.05, "*", 
                      ifelse(pu.ph.b$qval > 0.001 & pu.ph.b$qval <= 0.01, '**',
                             ifelse(pu.ph.b$qval <= 0.001, '***', '')))

#start plotting
ggplot(data = pu.ph.b, 
       mapping = aes(x = tax,y = beta,fill = comparison))+
  geom_stripped_cols(odd = "ivory2", even = 'white')+
  geom_bar(stat = 'identity',
           position=position_dodge())+
  geom_errorbar(mapping = aes(ymin = beta - se,
                              ymax = beta + se),
                width = 0.5,
                size = 0.15,
                position = position_dodge(width = 0.9)) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5)+
  facet_grid(.~sample, scales = "free_x")+
  labs(x=NULL, y="log(off or low/on and high)")+
  coord_flip()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        text = element_text(size=16),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"),
        strip.text.y = element_text(size = 15, color = "black", face = "bold"))+
  geom_text(mapping = aes(x = tax, y = beta + sign(beta), label = qval),
            position = position_dodge(width = 1),
            size = 3, color = 'red')+
  scale_fill_uchicago(alpha=0.8)+
  ylim(-4,4)
ggsave(paste(PFIG,'ancombc_genus_t.png',sep = '/'),width = 8,height = 1.2, device = 'png')

#produce an overarching table to document the ancombc results in taxa comparison
anc.res.taxa = bind_rows(lst(anc.res.ge.b, anc.res.ge.t), .id = 'id')
write.csv(anc.res.taxa, paste(FAS,'/ancombc_tax_results_genus.csv',sep=''))

#####################ASV level comparison########################
pu.ph <- pu
summary(rownames(pu.ph@otu_table) ==rownames(pu.ph@tax_table)) 

#subset "top" samples
pu.ph.s = subset_samples(pu.ph, depth.watercolumn %in% c('top')) 
pu.ph.s.ac = ancombc(phyloseq = pu.ph.s, formula = 'group',
                     p_adj_method = 'holm', zero_cut = 0.90, lib_cut = 0,
                     group = 'group', struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                     max_iter = 100, conserve = TRUE, alpha = 0.05, global = F)

#coefficient (beta) is the log-fold change in comparison with the reference group
res = pu.ph.s.ac$res

#gather beta (log change)
res.beta = res$beta
res.beta$ASV = rownames(res.beta)
res.beta = gather(res.beta,"comparison","beta",-ASV)
#gather se
res.se = res$se
res.se$ASV = rownames(res.se)
res.se = gather(res.se, "comparison", "se", -ASV)
#gather adjusted p value
res.qval = res$q_val
res.qval$ASV = rownames(res.qval)
res.qval = gather(res.qval,"comparison","qval",-ASV)

#sanity check
summary(res.beta$ASV == res.se$ASV & res.beta$ASV == res.qval$ASV)
summary(res.beta$comparison == res.se$comparison & res.beta$comparison == res.qval$comparison)

pu.ph.b = data.frame(res.beta,
                     se = res.se$se,
                     qval = res.qval$qval,
                     sample = 'top')

#rename ASV to taxonomy
tax = as.data.frame(tax_table(pu.ph.s))
tax$concatenated_tax = paste(tax$Phylum,tax$Class,tax$Order,tax$Family,tax$Genus, sep = '_')
pu.ph.b$tax = tax$concatenated_tax[match(pu.ph.b$ASV,rownames(tax))]
anc.res.asv.t <- pu.ph.b

#subset "bottom" samples
pu.ph.s = subset_samples(pu.ph, depth.watercolumn %in% c('bottom'))
pu.ph.s.ac = ancombc(phyloseq = pu.ph.s, formula = 'group',
                     p_adj_method = 'holm', zero_cut = 0.90, lib_cut = 0,
                     group = 'group', struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                     max_iter = 100, conserve = TRUE, alpha = 0.05, global = F)

#coefficient (beta) is the log-fold change in comparison with the reference group
res = pu.ph.s.ac$res

#gather beta (log change)
res.beta = res$beta
res.beta$ASV = rownames(res.beta)
res.beta = gather(res.beta,"comparison","beta",-ASV)
#gather se
res.se = res$se
res.se$ASV = rownames(res.se)
res.se = gather(res.se, "comparison", "se", -ASV)
#gather adjusted p value
res.qval = res$q_val
res.qval$ASV = rownames(res.qval)
res.qval = gather(res.qval,"comparison","qval",-ASV)

#sanity check
summary(res.beta$ASV == res.se$ASV & res.beta$ASV == res.qval$ASV)
summary(res.beta$comparison == res.se$comparison & res.beta$comparison == res.qval$comparison)

pu.ph.b = data.frame(res.beta, 
                     se = res.se$se,
                     qval = res.qval$qval,
                     sample = 'bottom') 

#rename ASV to taxonomy
tax = as.data.frame(tax_table(pu.ph.s))
tax$concatenated_tax = paste(tax$Phylum,tax$Class,tax$Order,tax$Family,tax$Genus, sep = '_')
pu.ph.b$tax = tax$concatenated_tax[match(pu.ph.b$ASV,rownames(tax))]
anc.res.asv.b <- pu.ph.b
anc.res.asv = bind_rows(list(anc.res.asv.b,anc.res.asv.t), .id = 'id')
write.csv(anc.res.asv, paste(FAS,'/ancombc_tax_results_ASV.csv',sep=''))

#subset stations RC1 and 2 and repeat the ASV level comparison
pu.ph <- pu
#subset RC1 and 2
pu.ph = subset_samples(pu.ph, station %in% c('RC1','RC2')) 

#subset "top" samples
pu.ph.s = subset_samples(pu.ph, depth.watercolumn %in% c('top'))
pu.ph.s.ac = ancombc(phyloseq = pu.ph.s, formula = 'group',
                     p_adj_method = 'holm', zero_cut = 0.90, lib_cut = 0,
                     group = 'group', struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                     max_iter = 100, conserve = TRUE, alpha = 0.05, global = F)

#coefficient (beta) is the log-fold change in comparison with the reference group
res = pu.ph.s.ac$res

#gather beta (log change)
res.beta = res$beta
res.beta$ASV = rownames(res.beta)
res.beta = gather(res.beta,"comparison","beta",-ASV)
#gather se
res.se = res$se
res.se$ASV = rownames(res.se)
res.se = gather(res.se, "comparison", "se", -ASV)
#gather adjusted p value
res.qval = res$q_val
res.qval$ASV = rownames(res.qval)
res.qval = gather(res.qval,"comparison","qval",-ASV)

#sanity check
summary(res.beta$ASV == res.se$ASV & res.beta$ASV == res.qval$ASV)
summary(res.beta$comparison == res.se$comparison & res.beta$comparison == res.qval$comparison)

pu.ph.b = data.frame(res.beta, 
                     se = res.se$se,
                     qval = res.qval$qval,
                     sample = 'top') 

#rename ASV to taxonomy
tax = as.data.frame(tax_table(pu.ph.s))
tax$concatenated_tax = paste(tax$Phylum,tax$Class,tax$Order,tax$Family,tax$Genus, sep = '_')
pu.ph.b$tax = tax$concatenated_tax[match(pu.ph.b$ASV,rownames(tax))]
anc.res.asv.t <- pu.ph.b

#subset "bottom" samples
pu.ph.s = subset_samples(pu.ph, depth.watercolumn %in% c('bottom')) 
pu.ph.s.ac = ancombc(phyloseq = pu.ph.s, formula = 'group',
                     p_adj_method = 'holm', zero_cut = 0.90, lib_cut = 0,
                     group = 'group', struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                     max_iter = 100, conserve = TRUE, alpha = 0.05, global = F)

#coefficient (beta) is the log-fold change in comparison with the reference group
res = pu.ph.s.ac$res

#gather beta (log change)
res.beta = res$beta
res.beta$ASV = rownames(res.beta)
res.beta = gather(res.beta,"comparison","beta",-ASV)
#gather se
res.se = res$se
res.se$ASV = rownames(res.se)
res.se = gather(res.se, "comparison", "se", -ASV)
#gather adjusted p value
res.qval = res$q_val
res.qval$ASV = rownames(res.qval)
res.qval = gather(res.qval,"comparison","qval",-ASV)

#sanity check
summary(res.beta$ASV == res.se$ASV & res.beta$ASV == res.qval$ASV)
summary(res.beta$comparison == res.se$comparison & res.beta$comparison == res.qval$comparison)

pu.ph.b = data.frame(res.beta, 
                     se = res.se$se,
                     qval = res.qval$qval,
                     sample = 'bottom') 

#rename ASV to taxonomy
tax = as.data.frame(tax_table(pu.ph.s))
tax$concatenated_tax = paste(tax$Phylum,tax$Class,tax$Order,tax$Family,tax$Genus, sep = '_')
pu.ph.b$tax = tax$concatenated_tax[match(pu.ph.b$ASV,rownames(tax))]
anc.res.asv.b <- pu.ph.b
anc.res.asv = bind_rows(list(anc.res.asv.b,anc.res.asv.t), .id = 'id')
#write results to files
write.csv(anc.res.asv, paste(FAS,'/ancombc_tax_results_ASV_RC1-2_subset.csv',sep=''))

#########################PICRUSt analysis###########################
############prepare the input ASV table##########
#renames Sample ID to meaningful sample names
(names(pu.asv) <- SIDconversion$SD[match(names(pu.asv),SIDconversion$SIDmetadata)])
pu.asv$ASVID <- rownames(pu.asv)
pu.asv = pu.asv[,c(ncol(pu.asv),1:ncol(pu.asv)-1)]
write.table(x = pu.asv, file = paste(FAS,"/unrarefied_asv_table_new_group.txt", sep = ''),
            row.names = F, col.names = T, sep = '\t',quote = F)

###############data analysis####################
#read in the unstratified picrust KO output
pcrst = read.csv('R_analysis_input/PICRUSt2_unstrat.tsv', sep = '\t',
                 check.names = F)
#set rownames
row.names(pcrst) = pcrst[,1]
pcrst = pcrst[,-1]
#explore the pcrst results
summary(pcrst)
str(pcrst)
summary(colSums(pcrst[,-1]))

#read customized kegg module table
km = read.csv('R_analysis_input/kegg_key_metabolic_processes.csv', check.names = F)
#parse km
km$KO_description = paste(km$KO,km$KO_description,sep = ' ')
#remove EC to shorten the KO_description
km$KO_description = gsub(";.*", "", km$KO_description)
#remove the duplicate genes (some genes that appear in more than one pathways)
km.short = distinct(subset(km, select = c(KO, KO_description)))

#use ANCOMBC for analysis
#create a phyloseq object
metadata.or$group = factor(metadata.or$group, levels = c('on and high', 'off or low'))
metadata.or$depth.watercolumn = factor(metadata.or$depth.watercolumn, levels = c('top','bottom'))
pk = phyloseq(otu_table(as.matrix(pcrst),taxa_are_rows = T),sample_data(metadata.or))

###########bottom sample analysis####################
pk.s = subset_samples(pk, depth.watercolumn %in% c('bottom')) 
pk.s.ac = ancombc(phyloseq = pk.s, formula = 'group',
                  p_adj_method = 'holm', zero_cut = 0.90, lib_cut = 0,
                  group = 'group', struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                  max_iter = 100, conserve = TRUE, alpha = 0.05, global = F)

#coefficient (beta) is the log-fold change in comparison with the reference group
res = pk.s.ac$res

#gather beta (log change)
res.beta = res$beta
res.beta$ASV = rownames(res.beta)
res.beta = gather(res.beta,"comparison","beta",-ASV)
#gather se
res.se = res$se
res.se$ASV = rownames(res.se)
res.se = gather(res.se, "comparison", "se", -ASV)
#gather adjusted p value
res.qval = res$q_val
res.qval$ASV = rownames(res.qval)
res.qval = gather(res.qval,"comparison","qval",-ASV)

#sanity check
summary(res.beta$ASV == res.se$ASV & res.beta$ASV == res.qval$ASV)
summary(res.beta$comparison == res.se$comparison & res.beta$comparison == res.qval$comparison)

pk.anc.s = data.frame(res.beta, 
                      se = res.se$se,
                      qval = res.qval$qval,
                      sample = 'bottom') 

pk.anc.s= merge(x = pk.anc.s,
                y = km.short,
                by.x = 'ASV',
                by.y = 'KO',
                all.y = T)
pk.anc.s = pk.anc.s[complete.cases(pk.anc.s),]
pk.anc.b <- pk.anc.s

###########top sample analysis####################
pk.s = subset_samples(pk, depth.watercolumn %in% c('top')) 
pk.s.ac = ancombc(phyloseq = pk.s, formula = 'group',
                  p_adj_method = 'holm', zero_cut = 0.90, lib_cut = 0,
                  group = 'group', struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                  max_iter = 100, conserve = TRUE, alpha = 0.05, global = F)

#coefficient (beta) is the log-fold change in comparison with the reference group
res = pk.s.ac$res

#gather beta (log change)
res.beta = res$beta
res.beta$ASV = rownames(res.beta)
res.beta = gather(res.beta,"comparison","beta",-ASV)
#gather se
res.se = res$se
res.se$ASV = rownames(res.se)
res.se = gather(res.se, "comparison", "se", -ASV)
#gather adjusted p value
res.qval = res$q_val
res.qval$ASV = rownames(res.qval)
res.qval = gather(res.qval,"comparison","qval",-ASV)

#sanity check
summary(res.beta$ASV == res.se$ASV & res.beta$ASV == res.qval$ASV)
summary(res.beta$comparison == res.se$comparison & res.beta$comparison == res.qval$comparison)

pk.anc.s = data.frame(res.beta, 
                      se = res.se$se,
                      qval = res.qval$qval,
                      sample = 'top') 

pk.anc.s= merge(x = pk.anc.s,
                y = km.short,
                by.x = 'ASV',
                by.y = 'KO',
                all.y = T)
pk.anc.s = pk.anc.s[complete.cases(pk.anc.s),]
pk.anc.t <- pk.anc.s

#visualize the results 
#select KOs that are at least significant in one comparison
pk.anc.selectedKO = unique(c(subset(pk.anc.t, qval <= 0.05)$ASV, subset(pk.anc.b, qval <= 0.05)$ASV))
#sanity check
summary(subset(pk.anc.t, qval <= 0.05)$ASV %in% pk.anc.selectedKO)
summary(subset(pk.anc.b, qval <= 0.05)$ASV %in% pk.anc.selectedKO)
#combine the results into a table
pk.anc.all = bind_rows(lst(subset(pk.anc.t, ASV %in% pk.anc.selectedKO),
                           subset(pk.anc.b, ASV %in% pk.anc.selectedKO)))
#order the KOs based on the pathway it belongs to
km.order = subset(km, select = c(KO_description,metabolism), km$KO %in% pk.anc.selectedKO)
km.order = distinct(km.order)

pk.anc.all$KO_description = factor(pk.anc.all$KO_description,
                                   levels = km.order$KO_description)

#manually add the fwdH into the df. It is missing here because no fwdh gene is inferred from the top samples
pk.anc.all = rbind(pk.anc.all,c('K00204','groupoff/low',0,0,1,'top','K00204 fwdH'))
pk.anc.all$beta = as.numeric(pk.anc.all$beta)
pk.anc.all$se = as.numeric(pk.anc.all$se)
#change qvalue to stars
pk.anc.all$qval = ifelse(pk.anc.all$qval > 0.01 & pk.anc.all$qval <= 0.05, "*", 
                         ifelse(pk.anc.all$qval > 0.001 & pk.anc.all$qval <= 0.01, '**',
                                ifelse(pk.anc.all$qval <= 0.001, '***', '')))

#plot
pk.anc.all %>% mutate(dpt = ifelse(sample == 'top','a. top', 'b. bottom')) -> pk.anc.all
pk.anc.all$dpt = factor(pk.anc.all$dpt,levels = c('a. top','b. bottom'))

ggplot(data = pk.anc.all, 
       mapping = aes(x = KO_description,y = beta,fill = comparison))+
  geom_stripped_cols(odd = "ivory2", even = 'white')+
  geom_bar(stat = 'identity',
           position=position_dodge())+
  geom_errorbar(mapping = aes(ymin = beta - se,
                              ymax = beta + se),
                width = 0.5,
                size = 0.15,
                position = position_dodge(width = 0.9)) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5)+
  facet_grid(.~dpt, scales = "free_x")+
  labs(x=NULL, y="log(off or low/on and high)")+
  coord_flip()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="white"),
        legend.position = 'none',
        strip.text.x = element_text(size = 10, color = "black", face = "bold"),
        strip.text.y = element_text(size = 10, color = "black", face = "bold"),)+
  geom_text(mapping = aes(x = KO_description, y = beta + sign(beta), label = qval),
            position = position_dodge(width = 1),
            size = 3, color = 'red')+
  scale_fill_lancet(alpha=0.8)+
  ylim(-2,2)
ggsave(paste(PFIG,"/",'ancombc_picrust.png',sep = ""),width = 5, height = 5, device = 'png')
write.csv(pk.anc.all,paste(FAS,'ancombc_picrust.csv',sep='/'))

