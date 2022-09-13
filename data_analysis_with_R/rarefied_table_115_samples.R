#Creator: Yue Zhang
#This scripts is used to 
# 1) analyze the subset of 115 samples by grouping them into 3 aerator states -- aerators "on", "off" and "altered".
# 2) analyze the control samples

###########Preparation################
#load required packages
packages_to_load = c('stringr','phyloseq','vegan','tidyr','ggplot2','qiime2R',
                     'dplyr','ANCOMBC','usedist','rioja','ggpubr','spaa','ape',
                     'cowplot','ggforestplot','GGally','scales','metagMisc','ggsci','rstatix','wesanderson')
lapply(packages_to_load,require, character.only = T)

#generating color palettes
col1 = wes_palette("Rushmore1")
col2 = wes_palette("GrandBudapest2")
col3 = wes_palette('Darjeeling2')
col21 <- rev(c("tomato1","darkblue","turquoise1","lightblue","darkred","mediumblue","purple","bisque",
               "greenyellow","yellow","violetred2","darkgreen","darkgoldenrod1","deeppink3","cadetblue4",
               "orchid2","seagreen3","purple4","dodgerblue2","red","gray27"))

###########load variables##########
PFIG = FAS = 'R_analysis_output/115_samples'
#dir.create(PFIG)

#set random seed
set.seed(666)

############metadata pre-processing###################
metadata = read.csv('R_analysis_input/metadata_simplified_115.csv')
#keep the column type row
metadata_typerow = metadata[1,]
metadata = metadata[-1,]
#create a new identifier column in the metadata
metadata$date = as.character(metadata$date)
metadata$date[metadata$date == '8'] = '08' #this is to avoid the issues when ordering them after converting them into characters in the next step
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
metadata %>% mutate(group = ifelse(date_numeric<=9 | date_numeric >= 29, 'on',
                                   ifelse(date_numeric >= 10 & date_numeric <= 22, 'off', 'altered'))) -> metadata

#combine with the column type row
metadata = bind_rows(metadata,metadata_typerow)
metadata$group[nrow(metadata)] = 'categorical'
metadata$date_numeric[nrow(metadata)] = 'numeric'
metadata$date[nrow(metadata)] = 'categorical'
metadata$date_time[nrow(metadata)] = 'categorical'
metadata$depth.station.date.time[nrow(metadata)] = 'categorical'
metadata = metadata[c(nrow(metadata),1:nrow(metadata)-1),]

write.table(metadata, paste(FAS,'metadata_115_preprocessed.csv',sep = '/'),sep = '\t', quote = F,row.names = F)

############input data using phyloseq package###############
#Create a phyloseq object
pr<-qza_to_phyloseq(
  features="R_analysis_input/rarefied_table.qza",
  tree="R_analysis_input/rooted-tree.qza",
  taxonomy="R_analysis_input/taxonomy.qza",
  metadata=paste(FAS,'metadata_115_preprocessed.csv',sep = '/'))
pr
asv = data.frame(otu_table(pr)) #note that the asv is the raw table which has not been normalized by either rows or columns
metadata= data.frame(sample_data(pr))
tax = data.frame(tax_table(pr))
metadata$group = factor(metadata$group, levels = c('on','altered','off'))

#####################ASV table pre-processing##############################
#order the asv table based on the total counts across all samples
asv$sum = apply(asv,1,sum)
asv.o = asv[order(asv$sum,decreasing = T),]
head(asv.o)

#remove the ghost ASVs
asv.o %>% dplyr::filter(sum >0) %>% as.data.frame -> asv.of 
nrow(asv.of) #the total number of ASVs after rarefaction = 24459
sum(asv.of) #the total number of reads after rarefaction = 10906600

#simplify the ASV ID
nums <- sprintf('%0.5d', 1:nrow(asv.of)) ##creating string "00001", "00002"。。。
ASVnames <- paste("ASV",nums, sep="")
ASVconversion <- data.frame(ASV_ID = rownames(asv.of),ASVnames, stringsAsFactors = F) #make df with two columns
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

#total number of ASVs AND reads of the rarefied table
sum(asv.of)
nrow(asv.of)
#TOTAL NUMBER OF ASV = 24459
#TOTAL NUMBER OF READS = 10906600

#calculate frequency -- normalize by sample
asv.off = decostand(asv.of,2,method = 'total')
colSums(asv.off) #col sum = 1

#############taxonomic table pre-processing################
#remove the ghost ASVs in the taxonomy table
tax = subset(tax, rownames(tax) %in% ASVconversion$ASV_ID)
rownames(tax) = ASVconversion$ASVnames[match(rownames(tax),ASVconversion$ASV_ID)]
tax.o = tax[order(rownames(tax)),]

######################Taxonomic bar plot#####################
#######################Phylum level bar plot###############
#aggregate the ASV table to the phylum level
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

#retain the 9 most abundant Phyla; calculate the sum of all the others
ft = 9 #features to retain

asv.offp = asv.offp[order(rowSums(asv.offp), decreasing = T),]  
asv.offps = bind_rows(asv.offp[c(1:ft),],colSums(asv.offp[c(c(ft+1):nrow(asv.offp)),]))
rownames(asv.offps)[c(ft+1)] = 'All others'
summary(colSums(asv.offps))

#stack asv.offps
rownames(asv.offps) -> asv.offps$Phylum
asv.offps.s = gather(asv.offps,"SD","frequency",-Phylum)

#merge the stacked ASV table with a subset of metadata to include some categorization columns
metadata.o %>%
  dplyr::select(`depth.station.date.time`,`depth.watercolumn`,`group`,`date_time`,`station`) %>%
  merge(x = ., y = asv.offps.s, all.x = T, by.x = "depth.station.date.time", by.y = "SD") -> asv.offps.s

#coerce the phylum column to ordered factors for plotting
asv.offps.s$Phylum = factor(asv.offps.s$Phylum, levels = c(rownames(asv.offps)))

#make bar plot
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

#####################ASV bar plot for cyanobacteria########################
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

#stack the ordered ASV table
tax.om.cyaso$Class <- rownames(tax.om.cyaso)
tax.om.cyaso.s = gather(tax.om.cyaso,"SD","frequency",-Class)

#merge the stacked table with a subset of metadata to include some categorization columns
metadata.o %>%
  dplyr::select(`depth.station.date.time`,`depth.watercolumn`,`group`,`date_time`,`station`) %>%
  merge(x = ., y = tax.om.cyaso.s, all.x = T, by.x = "depth.station.date.time", by.y = "SD") -> tax.om.cyaso.s

#plotting
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

#######################alpha diversity analysis############################
#shannon diversity matrix
alpha = read_qza("R_analysis_input/shannon_vector.qza")
alpha <- alpha$data %>% tibble::rownames_to_column("SampleID") 
alpha = merge(x = alpha,
              y = metadata.o,
              by.x = 'SampleID',
              by.y = 0,
              all.Y = T)

#pairwise wilcox test with adjusting p values
stat.test <- subset(alpha, depth.watercolumn == 'top') %>%
  wilcox_test(shannon ~ group) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "group", dodge = 0.8)

(p1t=ggboxplot(subset(alpha, depth.watercolumn == 'top'), x = 'group', y = 'shannon',
               add = "jitter", fill = 'group', size = 1)+
    xlab("Phase") + ylab('Shannon index') + ggtitle('Top') +
    coord_cartesian(ylim = c(2,12))+
    theme(text = element_text(size = 22), legend.position = 'none')+
    scale_fill_lancet(alpha=0.8) +
    stat_compare_means(method = "anova",label.x = 0.7,label.y = 11.5,size = 6)+
    stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0.05, step.increase = 0.15, size = 5)) # Add the pairwise comparison p-values

stat.test <- subset(alpha, depth.watercolumn == 'bottom') %>%
  wilcox_test(shannon ~ group) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "group", dodge = 0.8)

(p1b=ggboxplot(subset(alpha, depth.watercolumn == 'bottom'), x = 'group', y = 'shannon',
               add = "jitter", fill = 'group', size = 1)+
    xlab("Phase") + ylab('Shannon index') + ggtitle('Bottom')+
    coord_cartesian(ylim = c(2,12))+
    theme(text = element_text(size = 22), legend.position = 'none')+
    scale_fill_jama(alpha=0.8) +
    stat_compare_means(method = "anova",label.x = 0.7,label.y = 2,size = 6)+
    stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0.05, step.increase = 0.15, size = 5)) # Add pairwise comparisons p-value

#Faith PD matrix
alpha = read_qza("R_analysis_input/faith_pd_vector.qza")
alpha <- alpha$data %>% tibble::rownames_to_column("SampleID") 
alpha = merge(x = alpha,
              y = metadata.o,
              by.x = 'SampleID',
              by.y = 0,
              all.Y = T)

stat.test <- subset(alpha, depth.watercolumn == 'top') %>%
  wilcox_test(faith_pd ~ group) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "group", dodge = 0.8)

(p2t=ggboxplot(subset(alpha, depth.watercolumn == 'top'), x = 'group', y = 'faith_pd',
               add = "jitter", fill = 'group', size = 1)+
    xlab("Phase") + ylab('Faith PD') + ggtitle('Top') +
    coord_cartesian(ylim = c(50,350))+
    scale_fill_lancet(alpha=0.8) +
    theme(text = element_text(size = 22), legend.position = 'none')+
    stat_compare_means(method = "anova",label.x = 0.7,label.y = 320,size = 6)+
    stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0.05, step.increase = 0.15, size = 5)) # Add pairwise comparisons p-value

stat.test <- subset(alpha, depth.watercolumn == 'bottom') %>%
  wilcox_test(faith_pd ~ group) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "group", dodge = 0.8)

(p2b=ggboxplot(subset(alpha, depth.watercolumn == 'bottom'), x = 'group', y = 'faith_pd',
               add = "jitter", fill = 'group', size = 1)+
    xlab("Phase") + ylab('Faith PD') + ggtitle('Bottom') +
    coord_cartesian(ylim = c(50,350))+
    theme(text = element_text(size = 22), legend.position = 'none')+
    scale_fill_jama(alpha=0.8) +
    stat_compare_means(method = "anova",label.x = 0.7,label.y = 50,size = 6)+
    stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0.05, step.increase = 0.1, size = 5)) # Add pairwise comparisons p-value

#evenness matrix
alpha = read_qza("R_analysis_input/evenness_vector.qza")
alpha <- alpha$data %>% tibble::rownames_to_column("SampleID") 
alpha = merge(x = alpha,
              y = metadata.o,
              by.x = 'SampleID',
              by.y = 0,
              all.Y = T)

stat.test <- subset(alpha, depth.watercolumn == 'top') %>%
  wilcox_test(pielou_e ~ group) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "group", dodge = 0.8)

(p3t=ggboxplot(subset(alpha, depth.watercolumn == 'top'), x = 'group', y = 'pielou_e',
               add = "jitter", fill = 'group', size = 1)+
    xlab("Phase") + ylab('Pielou evenness') + ggtitle('Top') +
    coord_cartesian(ylim = c(0.2,1))+
    theme(text = element_text(size = 22), legend.position = 'none')+
    scale_fill_lancet(alpha=0.8) +
    stat_compare_means(method = "anova",label.x = 0.7,label.y = 0.95,size = 6)+
    stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0.05, step.increase = 0.15, size = 5)) # Add the p-values for pairwise comparison

stat.test <- subset(alpha, depth.watercolumn == 'bottom') %>%
  wilcox_test(pielou_e ~ group) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "group", dodge = 0.8)

(p3b=ggboxplot(subset(alpha, depth.watercolumn == 'bottom'), x = 'group', y = 'pielou_e',
               add = "jitter", fill = 'group', size = 1)+
    xlab("Phase") + ylab('Pielou evenness') + ggtitle('Bottom') + 
    coord_cartesian(ylim = c(0.2,1))+
    theme(text = element_text(size = 22), legend.position = 'none')+
    scale_fill_jama(alpha=0.8) +
    stat_compare_means(method = "anova",label.x = 0.7,label.y = 0.2,size = 6)+
    stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0.05, step.increase = 0.2, size = 5)) # Add pairwise comparison p-values


#observed otus
alpha = read_qza("R_analysis_input/observed_otus_vector.qza")
alpha <- alpha$data %>% tibble::rownames_to_column("SampleID") 
alpha = merge(x = alpha,
              y = metadata.o,
              by.x = 'SampleID',
              by.y = 0,
              all.Y = T)

stat.test <- subset(alpha, depth.watercolumn == 'top') %>%
  wilcox_test(observed_otus ~ group) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "group", dodge = 0.8)

(p4t=ggboxplot(subset(alpha, depth.watercolumn == 'top'), x = 'group', y = 'observed_otus',
               add = "jitter", fill = 'group', size = 1)+
    xlab("Phase") + ylab('Observed ASVs') + ggtitle('Top') + 
    coord_cartesian(ylim = c(0,5000))+
    theme(text = element_text(size = 22), legend.position = 'none')+
    scale_fill_lancet(alpha=0.8) +
    stat_compare_means(method = "anova",label.x = 0.7,label.y = 0.2,size = 6)+
    stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0.05, step.increase = 0.1, size = 5))

stat.test <- subset(alpha, depth.watercolumn == 'bottom') %>%
  wilcox_test(observed_otus ~ group) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "group", dodge = 0.8)

(p4b=ggboxplot(subset(alpha, depth.watercolumn == 'bottom'), x = 'group', y = 'observed_otus',
               add = "jitter", fill = 'group', size = 1)+
    xlab("Phase") + ylab('Observed ASVs') + ggtitle('Bottom') + 
    coord_cartesian(ylim = c(0,5000))+
    theme(text = element_text(size = 22), legend.position = 'none')+
    scale_fill_jama(alpha=0.8) +
    stat_compare_means(method = "anova",label.x = 0.7,label.y = 0.2,size = 6)+
    stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0.05, step.increase = 0.2, size = 5, y.position = 3000)) 

plot_grid(p4t,p2t,p3t,p1t,
          p4b,p2b,p3b,p1b,
          nrow = 2)
ggsave(paste(PFIG,"/",'alpha_diversity.png',sep = ""),width = 22, height = 9, device = 'png')

##############PCoA plot using the Weighted UniFrac distance###############
#coerce depth.watercolumn and group in metadata.o to factors
metadata.o$depth.watercolumn = factor(metadata.o$depth.watercolumn)
#read the distance matrix output by QIIME2
dm = read_qza("R_analysis_input/weighted_unifrac_distance_matrix.qza")
#select the rows we need
dm = dist_subset(dm$data,rownames(metadata.o))
#some distance measures may result in negative eigenvalues. In that case, add a correction:
pcoa.wunifrac = pcoa(dm,correction = 'cailliez')
#transform the pcoa.wunifrac object into a df for making it plottable with ggplot2
coordinates.pcoa.wunifrac = data.frame(pcoa.wunifrac$vectors,check.names = FALSE)

###############plot for all samples#############
###############by depth#############
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
         subtitle = paste('PERMANOVA = ', padonis$Adonis.tab[6][[1]],sep = ''))+
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
ggsave(paste(PFIG,"/",'pcoa_all_depth.png',sep = ""),width = 6, height = 4, device = 'png')

###############by phases#############
#test between on and off phases
padonis = adonis_pairwise(metadata.o,dm,'group',permut = 999, adj.meth = "holm")
padonis$Adonis.tab
padonis$Betadisper.tab
#plot
color = metadata.o$group
(pcoa_all_phase=ggplot()+
    geom_point(data=subset(coordinates.pcoa.wunifrac[,c(1,2)]),
               mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = color),
               alpha = 0.8, size = 3) +
    labs(title = 'All samples',
         subtitle = paste('on VS off = ', padonis$Adonis.tab[6][[1]][2], ', ',
                          'on VS altered = ', padonis$Adonis.tab[6][[1]][1], ', ',
                          'altered VS off = ', padonis$Adonis.tab[6][[1]][3], sep = ''))+
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
    scale_color_manual(values = col21[2:4]))
ggsave(paste(PFIG,"/",'pcoa_all_phase.png',sep = ""),width = 6, height = 4, device = 'png')

###############plot for top samples by phases#############
#subset
metadata.o.s = subset(metadata.o, depth.watercolumn == 'top')
#repeat the pcoa step
dm = read_qza("R_analysis_input/weighted_unifrac_distance_matrix.qza")
dm = dist_subset(dm$data,rownames(metadata.o.s))
pcoa.wunifrac = pcoa(dm,correction = 'cailliez')
coordinates.pcoa.wunifrac = data.frame(pcoa.wunifrac$vectors,check.names = FALSE)

#test between on and off phases
padonis = adonis_pairwise(metadata.o.s,dm,'group',permut = 999, adj.meth = "holm")
padonis$Adonis.tab
padonis$Betadisper.tab
#plot
color = metadata.o.s$group
(pcoa_top=ggplot()+
    geom_point(data=subset(coordinates.pcoa.wunifrac[,c(1,2)]),
               mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = color),
               alpha = 0.8, size = 3) +
    labs(title = 'Top',
         subtitle = paste('on VS off = ', padonis$Adonis.tab[6][[1]][2], ', ',
                          'on VS altered = ', padonis$Adonis.tab[6][[1]][1], ', ',
                          'altered VS off = ', padonis$Adonis.tab[6][[1]][3], sep = ''))+
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
    scale_color_manual(values = col21[2:4]))
ggsave(paste(PFIG,"/",'pcoa_top_phase.png',sep = ""),width = 6, height = 4, device = 'png')

###############plot for bottom samples by phase#############
#subset
metadata.o.s = subset(metadata.o, depth.watercolumn == 'bottom')
#repeat the pcoa step
dm = read_qza("R_analysis_input/weighted_unifrac_distance_matrix.qza")
dm = dist_subset(dm$data,rownames(metadata.o.s))
pcoa.wunifrac = pcoa(dm,correction = 'cailliez')
coordinates.pcoa.wunifrac = data.frame(pcoa.wunifrac$vectors,check.names = FALSE)

#test between on and off phases
padonis = adonis_pairwise(metadata.o.s,dm,'group',permut = 999, adj.meth = "holm")
padonis$Adonis.tab
padonis$Betadisper.tab
#plot
color = metadata.o.s$group
(pcoa_bottom=ggplot()+
    geom_point(data=subset(coordinates.pcoa.wunifrac[,c(1,2)]),
               mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = color),
               alpha = 0.8, size = 3) +
    labs(title = 'Bottom',
         subtitle = paste('on VS off = ', padonis$Adonis.tab[6][[1]][2], ', ',
                          'on VS altered = ', padonis$Adonis.tab[6][[1]][1], ', ',
                          'altered VS off = ', padonis$Adonis.tab[6][[1]][3], sep = ''))+
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
    scale_color_manual(values = col21[2:4]))
ggsave(paste(PFIG,"/",'pcoa_bottom_phase.png',sep = ""),width = 6, height = 4, device = 'png')

#################compare phylogenetic distance between sequencing run replicates and samples at different depths and stations###############
###########weighted unifrac distance analysis###################
#read the distance matrix output by QIIME2
dm_wu = read_qza("R_analysis_input/weighted_unifrac_distance_matrix.qza")

#select the rows we need
dm_wu = dist_subset(dm_wu$data,rownames(metadata.o))
#convert the distance object into a dataframe showing pairwise distance
dm_wu = dist2list(dm_wu)

#replace the indices to meaningful sample names
dm_wu$col = SIDconversion$SD[match(dm_wu$col,SIDconversion$SIDmetadata)]
dm_wu$row = SIDconversion$SD[match(dm_wu$row,SIDconversion$SIDmetadata)]
dm_wu$col = as.character(dm_wu$col)
dm_wu$row = as.character(dm_wu$row)

##calculate pairwise distances between depths
dm <- dm_wu
dm$col.suffix = substr(dm$col,3,nchar(dm$col))
dm$row.suffix = substr(dm$row,3,nchar(dm$row))
dm = subset(dm, dm$col.suffix == dm$row.suffix & dm$col != dm$row)
dm %>% distinct(`col.suffix`,.keep_all = T) -> dm
dm_depth = c('depth', mean(dm$value), sd(dm$value))

##calculate distances between sampling points at the same station and depth
dm <- dm_wu
dm$col.ds = substr(dm$col,1,5)
dm$row.ds = substr(dm$row,1,5)
dm = subset(dm, dm$col.ds == dm$row.ds & dm$col != dm$row)
dm_time = c('sampling time', mean(dm$value), sd(dm$value))

#calculate pairwise distances between stations at the same sampling time point and depth
dm <- dm_wu
dm$col.ddt = paste(substr(dm$col,1,1),substr(dm$col,7,nchar(dm$col)))
dm$row.ddt = paste(substr(dm$row,1,1),substr(dm$row,7,nchar(dm$row)))
dm = subset(dm, dm$col.ddt == dm$row.ddt & dm$col != dm$row)
dm_station = c('station', mean(dm$value), sd(dm$value))

##calculate distance between replicate samples of the two sequencing runs
#read the distance matrix output by QIIME2
dm_wu = read_qza("R_analysis_input/weighted_unifrac_distance_matrix.qza")
#read the metadata of the replicate samples
rplcs = read.csv('R_analysis_input/replicates_between_sequencing_runs.csv')
#select the distance we need
dm_wu = dist_subset(dm_wu$data,rplcs$SID)
#convert the distance object into a dataframe showing pairwise distance
dm_wu = dist2list(dm_wu)
#check if there is self-comparison (there shouldn't be)
subset(dm_wu, dm_wu$col == dm_wu$row)

#replace the indices to meaningful sample names
dm_wu$col = rplcs$Name[match(dm_wu$col,rplcs$SID)]
dm_wu$row = rplcs$Name[match(dm_wu$row,rplcs$SID)]

#retain the rows of the distance between sequencing replicates
dm <- dm_wu
dm = subset(dm,dm$col == dm$row)

#calculate mean and std
dm_replicate = c('replicates', mean(dm$value), sd(dm$value))

#summarize all numbers into one df
rplc_cp = data.frame(name = character(),
                     mean = numeric(),
                     sd = numeric())
rplc_cp[1,] = dm_depth
rplc_cp[2,] = dm_station
rplc_cp[3,] = dm_time
rplc_cp[4,] = dm_replicate

summary(rplc_cp)
rplc_cp$mean = as.numeric(rplc_cp$mean)
rplc_cp$sd = as.numeric(rplc_cp$sd)

#plot it!
rplc_cp$name = factor(rplc_cp$name, level = c('station','depth','sampling time','replicates'))
ggplot(rplc_cp) +
  geom_bar(aes(x = name, y = mean, fill = name), stat="identity", alpha = 0.6) +
  geom_errorbar(aes(x = name, ymin = mean - sd, ymax = mean + sd), width = 0.2)+
  scale_fill_npg()+
  xlab('')+
  ylab('Weighted UniFrac distance')+
  theme(axis.text.x = element_text(size=12, colour = "black"), 
        axis.text.y = element_text(size=15, colour = "black"), 
        axis.title=element_text(size=15),
        axis.line = element_line(size=0.5),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length.y = unit(.25, "cm"),
        axis.ticks.length.x = unit(.25, "cm"),
        plot.margin=unit(c(0.5,1,0.1,0.2),"cm"),
        legend.position = 'none')
ggsave(paste(PFIG,'/variability_between_replicates.png', sep = ''),width = 7,height = 4,device = 'png')

################descriptive statistics/rarefaction curves##################
#read the ASV/reads counts output by qimme2
freq = read.csv('R_analysis_input/combined_dada2_summarize.qzv.frequency_per_sample.csv')
#read the metadata file of all samples (160 samples sequenced for the project)
metadata160 = read.csv('R_analysis_input/metadata_160.csv')

metadata160 %>%
  subset(.,select = c(sample.id,sample.name, type), type == 'startblank' | type == 'negative') %>%
  merge(x = .,y = freq, by = 'sample.id', all.x = T) ->freq.neg

metadata160 %>%
  subset(.,select = c(sample.id,sample.name, type), sample.id %in% SIDconversion$SIDmetadata) %>%
  merge(x = .,y = freq, by = 'sample.id', all.x = T) ->freq.115

mean(freq.neg$frequency[freq.neg$type == 'negative'])
sd(freq.neg$frequency[freq.neg$type == 'negative'])

mean(freq.115$frequency)
sd(freq.115$frequency)
sum(freq.115$frequency)

#check the total number of reads in the asv table
pu_115<-qza_to_phyloseq(
  features="R_analysis_input/unrarefied_table.qza",
  tree="R_analysis_input/rooted-tree.qza",
  taxonomy="R_analysis_input/taxonomy.qza",
  paste(FAS,'metadata_115_preprocessed.csv',sep = '/'))
asv_115_unrarefied = data.frame(otu_table(pu_115))

#calculate the number of ASVs generated in the unrarefied table = 27303 (after removing the ghost ASVs)
asv_115_unrarefied.rmghost = asv_115_unrarefied[rowSums(asv_115_unrarefied)>0,]
nrow(asv_115_unrarefied.rmghost)

#####plot the rarefaction curve##########
#rarecurve(t(asv_115_unrarefied),step = 1000, label = F)
#ggsave(paste(PFIG,"/",'rarefaction.png',sep = ""),width = 4, height = 4, device = 'png')

#check the number of reads generated/passed filter
dada2.seq2 = read.csv('R_analysis_input/yue_rockcreek_16s_seq2_analysis-stats-dada2.tsv', sep = '\t')
dada2.seq1 = read.csv('R_analysis_input/yue_rockcreek_16s_seq1_177490_3_analysis-stats-dada2.tsv', sep = '\t')
dada2 = bind_rows(dada2.seq1,dada2.seq2)
dada2.115 = subset(dada2, sample.id %in% SIDconversion$SIDmetadata)

#total num. of reads generated = 35527527
sum(as.numeric(dada2.115$input))
#total num. of reads passed the filter = 31960245, matched the number in the asv table
sum(as.numeric(dada2.115$non.chimeric))

#########weighted unifrac distance pcoa plot to show differences between samples and controls############
metadata.115_n_controls = subset(metadata160, type == 'startblank' | type == 'zymo' | sample.id %in% rownames(metadata.o), select = c(sample.id,type))
metadata.115_n_controls$type[metadata.115_n_controls$type == 'startblank'] = 'negative'
metadata.115_n_controls$type[metadata.115_n_controls$type == 'zymo'] = 'positive'
#read the distance matrix output by QIIME2
dm = read_qza("R_analysis_input/weighted_unifrac_distance_matrix.qza")
#select the entries we need
dm = dist_subset(dm$data,metadata.115_n_controls$sample.id)
# Some distance measures may result in negative eigenvalues. In that case, add a correction:
pcoa.wunifrac = pcoa(dm,correction = 'cailliez')
#transform the pcoa.wunifrac object to a df for making it able to be plotted using ggplot2
coordinates.pcoa.wunifrac = data.frame(pcoa.wunifrac$vectors,check.names = FALSE)
#plot
color = metadata.115_n_controls$type
(pcoa_all_depth=ggplot()+
    geom_point(data=subset(coordinates.pcoa.wunifrac[,c(1,2)]),
               mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = color),
               alpha = 0.8, size = 3) + 
    xlab(paste('PC1 (',label_percent()(pcoa.wunifrac$values$Rel_corr_eig[1]),' variability explained)',sep = '')) +
    ylab(paste('PC2 (',label_percent()(pcoa.wunifrac$values$Rel_corr_eig[2]),' variability explained)',sep = '')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = 'black'),
          legend.title=element_blank(),
          legend.position = 'bottom',
          legend.text = element_text(size=15))+
    scale_color_jama())
ggsave(paste(PFIG,"/",'pcoa_controls.png',sep = ""),width = 4, height = 4, device = 'png')

###############Session information############
sessionInfo()
devtools::session_info()
