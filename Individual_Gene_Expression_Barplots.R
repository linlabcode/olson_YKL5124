### IMPORTING DATA ###
gene_data = read.csv('/Users/DanPark1116/Desktop/BCM_Materials/Year_Two/CDK7_Analysis/Initial_Output/cdk7_genes_fpkm.csv', stringsAsFactors = F)
gene_data = gene_data[!duplicated(gene_data$tracking_id),] #removes 2
rownames(gene_data) = gene_data[,1]
gene_data = gene_data[, -1]

### FILTERING FOR EXPRESSED GENES ###
gene_data$max = apply(gene_data, 1, max)
gene_data = gene_data[gene_data$max > 10, -43]

### ISOLATING DATA FOR DIFFERENT DOSES OF YKL ###
dose_data = gene_data[, c(1:6, 16:18, 22:27, 37:39)]
dose_data$gene_id = rownames(dose_data)

### CALCULATING AVERAGE EXPRESSION AT EACH DOSAGE ###
dose_data$WT_YKL_2_Avg = rowMeans(dose_data[, c(1:3)])
dose_data$WT_YKL_500_Avg = rowMeans(dose_data[, c(4:6)])
dose_data$WT_DMSO_Avg = rowMeans(dose_data[, c(7:9)])
dose_data$Mut_YKL_2_Avg = rowMeans(dose_data[, c(10:12)])
dose_data$Mut_YKL_500_Avg = rowMeans(dose_data[, c(13:15)])
dose_data$Mut_DMSO_Avg = rowMeans(dose_data[, c(16:18)])

dose_data = dose_data[, c(19:25)]

### CDCA3 INDIVIDUAL GENE EXPRESSION BARPLOT ###
CDCA3_data = dose_data[dose_data$gene_id == 'CDCA3', -1]
CDCA3_data = CDCA3_data[, c(3,2,1,6,5,4)]

dose_sd = gene_data[, c(1:6, 16:18, 22:27, 37:39)]
CDCA3_ykl_wt_high_sd = sd(dose_sd[rownames(dose_sd) == 'CDCA3', c(1:3)]) 
CDCA3_ykl_wt_high_sd #1.03
CDCA3_ykl_wt_low_sd = sd(dose_sd[rownames(dose_sd) == 'CDCA3', c(4:6)])
CDCA3_ykl_wt_low_sd #1.01
CDCA3_ykl_wt_dmso_sd = sd(dose_sd[rownames(dose_sd) == 'CDCA3', c(7:9)])
CDCA3_ykl_wt_dmso_sd #1.32
CDCA3_ykl_mut_high_sd = sd(dose_sd[rownames(dose_sd) == 'CDCA3', c(10:12)])
CDCA3_ykl_mut_high_sd #1.44
CDCA3_ykl_mut_low_sd = sd(dose_sd[rownames(dose_sd) == 'CDCA3', c(13:15)])
CDCA3_ykl_mut_low_sd #0.415
CDCA3_ykl_mut_dmso_sd = sd(dose_sd[rownames(dose_sd) == 'CDCA3', c(16:18)])
CDCA3_ykl_mut_dmso_sd #1.42

CDCA3_data = t(CDCA3_data)
CDCA3_SD = data.frame(1.32, 1.01, 1.03, 1.42, 0.415, 1.44)
colnames(CDCA3_SD) = c('WT_DMSO_Avg', 'WT_YKL_500_Avg', 'WT_YKL_2_Avg', 'Mut_DMSO_Avg', 'Mut_YKL_500_Avg', 'Mut_YKL_2_Avg')
CDCA3_SD = t(CDCA3_SD)

CDCA3_data_final = cbind(CDCA3_data, CDCA3_SD)
colnames(CDCA3_data_final) = c('Averages', 'SD')
CDCA3_data_final = as.data.frame(CDCA3_data_final)
CDCA3_data_final$Cell_Line = 'WT'
CDCA3_data_final[c(4:6), 3] = 'Mutant'
colnames(CDCA3_data_final)[3] = 'Cell Line'
CDCA3_data_final$Treatment = rownames(CDCA3_data_final)
CDCA3_data_final[1,4] = 'WT DMSO'
CDCA3_data_final[2,4] = 'WT YKL 1'
CDCA3_data_final[3,4] = 'WT YKL 2'
CDCA3_data_final[4,4] = 'Mut DMSO'
CDCA3_data_final[5,4] = 'Mut YKL 1'
CDCA3_data_final[6,4] = 'Mut YKL 2'

library(ggplot2)

CDCA3_data_final$Treatment = factor(CDCA3_data_final$Treatment,
                                    levels = c('WT DMSO', 'WT YKL 1', 'WT YKL 2', 'Mut DMSO', 'Mut YKL 1', 'Mut YKL 2'), ordered = T)

cdca3_plot = ggplot(CDCA3_data_final, aes(fill=`Cell Line`,x=Treatment, y=Averages)) +
  geom_bar(stat="identity", position = 'dodge') +
  geom_errorbar(aes(ymin = Averages - SD, ymax = Averages + SD), width=0.6) +
  theme_classic() +
  labs(x = 'Treatments and Dose (1 = Lower; 2 = Higher)', title = 'CDCA3 FPKM For Different YKL Dosages in WT/Mutant Cell Lines', y = 'FPKM Averages Across Triplicate') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c('grey', 'red')) +
  theme(text = element_text(size = 16))

ggsave('/Users/DanPark1116/Desktop/BCM_Materials/Year_Two/CDK7_Analysis/Figures/cdca3_fpkm_barplot_formatted.pdf', cdca3_plot, width = 25, height = 15, device = 'pdf', units = 'cm')

### PLK1 INDIVIDUAL GENE EXPRESSION BARPLOT ###
PLK1_data = dose_data[dose_data$gene_id == 'PLK1', -1]
PLK1_data = PLK1_data[, c(3,2,1,6,5,4)]

dose_sd = gene_data[, c(1:6, 16:18, 22:27, 37:39)]
PLK1_ykl_wt_high_sd = sd(dose_sd[rownames(dose_sd) == 'PLK1', c(1:3)]) 
PLK1_ykl_wt_high_sd #4.39
PLK1_ykl_wt_low_sd = sd(dose_sd[rownames(dose_sd) == 'PLK1', c(4:6)])
PLK1_ykl_wt_low_sd #1.39
PLK1_ykl_wt_dmso_sd = sd(dose_sd[rownames(dose_sd) == 'PLK1', c(7:9)])
PLK1_ykl_wt_dmso_sd #7.04
PLK1_ykl_mut_high_sd = sd(dose_sd[rownames(dose_sd) == 'PLK1', c(10:12)])
PLK1_ykl_mut_high_sd #10.12
PLK1_ykl_mut_low_sd = sd(dose_sd[rownames(dose_sd) == 'PLK1', c(13:15)])
PLK1_ykl_mut_low_sd #3.96
PLK1_ykl_mut_dmso_sd = sd(dose_sd[rownames(dose_sd) == 'PLK1', c(16:18)])
PLK1_ykl_mut_dmso_sd #10.86

PLK1_data = t(PLK1_data)
PLK1_SD = data.frame(7.04, 1.39, 4.39, 10.12, 3.96, 10.86)
colnames(PLK1_SD) = c('WT_DMSO_Avg', 'WT_YKL_500_Avg', 'WT_YKL_2_Avg', 'Mut_DMSO_Avg', 'Mut_YKL_500_Avg', 'Mut_YKL_2_Avg')
PLK1_SD = t(PLK1_SD)

PLK1_data_final = cbind(PLK1_data, PLK1_SD)
colnames(PLK1_data_final) = c('Averages', 'SD')
PLK1_data_final = as.data.frame(PLK1_data_final)
PLK1_data_final$Cell_Line = 'WT'
PLK1_data_final[c(4:6), 3] = 'Mutant'
colnames(PLK1_data_final)[3] = 'Cell Line'
PLK1_data_final$Treatment = rownames(PLK1_data_final)
PLK1_data_final[1,4] = 'WT DMSO'
PLK1_data_final[2,4] = 'WT YKL 1'
PLK1_data_final[3,4] = 'WT YKL 2'
PLK1_data_final[4,4] = 'Mut DMSO'
PLK1_data_final[5,4] = 'Mut YKL 1'
PLK1_data_final[6,4] = 'Mut YKL 2'

PLK1_data_final$Treatment = factor(PLK1_data_final$Treatment,
                                   levels = c('WT DMSO', 'WT YKL 1', 'WT YKL 2', 'Mut DMSO', 'Mut YKL 1', 'Mut YKL 2'), ordered = T)

plk1_plot = ggplot(PLK1_data_final, aes(fill=`Cell Line`,x=Treatment, y=Averages)) +
  geom_bar(stat="identity", position = 'dodge') +
  geom_errorbar(aes(ymin = Averages - SD, ymax = Averages + SD), width=0.6) +
  theme_classic() +
  labs(x = 'Treatments and Dose (1 = Lower; 2 = Higher)', title = 'PLK1 FPKM For Different YKL Dosages in WT/Mutant Cell Lines', y = 'FPKM Averages Across Triplicate') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c('grey', 'red')) +
  theme(text = element_text(size = 16))

ggsave('/Users/DanPark1116/Desktop/BCM_Materials/Year_Two/CDK7_Analysis/Figures/plk1_fpkm_barplot_formatted.pdf', plk1_plot, width = 25, height = 15, device = 'pdf', units = 'cm')





