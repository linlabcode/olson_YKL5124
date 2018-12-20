### CELL CYCLE GENES BOXPLOT ###

### IMPORTING DATA ###
cycle_genes = read.csv('/Users/DanPark1116/Desktop/BCM_Materials/Year_Two/CDK7_Analysis/CDK7_GRAY_cuffnorm/cycle_genes_e2f3.csv', stringsAsFactors = F)

### ISOLATE YKL DOSAGE DATA ###
cycle_genes_ykl_dosages_wt = cycle_genes[, c(1:6, 16:18)]
cycle_genes_ykl_dosages_mut = cycle_genes[, c(22:27, 37:39)]

### GENERATE AVERAGES AND BARPLOT ###
cycle_genes_ykl_dosages_wt$dmso_avg = rowMeans(cycle_genes_ykl_dosages_wt[, c(7:9)])
cycle_genes_ykl_dosages_wt$low_ykl_avg = rowMeans(cycle_genes_ykl_dosages_wt[, c(4:6)])
cycle_genes_ykl_dosages_wt$high_ykl_avg = rowMeans(cycle_genes_ykl_dosages_wt[, c(1:3)])
cycle_genes_ykl_dosages_wt_avgs = cycle_genes_ykl_dosages_wt[, c(10:12)]
cycle_genes_ykl_dosages_wt_avgs_norpl = cycle_genes_ykl_dosages_wt_avgs[!grepl('RPL', rownames(cycle_genes_ykl_dosages_wt_avgs)),]

cycle_genes_ykl_dosages_mut$dmso_avg = rowMeans(cycle_genes_ykl_dosages_mut[, c(7:9)])
cycle_genes_ykl_dosages_mut$low_ykl_avg = rowMeans(cycle_genes_ykl_dosages_mut[, c(4:6)])
cycle_genes_ykl_dosages_mut$high_ykl_avg = rowMeans(cycle_genes_ykl_dosages_mut[, c(1:3)])
cycle_genes_ykl_dosages_mut_avgs = cycle_genes_ykl_dosages_mut[, c(10:12)]
cycle_genes_ykl_dosages_mut_avgs_norpl = cycle_genes_ykl_dosages_mut_avgs[!grepl('RPL', rownames(cycle_genes_ykl_dosages_mut_avgs)),]

mean(cycle_genes_ykl_dosages_wt_avgs$dmso_avg) #149.53
mean(cycle_genes_ykl_dosages_wt_avgs$low_ykl_avg) #92.59
mean(cycle_genes_ykl_dosages_wt_avgs$high_ykl_avg) #82.03

mean(cycle_genes_ykl_dosages_mut_avgs$dmso_avg) #206.45
mean(cycle_genes_ykl_dosages_mut_avgs$low_ykl_avg) #145.20
mean(cycle_genes_ykl_dosages_mut_avgs$high_ykl_avg) #142.00

mean(cycle_genes_ykl_dosages_wt_avgs_norpl$dmso_avg) #90.5
sd(cycle_genes_ykl_dosages_wt_avgs_norpl$dmso_avg) #80.9
mean(cycle_genes_ykl_dosages_wt_avgs_norpl$low_ykl_avg) #55.9
sd(cycle_genes_ykl_dosages_wt_avgs_norpl$low_ykl_avg) #49.04
mean(cycle_genes_ykl_dosages_wt_avgs_norpl$high_ykl_avg) #45.8
sd(cycle_genes_ykl_dosages_wt_avgs_norpl$high_ykl_avg) #37.9

### BOXPLOT COMPARISON OF CYCLE GENES IN WT/MUT ###
cycle_dmso_wt = cycle_genes_ykl_dosages_wt_avgs_norpl[, c(1,2)]
colnames(cycle_dmso_wt) = c('Expression', 'Treatment')
cycle_dmso_wt$Treatment = 'WT_DMSO'

cycle_ykl_low_wt = cycle_genes_ykl_dosages_wt_avgs_norpl[, c(2,3)]
colnames(cycle_ykl_low_wt) = c('Expression', 'Treatment')
cycle_ykl_low_wt$Treatment = 'WT_YKL_Low'

cycle_ykl_high_wt = cycle_genes_ykl_dosages_wt_avgs_norpl[, c(3,2)]
colnames(cycle_ykl_high_wt) = c('Expression', 'Treatment')
cycle_ykl_high_wt$Treatment = 'WT_YKL_High'

cycle_dmso_mut = cycle_genes_ykl_dosages_mut_avgs_norpl[, c(1,2)]
colnames(cycle_dmso_mut) = c('Expression', 'Treatment')
cycle_dmso_mut$Treatment = 'Mut_DMSO'

cycle_ykl_low_mut = cycle_genes_ykl_dosages_mut_avgs_norpl[, c(2,3)]
colnames(cycle_ykl_low_mut) = c('Expression', 'Treatment')
cycle_ykl_low_mut$Treatment = 'Mut_YKL_Low'

cycle_ykl_high_mut = cycle_genes_ykl_dosages_mut_avgs_norpl[, c(3,2)]
colnames(cycle_ykl_high_mut) = c('Expression', 'Treatment')
cycle_ykl_high_mut$Treatment = 'Mut_YKL_High'

cycle_wt_mut_dosage_rbind = rbind(cycle_dmso_wt, cycle_ykl_low_wt, cycle_ykl_high_wt, cycle_dmso_mut, cycle_ykl_low_mut, cycle_ykl_high_mut)
cycle_df1 = data.frame(a = c(1,1,2,2), b = c(275, 300, 300, 275))
cycle_df2 = data.frame(a = c(1,1,3,3), b = c(325, 350, 350, 325))
cycle_df3 = data.frame(a = c(4,4,5,5), b = c(300, 325, 325, 300))
cycle_df4 = data.frame(a = c(4,4,6,6), b = c(350, 375, 375, 350))

cycle_wt_mut_dosage_rbind$Treatment = factor(cycle_wt_mut_dosage_rbind$Treatment, levels = c('WT_DMSO', 'WT_YKL_Low', 'WT_YKL_High', 'Mut_DMSO', 'Mut_YKL_Low', 'Mut_YKL_High'), ordered = T)

library(ggplot2)
ggplot(cycle_wt_mut_dosage_rbind, aes(x=factor(Treatment), y=Expression)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(data = cycle_df1, aes(x = a, y = b)) + annotate('text', x = 1.5, y = 305, label = '**', size = 8)+
  geom_line(data = cycle_df2, aes(x = a, y = b)) + annotate('text', x = 2.0, y = 355, label = '***', size = 8)+
  geom_line(data = cycle_df3, aes(x = a, y = b)) + annotate('text', x = 4.5, y = 335, label = 'ns', size = 6)+
  geom_line(data = cycle_df4, aes(x = a, y = b)) + annotate('text', x = 5.0, y = 385, label = 'ns', size = 6)+
  labs(x = 'Cell Line and Treatment') +
  ggtitle('Expression Levels for Cell Cycle Genes in WT/Mutant Cell Lines (n = 53)') +
  #geom_hline(yintercept = 0, linetype = 'dotted')+
  theme(text = element_text(size = 15)) +
  ylim(0, 400)

wilcox.test(cycle_genes_ykl_dosages_wt_avgs_norpl$dmso_avg, cycle_genes_ykl_dosages_wt_avgs_norpl$low_ykl_avg) #0.0105
wilcox.test(cycle_genes_ykl_dosages_wt_avgs_norpl$dmso_avg, cycle_genes_ykl_dosages_wt_avgs_norpl$high_ykl_avg) #0.00059
wilcox.test(cycle_genes_ykl_dosages_wt_avgs_norpl$low_ykl_avg, cycle_genes_ykl_dosages_wt_avgs_norpl$high_ykl_avg) #0.309

wilcox.test(cycle_genes_ykl_dosages_mut_avgs_norpl$dmso_avg, cycle_genes_ykl_dosages_mut_avgs_norpl$low_ykl_avg) #0.32
wilcox.test(cycle_genes_ykl_dosages_mut_avgs_norpl$dmso_avg, cycle_genes_ykl_dosages_mut_avgs_norpl$high_ykl_avg) #0.29
wilcox.test(cycle_genes_ykl_dosages_mut_avgs_norpl$low_ykl_avg, cycle_genes_ykl_dosages_mut_avgs_norpl$high_ykl_avg) #0.98



