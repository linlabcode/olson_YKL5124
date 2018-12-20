### YKL UPREGULATED RADAR PLOT ###

### INSTALLING PACKAGES FOR RADAR PLOTS ###
install.packages('fmsb')
library(fmsb)

### IMPORTING CUFFNORM DATA ###
setwd('/Users/DanPark1116/Desktop/BCM_Materials/Year_Two/CDK7_Analysis/CDK7_GRAY_cuffnorm/log2fc_pvals/')
ykl_dmso_wt_low = read.delim('CDK7_GRAY_HAP1_WT_DMSO_R_vs_HAP1_WT_YKL_LO_R_exprs_matrix.txt', stringsAsFactors = F)
thz531_dmso_wt = read.delim('CDK7_GRAY_HAP1_WT_531_R_vs_HAP1_WT_DMSO_R_exprs_matrix.txt', stringsAsFactors = F)
thz531_dmso_wt$LOG2_FOLD_CHANGE = thz531_dmso_wt$LOG2_FOLD_CHANGE * -1
thz1_dmso_wt = read.delim('CDK7_GRAY_HAP1_WT_DMSO_R_vs_HAP1_WT_THZ1_R_exprs_matrix.txt', stringsAsFactors = F)
combo_dmso_wt_low = read.delim('CDK7_GRAY_HAP1_WT_COMBO_LO_R_vs_HAP1_WT_DMSO_R_exprs_matrix.txt', stringsAsFactors = F)
combo_dmso_wt_low$LOG2_FOLD_CHANGE = combo_dmso_wt_low$LOG2_FOLD_CHANGE * -1

ykl_dmso_mut_low = read.delim('CDK7_GRAY_HAP1_MUT_DMSO_R_vs_HAP1_MUT_YKL_LO_R_exprs_matrix.txt', stringsAsFactors = F)
thz531_dmso_mut = read.delim('CDK7_GRAY_HAP1_MUT_531_R_vs_HAP1_MUT_DMSO_R_exprs_matrix.txt', stringsAsFactors = F)
thz531_dmso_mut$LOG2_FOLD_CHANGE = thz531_dmso_mut$LOG2_FOLD_CHANGE * -1
thz1_dmso_mut = read.delim('CDK7_GRAY_HAP1_MUT_DMSO_R_vs_HAP1_MUT_THZ1_R_exprs_matrix.txt', stringsAsFactors = F)
combo_dmso_mut_low = read.delim('CDK7_GRAY_HAP1_MUT_COMBO_LO_R_vs_HAP1_MUT_DMSO_R_exprs_matrix.txt', stringsAsFactors = F)
combo_dmso_mut_low$LOG2_FOLD_CHANGE = combo_dmso_mut_low$LOG2_FOLD_CHANGE * -1

### ISOLATING GENES THAT ARE UPREGULATED IN YKL VS DMSO ###
ykl_dmso_wt_upreg = ykl_dmso_wt_low[ykl_dmso_wt_low$LOG2_FOLD_CHANGE > 1,]
ykl_dmso_wt_upreg = ykl_dmso_wt_upreg[ykl_dmso_wt_upreg$P_VALUE < 0.05,]

ykl_dmso_wt_upreg_background = ykl_dmso_wt_low[!rownames(ykl_dmso_wt_low) %in% rownames(ykl_dmso_wt_upreg), ]

### GENERATING RADAR PLOTS FOR UPREGULATED GENES
ykl_dmso_wt_upreg_genes = rownames(ykl_dmso_wt_upreg)

thz531_dmso_wt_upreg = thz531_dmso_wt[rownames(thz531_dmso_wt) %in% ykl_dmso_wt_upreg_genes,]
thz1_dmso_wt_upreg = thz1_dmso_wt[rownames(thz1_dmso_wt) %in% ykl_dmso_wt_upreg_genes,]
combo_dmso_wt_upreg = combo_dmso_wt_low[rownames(combo_dmso_wt_low) %in% ykl_dmso_wt_upreg_genes,]

mean(ykl_dmso_wt_upreg$LOG2_FOLD_CHANGE) #1.54
mean(thz531_dmso_wt_upreg$LOG2_FOLD_CHANGE) #-0.02
mean(thz1_dmso_wt_upreg$LOG2_FOLD_CHANGE) #0.74
mean(combo_dmso_wt_upreg$LOG2_FOLD_CHANGE) #1.27

ykl_upreg_log2fc_values = data.frame(c(1.5, -0.02, 0.74, 1.3))
ykl_upreg_log2fc_values = t(ykl_upreg_log2fc_values)
colnames(ykl_upreg_log2fc_values) = c('YKL', 'THZ531', 'THZ1', 'COMBO')
ykl_upreg_log2fc_values = as.data.frame(ykl_upreg_log2fc_values)
rownames(ykl_upreg_log2fc_values) = 'Log2FC'

ykl_upreg_log2fc_values[2,] = -0.2
ykl_upreg_log2fc_values[3,] = 1.6

ykl_upreg_log2fc_values = ykl_upreg_log2fc_values[c(3,2,1),]

### CREATING RADAR PLOT FOR GENES UPREGULATED WITH YKL ###
par(mar=c(0.5, 0.5, 0.5, 0.5))
radarchart( ykl_upreg_log2fc_values  , axistype=1 ,
            
            #custom polygon
            pcol=rgb(1, 0.2, 0.2, 0.8) , pfcol=rgb(1, 0.2, 0.2, 0.5) , plwd=4 ,
            
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1.6,0.4), cglwd=0.8
            
            #custom labels
            #vlcex=0.8
)

par(new = T)

ykl_dmso_mut_upreg = ykl_dmso_mut[rownames(ykl_dmso_mut) %in% ykl_dmso_wt_upreg_genes,]
thz531_dmso_mut_upreg = thz531_dmso_mut[rownames(thz531_dmso_mut) %in% ykl_dmso_wt_upreg_genes,]
thz1_dmso_mut_upreg = thz1_dmso_mut[rownames(thz1_dmso_mut) %in% ykl_dmso_wt_upreg_genes,]
combo_dmso_mut_upreg = combo_dmso_mut[rownames(combo_dmso_mut) %in% ykl_dmso_wt_upreg_genes,]

mean(ykl_dmso_mut_upreg$LOG2_FOLD_CHANGE) #0.12
mean(thz531_dmso_mut_upreg$LOG2_FOLD_CHANGE) #-0.14
mean(thz1_dmso_mut_upreg$LOG2_FOLD_CHANGE) #-0.10
mean(combo_dmso_mut_upreg$LOG2_FOLD_CHANGE) #-0.04

ykl_upreg_mut_log2fc_values = data.frame(c(0.12, -0.14, -0.10, -0.04))
ykl_upreg_mut_log2fc_values = t(ykl_upreg_mut_log2fc_values)
colnames(ykl_upreg_mut_log2fc_values) = c('YKL', 'THZ531', 'THZ1', 'COMBO')
ykl_upreg_mut_log2fc_values = as.data.frame(ykl_upreg_mut_log2fc_values)
rownames(ykl_upreg_mut_log2fc_values) = 'Log2FC'

ykl_upreg_mut_log2fc_values[2,] = -0.2
ykl_upreg_mut_log2fc_values[3,] = 1.6

ykl_upreg_mut_log2fc_values = ykl_upreg_mut_log2fc_values[c(3,2,1),]

radarchart( ykl_upreg_mut_log2fc_values  , axistype=1 ,
            
            #custom polygon
            pcol=rgb(0, 0, 0, 0.3) , pfcol=rgb(0,0,0,0.5) , plwd=4 ,
            
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1.6,0.4), cglwd=0.8
            
            #custom labels
            #vlcex=0.8
)
legend(1.0,-0.5,
       legend=c("WT","Mut"),
       pch=c(15,15),
       col=c("red","grey"))







