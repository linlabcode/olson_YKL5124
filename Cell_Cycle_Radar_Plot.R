### RADAR PLOT FOR CELL CYCLE GENES ###

### INSTALLING PACKAGES FOR RADAR PLOTS ###
install.packages('fmsb')
library(fmsb)

### IMPORTING EXPRESSION DATA ###
gene_data = read.csv('/Users/DanPark1116/Desktop/BCM_Materials/Year_Two/CDK7_Analysis/Initial_Output/cdk7_genes_fpkm.csv', stringsAsFactors = F)
gene_data = gene_data[!duplicated(gene_data$tracking_id),] #removes 2
rownames(gene_data) = gene_data[,1]
gene_data = gene_data[, -1]

### FILTERING FOR EXPRESSED GENES ###
gene_data$max = apply(gene_data, 1, max)
gene_data = gene_data[gene_data$max > 10, -43]

### ISOLATING ENRICHED CELL CYCLE GENES ###
cell_cycle_genes = read.delim('/Users/DanPark1116/Desktop/KONG_E2F3_TARGETS.txt', stringsAsFactors = F)
cell_cycle_genes = cell_cycle_genes[cell_cycle_genes$CORE.ENRICHMENT == 'Yes',]
cell_cycle_genes = cell_cycle_genes[,2]

cycle_genes = gene_data[rownames(gene_data) %in% cell_cycle_genes,]

cycle_genes_ykl = cycle_genes[, c(4:6, 16:18)]
cycle_genes_thz531 = cycle_genes[, c(7:9, 16:18)]
cycle_genes_thz1 = cycle_genes[, c(19:21, 16:18)]
cycle_genes_combo = cycle_genes[, c(13:15, 16:18)]

### LOOKING AT LOG2FC FOR CYCLE GENES ACROSS ALL TREATMENTS - RADAR PLOT ###
cycle_genes_ykl$ykl_avg = rowMeans(cycle_genes_ykl[, c(1:3)])
cycle_genes_ykl$dmso_avg = rowMeans(cycle_genes_ykl[, c(4:6)])
cycle_genes_ykl$log2fc = log2(cycle_genes_ykl$ykl_avg / cycle_genes_ykl$dmso_avg)
mean(cycle_genes_ykl$log2fc) #-0.69

cycle_genes_thz531$thz531_avg = rowMeans(cycle_genes_thz531[, c(1:3)])
cycle_genes_thz531$dmso_avg = rowMeans(cycle_genes_thz531[, c(4:6)])
cycle_genes_thz531$log2fc = log2(cycle_genes_thz531$thz531_avg / cycle_genes_thz531$dmso_avg)
mean(cycle_genes_thz531$log2fc) #-0.017

cycle_genes_thz1$thz1_avg = rowMeans(cycle_genes_thz1[, c(1:3)])
cycle_genes_thz1$dmso_avg = rowMeans(cycle_genes_thz1[, c(4:6)])
cycle_genes_thz1$log2fc = log2(cycle_genes_thz1$thz1_avg / cycle_genes_thz1$dmso_avg)
mean(cycle_genes_thz1$log2fc) #-0.36

cycle_genes_combo$combo_avg = rowMeans(cycle_genes_combo[, c(1:3)])
cycle_genes_combo$dmso_avg = rowMeans(cycle_genes_combo[, c(4:6)])
cycle_genes_combo$log2fc = log2(cycle_genes_combo$combo_avg / cycle_genes_combo$dmso_avg)
mean(cycle_genes_combo$log2fc) #-0.587

cycle_genes_log2fc_values = data.frame(c(-0.69, -0.02, -0.36, -0.59))
cycle_genes_log2fc_values = t(cycle_genes_log2fc_values)
colnames(cycle_genes_log2fc_values) = c('YKL', 'THZ531', 'THZ1', 'COMBO')
cycle_genes_log2fc_values = as.data.frame(cycle_genes_log2fc_values)
rownames(cycle_genes_log2fc_values) = 'Log2FC'

cycle_genes_log2fc_values[2,] = 0
cycle_genes_log2fc_values[3,] = -0.8

cycle_genes_log2fc_values = abs(cycle_genes_log2fc_values)
cycle_genes_log2fc_values = cycle_genes_log2fc_values[c(3,2,1),]

### CYCLE GENES RADAR PLOT - YKL WT/MUT ###
par(mar=c(0.5, 0.5, 0.5, 0.5))
radarchart( cycle_genes_log2fc_values  , axistype=1 ,
            
            #custom polygon
            pcol=rgb(1, 0.2, 0.2, 0.8) , pfcol=rgb(1, 0.2, 0.2, 0.5) , plwd=4 ,
            
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1,0.2), cglwd=0.8
            
            #custom labels
            #vlcex=0.8
)
par(new = T)
cycle_genes_mut_log2fc_values = data.frame(c(-0.29, -0.21, -0.01, -0.25))
cycle_genes_mut_log2fc_values = t(cycle_genes_mut_log2fc_values)
colnames(cycle_genes_mut_log2fc_values) = c('YKL', 'THZ531', 'THZ1', 'COMBO')
cycle_genes_mut_log2fc_values = as.data.frame(cycle_genes_mut_log2fc_values)
rownames(cycle_genes_mut_log2fc_values) = 'Log2FC'

cycle_genes_mut_log2fc_values[2,] = 0
cycle_genes_mut_log2fc_values[3,] = -0.8

cycle_genes_mut_log2fc_values = abs(cycle_genes_mut_log2fc_values)  
cycle_genes_mut_log2fc_values = cycle_genes_mut_log2fc_values[c(3,2,1),]

radarchart( cycle_genes_mut_log2fc_values  , axistype=1 ,
            
            #custom polygon
            pcol=rgb(0, 0, 0, 0.3) , pfcol=rgb(0,0,0,0.5) , plwd=4 ,
            
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1,0.2), cglwd=0.8
            
            #custom labels
            #vlcex=0.8
)


legend(1.0,-0.5,
       legend=c("WT","Mut"),
       pch=c(15,15),
       col=c("red","grey"))


