### INDIVIDUAL GENE RADAR PLOTS ###

### INSTALLING PACKAGES FOR RADAR PLOTS ###
install.packages('fmsb')
library(fmsb)

### IMPORTING DATA ###
gene_data = read.csv('/Users/DanPark1116/Desktop/BCM_Materials/Year_Two/CDK7_Analysis/Initial_Output/cdk7_genes_fpkm.csv', stringsAsFactors = F)
gene_data = gene_data[!duplicated(gene_data$tracking_id),] #removes 2
rownames(gene_data) = gene_data[,1]
gene_data = gene_data[, -1]

### FILTERING FOR EXPRESSED GENES ###
gene_data$max = apply(gene_data, 1, max)
gene_data = gene_data[gene_data$max > 10, -43]

gene_data_mut = gene_data[, c(22:42)]

### GENERATING GROUPS FOR EACH TREATMENT AND DMSO ###
ykl_dmso = gene_data[, c(4:6, 16:18)]
ykl_high_dmso = gene_data[, c(1:3, 16:18)]
thz531_dmso = gene_data[, c(7:9, 16:18)]
combo_dmso = gene_data[, c(13:15, 16:18)]
combo_high_dmso = gene_data[, c(10:12, 16:18)]
thz1_dmso = gene_data[, c(19:21, 16:18)]

ykl_mut_dmso = gene_data_mut[, c(4:6, 16:18)]
thz531_mut_dmso = gene_data_mut[, c(7:9, 16:18)]
combo_mut_dmso = gene_data_mut[, c(13:15, 16:18)]
thz1_mut_dmso = gene_data_mut[, c(19:21, 16:18)]

### CALCULATING LOG2FC FOR ALL TREATMENTS AND ALL GENES ###
ykl_dmso$dmso_avg = rowMeans(ykl_dmso[, c(4:6)])
ykl_dmso$ykl_avg = rowMeans(ykl_dmso[, c(1:3)])
ykl_dmso$log2fc = log2(ykl_dmso$ykl_avg / ykl_dmso$dmso_avg)

ykl_mut_dmso$dmso_avg = rowMeans(ykl_mut_dmso[, c(4:6)])
ykl_mut_dmso$ykl_avg = rowMeans(ykl_mut_dmso[, c(1:3)])
ykl_mut_dmso$log2fc = log2(ykl_mut_dmso$ykl_avg / ykl_mut_dmso$dmso_avg)

thz531_dmso$dmso_avg = rowMeans(thz531_dmso[, c(4:6)])
thz531_dmso$thz531_avg = rowMeans(thz531_dmso[, c(1:3)])
thz531_dmso$log2fc = log2(thz531_dmso$thz531_avg / thz531_dmso$dmso_avg)

thz531_mut_dmso$dmso_avg = rowMeans(thz531_mut_dmso[, c(4:6)])
thz531_mut_dmso$thz531_avg = rowMeans(thz531_mut_dmso[, c(1:3)])
thz531_mut_dmso$log2fc = log2(thz531_mut_dmso$thz531_avg / thz531_mut_dmso$dmso_avg)

thz1_dmso$dmso_avg = rowMeans(thz1_dmso[, c(4:6)])
thz1_dmso$thz1_avg = rowMeans(thz1_dmso[, c(1:3)])
thz1_dmso$log2fc = log2(thz1_dmso$thz1_avg / thz1_dmso$dmso_avg)

thz1_mut_dmso$dmso_avg = rowMeans(thz1_mut_dmso[, c(4:6)])
thz1_mut_dmso$thz1_avg = rowMeans(thz1_mut_dmso[, c(1:3)])
thz1_mut_dmso$log2fc = log2(thz1_mut_dmso$thz1_avg / thz1_mut_dmso$dmso_avg)

combo_dmso$dmso_avg = rowMeans(combo_dmso[, c(4:6)])
combo_dmso$combo_avg = rowMeans(combo_dmso[, c(1:3)])
combo_dmso$log2fc = log2(combo_dmso$combo_avg / combo_dmso$dmso_avg)

combo_mut_dmso$dmso_avg = rowMeans(combo_mut_dmso[, c(4:6)])
combo_mut_dmso$combo_avg = rowMeans(combo_mut_dmso[, c(1:3)])
combo_mut_dmso$log2fc = log2(combo_mut_dmso$combo_avg / combo_mut_dmso$dmso_avg)

### GENERATING LOG2FC RADAR PLOTS FOR INDIVIDUAL GENES - PSMC3 ###
psmc3_ykl_wt = ykl_dmso[rownames(ykl_dmso) == 'PSMC3',] #-0.11
psmc3_ykl_mut = ykl_mut_dmso[rownames(ykl_mut_dmso) == 'PSMC3',] #-0.15

psmc3_thz531_wt = thz531_dmso[rownames(thz531_dmso) == 'PSMC3',] #-0.38
psmc3_thz531_mut = thz531_mut_dmso[rownames(thz531_mut_dmso) == 'PSMC3',] #-0.37

psmc3_thz1_wt = thz1_dmso[rownames(thz1_dmso) == 'PSMC3',] #-0.05
psmc3_thz1_mut = thz1_mut_dmso[rownames(thz1_mut_dmso) == 'PSMC3',] #-0.42

psmc3_combo_wt = combo_dmso[rownames(combo_dmso) == 'PSMC3',] #-0.34
pscm3_combo_mut = combo_mut_dmso[rownames(combo_mut_dmso) == 'PSMC3',] #-0.32

psmc3_wt_log2fc = data.frame(c(0.11, 0.38, 0.05, 0.34))
psmc3_wt_log2fc = t(psmc3_wt_log2fc)
colnames(psmc3_wt_log2fc) = c('YKL', 'THZ531', 'THZ1', 'COMBO')
psmc3_wt_log2fc = as.data.frame(psmc3_wt_log2fc)
rownames(psmc3_wt_log2fc) = 'Log2FC'

psmc3_wt_log2fc[2,] = 0
psmc3_wt_log2fc[3,] = 0.6

psmc3_wt_log2fc = psmc3_wt_log2fc[c(3,2,1),]

psmc3_mut_log2fc = data.frame(c(0.15, 0.37, 0.42, 0.32))
psmc3_mut_log2fc = t(psmc3_mut_log2fc)
colnames(psmc3_mut_log2fc) = c('YKL', 'THZ531', 'THZ1', 'COMBO')
psmc3_mut_log2fc = as.data.frame(psmc3_mut_log2fc)
rownames(psmc3_mut_log2fc) = 'Log2FC'

psmc3_mut_log2fc[2,] = 0
psmc3_mut_log2fc[3,] = 0.6

psmc3_mut_log2fc = psmc3_mut_log2fc[c(3,2,1),]

### CREATING RADAR PLOT FOR PSMC3 ###
par(mar=c(0.5, 0.5, 0.5, 0.5))
radarchart( psmc3_wt_log2fc  , axistype=1 ,
            
            #custom polygon
            pcol=rgb(1, 0.2, 0.2, 0.8) , pfcol=rgb(1, 0.2, 0.2, 0.5) , plwd=4 ,
            
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,0.6,0.15), cglwd=0.8
            
            #custom labels
            #vlcex=0.8
)
par(new = T)
radarchart( psmc3_mut_log2fc  , axistype=1 ,
            
            #custom polygon
            pcol=rgb(0, 0, 0, 0.3) , pfcol=rgb(0,0,0,0.5) , plwd=4 ,
            
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,0.6,0.15), cglwd=0.8
            
            #custom labels
            #vlcex=0.8
)

legend(1.0,-0.5,
       legend=c("WT","Mut"),
       pch=c(15,15),
       col=c("red","grey"))

### GENERATING LOG2FC RADAR PLOTS FOR INDIVIDUAL GENES - NGFRAP1 ###
NGFRAP1_ykl_wt = ykl_dmso[rownames(ykl_dmso) == 'NGFRAP1',] #-0.35
NGFRAP1_ykl_mut = ykl_mut_dmso[rownames(ykl_mut_dmso) == 'NGFRAP1',] #-0.31

NGFRAP1_thz531_wt = thz531_dmso[rownames(thz531_dmso) == 'NGFRAP1',] #-0.49
NGFRAP1_thz531_mut = thz531_mut_dmso[rownames(thz531_mut_dmso) == 'NGFRAP1',] #-0.42

NGFRAP1_thz1_wt = thz1_dmso[rownames(thz1_dmso) == 'NGFRAP1',] #-0.22
NGFRAP1_thz1_mut = thz1_mut_dmso[rownames(thz1_mut_dmso) == 'NGFRAP1',] #-0.22

NGFRAP1_combo_wt = combo_dmso[rownames(combo_dmso) == 'NGFRAP1',] #-0.50
NGFRAP1_combo_mut = combo_mut_dmso[rownames(combo_mut_dmso) == 'NGFRAP1',] #-0.44

NGFRAP1_wt_log2fc = data.frame(c(0.35, 0.49, 0.22, 0.50))
NGFRAP1_wt_log2fc = t(NGFRAP1_wt_log2fc)
colnames(NGFRAP1_wt_log2fc) = c('YKL', 'THZ531', 'THZ1', 'COMBO')
NGFRAP1_wt_log2fc = as.data.frame(NGFRAP1_wt_log2fc)
rownames(NGFRAP1_wt_log2fc) = 'Log2FC'

NGFRAP1_wt_log2fc[2,] = 0
NGFRAP1_wt_log2fc[3,] = 0.6

NGFRAP1_wt_log2fc = NGFRAP1_wt_log2fc[c(3,2,1),]

NGFRAP1_mut_log2fc = data.frame(c(0.31, 0.42, 0.22, 0.44))
NGFRAP1_mut_log2fc = t(NGFRAP1_mut_log2fc)
colnames(NGFRAP1_mut_log2fc) = c('YKL', 'THZ531', 'THZ1', 'COMBO')
NGFRAP1_mut_log2fc = as.data.frame(NGFRAP1_mut_log2fc)
rownames(NGFRAP1_mut_log2fc) = 'Log2FC'

NGFRAP1_mut_log2fc[2,] = 0
NGFRAP1_mut_log2fc[3,] = 0.6

NGFRAP1_mut_log2fc = NGFRAP1_mut_log2fc[c(3,2,1),]

### CREATING RADAR PLOT FOR NGFRAP1 ###
par(mar=c(0.5, 0.5, 0.5, 0.5))
radarchart( NGFRAP1_wt_log2fc  , axistype=1 ,
            
            #custom polygon
            pcol=rgb(1, 0.2, 0.2, 0.8) , pfcol=rgb(1, 0.2, 0.2, 0.5) , plwd=4 ,
            
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,0.6,0.15), cglwd=0.8
            
            #custom labels
            #vlcex=0.8
)
par(new = T)
radarchart( NGFRAP1_mut_log2fc  , axistype=1 ,
            
            #custom polygon
            pcol=rgb(0, 0, 0, 0.3) , pfcol=rgb(0,0,0,0.5) , plwd=4 ,
            
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,0.6,0.15), cglwd=0.8
            
            #custom labels
            #vlcex=0.8
)

legend(1.0,-0.5,
       legend=c("WT","Mut"),
       pch=c(15,15),
       col=c("red","grey"))





