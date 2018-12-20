### ANALYZING COUNTS AND GENERATING SCATTERS ###

### IMPORTING AND FILTERING DATA ###
counts_data = read.delim('/Users/DanPark1116/Desktop/BCM_Materials/Year_Two/CDK7_Analysis/CDK7_GRAY_cuffnorm/DESeq_Counts/cdk_counts.txt', stringsAsFactors = F)
counts_data$max = apply(counts_data, 1, max)
counts_data = counts_data[counts_data$max > 10, -43]

### ISOLATING RELEVANT DATA ###
ykl_dmso_counts = counts_data[, c(6,8,12,22,26,31,32,33,34,35,37,41)]
ykl_dmso_counts = ykl_dmso_counts[, c(7,11,6,2,8,3,9,10,1,4,12,5)]

### CALCULATING COUNT AVERAGES FOR DMSO AND YKL ###
ykl_dmso_counts_avgs = ykl_dmso_counts
ykl_dmso_counts_avgs$wt_dmso = rowMeans(ykl_dmso_counts_avgs[, c(1:3)])
ykl_dmso_counts_avgs$wt_ykl = rowMeans(ykl_dmso_counts_avgs[, c(4:6)])
ykl_dmso_counts_avgs$mut_dmso = rowMeans(ykl_dmso_counts_avgs[, c(7:9)])
ykl_dmso_counts_avgs$mut_ykl = rowMeans(ykl_dmso_counts_avgs[, c(10:12)])
ykl_dmso_counts_avgs = ykl_dmso_counts_avgs[, c(13:16)]

### PLOTTING COUNTS FOR DMSO TREATMENT ON WT/MUTANT CELLS ###
plot(ykl_dmso_counts_avgs$wt_dmso, ykl_dmso_counts_avgs$mut_dmso, xlim = c(0, 20000), ylim = c(0, 20000), cex = 0.5,
     xlab = 'DMSO WT Counts', ylab = 'DMSO Mutant Counts', main = 'DMSO Counts Comparison Between WT and Mutant')
abline(a = 0, b = 1, col = 'red', lty = 2)

### PLOTTING COUNTS FOR YKL TREATMENT ON WT/MUTANT CELLS ###
plot(ykl_dmso_counts_avgs$wt_ykl, ykl_dmso_counts_avgs$mut_ykl, xlim = c(0, 20000), ylim = c(0, 20000), cex = 0.5,
     xlab = 'YKL WT Counts', ylab = 'YKL Mutant Counts', main = 'YKL Counts Comparison Between WT and Mutant')
abline(a = 0, b = 1, col = 'red', lty = 2)



