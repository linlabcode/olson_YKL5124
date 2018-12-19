###python dependencies
import pandas as pd 
import seaborn as sns
from statsmodels.stats.multitest import multipletests
import numpy as np
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt

###file locations
projectDir = "/path/to/projectDir/"
cuffnorm = "{}/cufflinks_formatted/CDK7_GRAY_cuffnorm/output".format(projectDir)
ykl_dmso_wt_low = "{}/CDK7_GRAY_HAP1_WT_DMSO_R_vs_HAP1_WT_YKL_LO_R_exprs_matrix.txt".format(cuffnorm)
thz531_dmso_wt_low = "{}/CDK7_GRAY_HAP1_WT_531_R_vs_HAP1_WT_DMSO_R_exprs_matrix.txt".format(cuffnorm)
thz1_dmso_wt_low = "{}/CDK7_GRAY_HAP1_WT_DMSO_R_vs_HAP1_WT_THZ1_R_exprs_matrix.txt".format(cuffnorm)
combo_dmso_wt_low = "{}/CDK7_GRAY_HAP1_WT_COMBO_LO_R_vs_HAP1_WT_DMSO_R_exprs_matrix.txt".format(cuffnorm)
fpkmTable = "{}/CDK7_GRAY_all_fpkm_exprs_raw.txt".format(cuffnorm)


def filterGenes(dataframe): 
    gene_list = dataframe.index.tolist()
    gene_list = [gene for gene in gene_list if 'HIST' not in gene and 'MIR' not in gene and 'SNO' not in gene]
    dataframe_filtered = dataframe.loc[gene_list, :]
    return dataframe_filtered

def getDownGenes(dataframe): 
    dataframe_de = dataframe[(dataframe.LOG2_FOLD_CHANGE < -0.5) & (dataframe.P_VALUE < 0.05)]
    downGenes = dataframe_de.index.tolist()
    return downGenes

def getJaccardAndFingerPrints(conditionA, conditionB): 
    '''
    Function to compute the jaccard similarity between two sets of genes and generate binary vectors for 
    statistical evaluation
    '''
    conditionA = set(conditionA)
    conditionB = set(conditionB)
    union=conditionA.union(conditionB)
    intersection = conditionA.intersection(conditionB)
    jaccardSimilarity = len(intersection)/len(union)
    condA_notB = conditionA.difference(conditionB)
    condB_notA = conditionB.difference(conditionA)
    ###computer fischer's exact test 
    contingency_table = [[len(union), len(condA_notB)], [len(condB_notA), len(intersection)]]
    oddratio, pval = fisher_exact(contingency_table, alternative='greater')
    #fingerprintA = np.zeros(len(union))
    #fingerprintB = np.zeros(len(union))
    #for idx,gene in enumerate(list(union)): 
    #    if gene in conditionA: fingerprintA[idx] = 1
    #    if gene in conditionB: fingerprintB[idx] = 1
    
    return jaccardSimilarity, pval

###load up dataframes
ykl_dmso_wt_low_df = pd.read_csv(ykl_dmso_wt_low, sep="\t")
thz531_dmso_wt_low_df = pd.read_csv(thz531_dmso_wt_low, sep="\t")
thz1_dmso_wt_low_df = pd.read_csv(thz1_dmso_wt_low, sep="\t")
combo_dmso_wt_low_df = pd.read_csv(combo_dmso_wt_low, sep="\t")

###correcting log2 fold change based on order analysis was done 
thz531_dmso_wt_low_df['LOG2_FOLD_CHANGE'] = thz531_dmso_wt_low_df.LOG2_FOLD_CHANGE*-1
combo_dmso_wt_low_df['LOG2_FOLD_CHANGE'] = combo_dmso_wt_low_df.LOG2_FOLD_CHANGE*-1

dfConditions = ['Combo', 'THZ1', 'THZ531', 'YKL-5-124']
dfList = [combo_dmso_wt_low_df,thz1_dmso_wt_low_df,thz531_dmso_wt_low_df,ykl_dmso_wt_low_df,]

##filter genes
dfFilteredList = [filterGenes(dataframe) for dataframe in dfList]


####obtaining list of down-regulated genes for treatment vs wt 

downGeneLists = [getDownGenes(dataframe) for dataframe in dfAdjList]

##jaccard table matrix
jaccardMatrix = np.zeros((len(downGeneLists), len(downGeneLists)))
jaccardPValMatrix = np.zeros((len(downGeneLists), len(downGeneLists)))
for row,conditionA in enumerate(downGeneLists): 
    for column,conditionB in enumerate(downGeneLists):
        #jaccard, fingerprintA, fingerprintB = getJaccardAndFingerPrints(conditionA, conditionB)
        jaccard, pval = getJaccardAndFingerPrints(conditionA, conditionB)
        jaccardMatrix[row,column] = jaccard
        jaccardPValMatrix[row,column] = pval


##jaccard figure 
jaccDF = pd.DataFrame(jaccardMatrix, columns=dfConditions, index=dfConditions)
f, ax = plt.subplots()
sns.heatmap(jaccDF, annot = True, linewidths = 0.5, ax=ax, cmap="Blues")
ax.set_ylabel('')
plt.tight_layout()
f.savefig("{}/figures/wt_low_jaccard_annotHeatmap_log2fc_less.pdf".format(projectDir))