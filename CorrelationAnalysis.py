# Generate heatmaps from TCGA mRNA expression data
# to compare expression of DAB2 with target genes

import pandas as pd
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
import sys


## FIXME ADD COMMAND LINE ARG SUPPORT
## ie -i input_df -t target_gene -o output directory -p project_name

# may need to sort values by target gene:
# df_targets.sort_values(['SampleType', GENE_OF_INTEREST],ascending=False, inplace=True)
# df_targets_log2.sort_values(['SampleType', GENE_OF_INTEREST],ascending=False, inplace=True)



# Get project name as first command line argument
PROJECT = str(sys.argv[1]) 

# Build path to input/output files
PROJECT_PATH = 'TCGA_Analysis/%s' % PROJECT
PICKLE_INFILE = '%s/df_%s_DAB2_targets_log2.pickle' % (PROJECT_PATH, PROJECT)

HEATMAP_OUTFILE = '%s/%s_DAB2_heatmap.pdf' % (PROJECT_PATH, PROJECT)
CORRELATION_OUTFILE = '%s/%s_DAB2_correlations.txt' % (PROJECT_PATH, PROJECT)
ANAL_LOG_OUTFILE = '%s/%s_DAB2_analysis_log.txt' % (PROJECT_PATH, PROJECT)

lf = open(ANAL_LOG_OUTFILE, 'w')
print 'Writing analysis log file to: %s' % ANAL_LOG_OUTFILE

log_msg = 'TCGA Project: %s' % PROJECT
print log_msg
lf.write(log_msg+'\n')

df = pd.read_pickle(PICKLE_INFILE)
log_msg = 'Loaded dataframe from file: %s' % PICKLE_INFILE
print log_msg
lf.write(log_msg+'\n')

# 11 - Normal tissue -> -2
# 01 - Solid tumor -> 0
# 02 - Recurrent solid tumor -> 1
# 06 - Mets -> 2

def convert_type(x):
	x = int(x)
    	if x == 11: # If sample type is 11 (normal tissue) convert to -2
        	return -2
    	if x == 1: # If sample type is 01 (solid tumor tissue) convert to 0
        	return 0
    	if x == 2: # If sample type is 02 (recurrent solid tumor) convert to 1
        	return 1
    	if x == 6: # If sample type is 06 (tumor mets) convert to 2
        	return 2


df['TissueType'] = df['SampleType'].apply(lambda x: convert_type(x))

# Drop metastatic samples
total_n = df.shape[0]
norm_n = df[df['TissueType'] == -2].shape[0]
solid_n = df[df['TissueType'] == 0].shape[0]
mets_n = df[df['TissueType'] == 2].shape[0]
recurr_n = df[df['TissueType'] == 1].shape[0]


log_msg = 'Total cases: %i\n\tSolid Tumor:%i\n\tMetastatic:%i\n\tRecurrent Solid Tumor:%i\n\tNormal Tissue:%i' % (total_n, solid_n, mets_n, recurr_n, norm_n)
print log_msg
lf.write(log_msg+'\n')

#log_msg = 'Dropping metastatic samples types...'
#print log_msg
#lf.write(log_msg+'\n')
#pre_drop = df.shape[0]
#df = df[df['TissueType'] != -3]
#log_msg = 'Dropped %i metastatic cases' % (pre_drop- df.shape[0])
#print log_msg
#lf.write(log_msg+'\n')

df_vals = df.drop(['PtID','TumorStage','SampleType'],axis=1)

df_zscores = pd.DataFrame()
for column in df_vals:
	df_zscores[column] = stats.zscore(df_vals[column])
#df_zscores['TissueType'] = df['TissueType'].values
df_zscores.sort_values('DAB2', ascending=False, inplace=True)

# Compute pearson and spearman correlation coefficients for DAB2 and all target genes
# Write results to tab-delimited txt file
log_msg = 'Computing correlation coefficients...'
print log_msg
lf.write(log_msg+'\n')

f = open(CORRELATION_OUTFILE, 'w')
x = list(df_zscores.DAB2)
header = 'Gene\t PearsonR (p)\t SpearmanR (p)\n'
f.write(header)
for column in df_zscores:
    y = list(df_zscores[column])
    pearsonR, pearsonP = stats.pearsonr(x, y)
    spearmanR, spearmanP = stats.spearmanr(x, y)
    line = '%s\t%.3f (%.3E)\t%.3f (%.3E)' % (column, pearsonR, pearsonP, spearmanR, spearmanP)
    print line
    f.write(line+'\n')
f.close()
log_msg = 'Saved correlation data to: %s' % CORRELATION_OUTFILE
print log_msg
lf.write(log_msg+'\n')

z_cutoff_high = 5
log_msg = 'Filtering out cases with z-score > %i' % z_cutoff_high
print log_msg
lf.write(log_msg+'\n')

pre_filtered_n = len(df_zscores)
log_msg = 'Number of samples (before Z-cutoff): %d' % (pre_filtered_n)
print log_msg
lf.write(log_msg+'\n')

for column in df_zscores:
    df_zscores.drop(df_zscores[df_zscores[column] > z_cutoff_high].index, inplace=True)
    
post_filtered_n = len(df_zscores)
pct_lost = ( (pre_filtered_n - post_filtered_n)/float(pre_filtered_n)) * 100
log_msg = 'Number of samples (after Z-cutoff): %d' % (post_filtered_n)
print log_msg
lf.write(log_msg+'\n')

log_msg = 'Dropped %d samples (%f%% of total dataset)' % ((-1*(post_filtered_n - pre_filtered_n)), pct_lost)
print log_msg
lf.write(log_msg+'\n')

FIG_HEIGHT = int(post_filtered_n/35)
plt.figure(figsize=(25,FIG_HEIGHT))
sns.heatmap(df_zscores,yticklabels=False, xticklabels=True)
#plt.show()
plt.savefig(HEATMAP_OUTFILE, bbox_inches = 'tight',dpi=300)
log_msg = 'Saved heatmap as %s' % HEATMAP_OUTFILE
print log_msg
lf.write(log_msg+'\n')

print 'DONE!'
lf.write('DONE!')
lf.close()
