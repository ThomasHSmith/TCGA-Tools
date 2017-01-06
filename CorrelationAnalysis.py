# Generate heatmaps from TCGA mRNA expression data
# to compare expression of target gene with other genes in input dataset

import sys, getopt

def main(argv):
	HELP_MSG = 'CorrelationAnalysis.py -i <df_log2_targets.pickle> -t <Target Gene (ie GAPDH)> -o <Output dir> -p <Project name>'

	try:
		opts, args = getopt.getopt(argv, 'i:t:o:p:h',['input_df=', 'target_gene=', 'output_dir=', 'project_name='])
	# Print help msg if wrong cl args
	except getopt.GetoptError:
		print HELP_MSG
		sys.exit(2)
	for o, arg in opts:
		if o == '-h':
			print HELP_MSG
			sys.exit(2)
		elif o in ('-i','--input_df'):
			PICKLE_INFILE = arg
		elif o in ('-t','--target_gene'):
			TARGET_GENE = arg
		elif o in ('-o','--output_dir'):
			OUT_PATH = arg
		elif o in ('-p','--project_name'):
			PROJECT_NAME = arg 

# may need to sort values by target gene:
# df_targets.sort_values(['SampleType', GENE_OF_INTEREST],ascending=False, inplace=True)
# df_targets_log2.sort_values(['SampleType', GENE_OF_INTEREST],ascending=False, inplace=True)
	import pandas as pd
	import seaborn as sns
	from scipy import stats
	import matplotlib.pyplot as plt
	print '\n\t***************'

	def log_echo(msg):
		print '\t' + msg
		lf.write(msg+'\n')
	
	# Build path to input/output files
	HEATMAP_OUTFILE = '%s/%s_%s_heatmap.pdf' % (OUT_PATH, PROJECT_NAME, TARGET_GENE)
	CORRELATION_OUTFILE = '%s/%s_%s_correlations.txt' %  (OUT_PATH, PROJECT_NAME, TARGET_GENE)
	ANAL_LOG_OUTFILE = '%s/%s_%s_analysis_log.txt' %  (OUT_PATH, PROJECT_NAME, TARGET_GENE)

	lf = open(ANAL_LOG_OUTFILE, 'w')
	log_echo('Writing analysis log file to: %s' % (ANAL_LOG_OUTFILE.split('/')[-1]) )
	log_echo('TCGA Project: %s' % PROJECT_NAME)
	
	df = pd.read_pickle(PICKLE_INFILE)
	log_echo('Loaded dataframe from file: %s' % PICKLE_INFILE.split('/')[-1])

	# Give Tissue type numerical value for heatmap
	# 11 - Normal tissue -> -2
	# 01 - Solid tumor -> 0
	# 02 - Recurrent solid tumor -> 1
	# 06 - Mets -> 2
	vals_dict = {'11':-2, '01':0, '02':1, '06':2}
	df['TissueType'] = df['SampleType'].apply(lambda x: vals_dict[x])

	# Drop metastatic samples
	total_n = df.shape[0]
	norm_n = df[df['TissueType'] == -2].shape[0]
	solid_n = df[df['TissueType'] == 0].shape[0]
	mets_n = df[df['TissueType'] == 2].shape[0]
	recurr_n = df[df['TissueType'] == 1].shape[0]


	log_echo('\n\tTotal cases: %i\n\tSolid Tumor:%i\n\tMetastatic:%i\n\tRecurrent Solid Tumor:%i\n\tNormal Tissue:%i\n' % (total_n, solid_n, mets_n, recurr_n, norm_n))

## FIXME ADD OPTION TO DROP SUBTYPES
#log_msg = 'Dropping metastatic samples types...'
#print log_msg
#lf.write(log_msg+'\n')
#pre_drop = df.shape[0]
#df = df[df['TissueType'] != -3]
#log_msg = 'Dropped %i metastatic cases' % (pre_drop- df.shape[0])
#print log_msg
#lf.write(log_msg+'\n')

	## FIXME Add option to include SampleType in correlation calculations?
	# Grab all non-categorical data
	df_vals = df.drop(['PtID','TumorStage','SampleType'],axis=1)
	
	df_zscores = df_vals.apply(lambda x: stats.zscore(x))
	#df_zscores['TissueType'] = df['TissueType'].values

	## FIXME Add option to inverse sorting
	df_zscores.sort_values(TARGET_GENE, ascending=False, inplace=True)

	# Compute pearson and spearman correlation coefficients for target gene and all target genes
	# Write results to tab-delimited txt file
	log_echo('Computing correlation coefficients...')

	f = open(CORRELATION_OUTFILE, 'w')
	x = list(df_zscores[TARGET_GENE])
	header = 'Gene\t PearsonR (p)\t SpearmanR (p)'
	log_echo(header)
	f.write(header + '\n')
	for column in df_zscores:
	    y = list(df_zscores[column])
	    pearsonR, pearsonP = stats.pearsonr(x, y)
	    spearmanR, spearmanP = stats.spearmanr(x, y)
	    line = '%s\t%.3f (%.3E)\t%.3f (%.3E)' % (column, pearsonR, pearsonP, spearmanR, spearmanP)
	    log_echo(line)	   
	    f.write(line+'\n')
	f.close()
	log_echo('Saved correlation data to: %s\n' % CORRELATION_OUTFILE.split('/')[-1])

## FIXME GET Z_SCORE CUTOFF FROM COMMAND LINE
## FIXME APPLY Z-CUTOFF BEFORE CALC CORRELATION COEFF
	z_cutoff = 5
	log_echo('Filtering out cases with z-score > %i' % z_cutoff)
	pre_filtered_n = len(df_zscores)
	log_echo('Number of samples (before Z-cutoff): %d' % (pre_filtered_n))
	for column in df_zscores:
	    df_zscores.drop(df_zscores[df_zscores[column] > z_cutoff].index, inplace=True)
	post_filtered_n = len(df_zscores)
	pct_lost = ( (pre_filtered_n - post_filtered_n)/float(pre_filtered_n)) * 100
	log_echo('Number of samples (after Z-cutoff): %d' % (post_filtered_n))
	log_echo('Dropped %d samples (%f%% of total dataset)\n' % ((-1*(post_filtered_n - pre_filtered_n)), pct_lost))
	
	## FIXME FIG HEIGHT MIN
	# Generate heatmap
	FIG_HEIGHT = 5
	if post_filtered_n > 70:
		FIG_HEIGHT = int(post_filtered_n/35)
	plt.figure(figsize=(25,FIG_HEIGHT))
	sns.heatmap(df_zscores,yticklabels=False, xticklabels=True)
	#plt.show()
	plt.savefig(HEATMAP_OUTFILE, bbox_inches = 'tight',dpi=300)
	log_echo('Saved heatmap as %s' % HEATMAP_OUTFILE.split('/')[-1])
	log_echo('Done!')	
	lf.close()

if __name__ == '__main__':
	main(sys.argv[1:])
