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

	import pandas as pd
	from scipy import stats
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
	log_echo('Gene of interest: %s' % TARGET_GENE)
	df = pd.read_pickle(PICKLE_INFILE)
	log_echo('Loaded dataframe from file: %s' % PICKLE_INFILE.split('/')[-1])

	vals_dict = {'11':-2, '01':0, '02':1, '06':2}
	df['TissueType'] = df['SampleType'].apply(lambda x: vals_dict[x])

	total_n = df.shape[0]
	norm_n = df[df['TissueType'] == -2].shape[0]
	solid_n = df[df['TissueType'] == 0].shape[0]
	mets_n = df[df['TissueType'] == 2].shape[0]
	recurr_n = df[df['TissueType'] == 1].shape[0]


	log_echo('\n\tTotal cases: %i\n\tSolid Tumor:%i\n\tMetastatic:%i\n\tRecurrent Solid Tumor:%i\n\tNormal Tissue:%i\n' % (total_n, solid_n, mets_n, recurr_n, norm_n))

	
	## FIXME ADD OPTION TO DROP SUBTYPES

	#pre_drop = df.shape[0]
	#df = df[df['TissueType'] != -3]
	#log_msg = 'Dropped %i metastatic cases' % (pre_drop- df.shape[0])
	#log_echo(log_msg)

	## FIXME Add option to include SampleType in correlation calculations?
	# Grab all non-categorical data
	df_vals = df.drop(['PtID','TumorStage','SampleType'],axis=1)
	
	## FIXME Add option for zero-filtering
#	L = []
#	zero_cutoff = 50
#	for column in df_vals:
#    		if (df_vals[df_vals[column]== 0].shape[0]) > zero_cutoff:
#        		L.append(column)	
#	df_vals.drop(L, inplace=True, axis=1)
#	log_echo('Dropped %d targets from list with zero values for > %d samples' % (len(L), zero_cutoff))
	
	z_cutoff = 5
	log_echo('Filtering out cases with z-score > %i' % z_cutoff)
	pre_filtered_n = len(df_zscores)
	log_echo('Number of samples (before Z-cutoff): %d' % (pre_filtered_n))
	df_zscores = pd.DataFrame()
	for column in df_vals:
   		df_zscores[column] = stats.zscore(df_vals[column])
	
	df_zscores = pd.DataFrame()
	for column in df:
    		df_zscores[column] = stats.zscore(df[column])
	df_zscores['TissueType'] = df['TissueType'].values
	for column in df_zscores:
    		df_zscores.drop(df_zscores[df_zscores[column] > z_cutoff_high].index, inplace=True)
	
	post_filtered_n = len(df_zscores)
	pct_lost = ( (pre_filtered_n - post_filtered_n)/float(pre_filtered_n)) * 100
	log_echo('Number of samples (after Z-cutoff): %d' % (post_filtered_n))
	log_echo('Dropped %d samples (%f%% of total dataset)\n' % ((-1*(post_filtered_n - pre_filtered_n)), pct_lost))
	

	# Move target gene to first position
	target_col = df_zscores[TARGET_GENE]
	df_zscores.drop(labels=[TARGET_GENE], axis=1, inplace=True)
	df_zscores.insert(0, TARGET_GENE, target_col)

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
	df_filtered = pd.DataFrame()
	for column in df_zscores:
		df_filtered[column] = df_zscores[df_zscores[column] <= z_cutoff].index
		#df_zscores.drop(df_zscores[df_zscores[column] > z_cutoff].index, inplace=True)
	#post_filtered_n = len(df_zscores)
	## FIXME FIG HEIGHT MIN
	# Generate heatmap

	import seaborn as sns
	import matplotlib.pyplot as plt

	FIG_HEIGHT = 5
	if post_filtered_n > 70:
		FIG_HEIGHT = int(post_filtered_n/45)
	plt.figure(figsize=(25,FIG_HEIGHT))
	sns.heatmap(df_zscores,yticklabels=False, xticklabels=True)
	#plt.show()
	plt.savefig(HEATMAP_OUTFILE, bbox_inches = 'tight',dpi=300)
	log_echo('Saved heatmap as %s' % HEATMAP_OUTFILE.split('/')[-1])
	log_echo('Done!')	
	lf.close()

if __name__ == '__main__':
	main(sys.argv[1:])
