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


	""" Sample Info """	
	# Output sample types (ie Normal tissue, solid tumor, metastatic)
	msg = ''
	for item in df.SampleType.unique():
	    msg+='%s: %d' % (item, sum(df.SampleType==item))
	msg+='Total: %d' % len(df.SampleType)
	log_echo(msg)

	## FIXME Add opt to drop subtypes before pearsons calc
	#pre_drop = df.shape[0]
	#df = df[df['TissueType'] != -3]
	#log_msg = 'Dropped %i metastatic cases' % (pre_drop- df.shape[0])
	#log_echo(log_msg)
	
	## FIXME Add option to include SampleType in correlation calculations?
	# Grab all non-categorical data


	""" Correlation calculations """
	# Compute pearson and spearman correlation coefficients for target gene and all target genes
	# Write results to tab-delimited txt file

	df_corr = df.select_dtypes(exclude=['object'])
	target_col = df_corr[TARGET_GENE]
	df_corr.drop(labels=[TARGET_GENE], axis=1, inplace=True)
	df_corr.insert(0, TARGET_GENE, target_col)

	## FIXME Make this a function, then add opt to choose whether to calc correlations	
	log_echo('Computing correlation coefficients...')
	f = open(CORRELATION_OUTFILE, 'w')
	x = list(df_corr[TARGET_GENE])
	header = 'Gene\t PearsonR (p)\t SpearmanR (p)'
	log_echo(header)
	f.write(header + '\n')
	for column in df_corr:
	    y = list(df_corr[column])
	    pearsonR, pearsonP = stats.pearsonr(x, y)
	    spearmanR, spearmanP = stats.spearmanr(x, y)
	    line = '%s\t%.3f (%.3E)\t%.3f (%.3E)' % (column, pearsonR, pearsonP, spearmanR, spearmanP)
	    log_echo(line)	   
	    f.write(line+'\n')
	f.close()
	log_echo('Saved correlation data to: %s\n' % CORRELATION_OUTFILE.split('/')[-1])

	""" Calculate z-scores """
	# Create new df for z-scores
	z_cutoff = 5 ## FIXME DEFAULT Z CUTOFF
	df_z = df.select_dtypes(exclude=['object']).apply(stats.zscore)	
	SampleToInt = {'NormalControl':-3, 'SolidTumor':0, 'Metastatic':0}
	df_z.insert(len(df_z.columns), 'SampleType', df['SampleType'].apply(lambda x: SampleToInt[x]))
	target_col = df_z[TARGET_GENE]
	df_z.drop(labels=[TARGET_GENE], axis=1, inplace=True)
	df_z.insert(0, TARGET_GENE, target_col)
	df_z.sort_values(TARGET_GENE, ascending=False, inplace=True)
	## FIXME Add option to inverse sorting
	
	# Filter out samples with z-scores above or below threshold
	log_echo('Filtering out cases with (z-score > %i) or (z-score < -%i)' % (z_cutoff, z_cutoff))
	pre_filtered_n = len(df_z)
	log_echo('Number of samples before z-cutoff: %d' % (pre_filtered_n))
	df_z = df_z[~(df_z > 5).any(axis=1)]
	df_z = df_z[~(df_z < -5).any(axis=1)]
	post_filtered_n = len(df_z)
	pct_lost = ( (pre_filtered_n - post_filtered_n)/float(pre_filtered_n)) * 100
	log_echo('Number of samples after z-cutoff: %d' % (post_filtered_n))
	log_echo('Dropped %d samples (%f%% of total dataset)\n' % ((-1*(post_filtered_n - pre_filtered_n)), pct_lost))
	
	""" Generate heatmap """
	import seaborn as sns
	import matplotlib.pyplot as plt

	FIG_HEIGHT = 5
	if post_filtered_n > 70:
		FIG_HEIGHT = int(post_filtered_n/45)
	plt.figure(figsize=(25,FIG_HEIGHT))
	sns.heatmap(df_z,yticklabels=False, xticklabels=True)
	#plt.show()
	plt.savefig(HEATMAP_OUTFILE, bbox_inches = 'tight',dpi=300)
	log_echo('Saved heatmap as %s' % HEATMAP_OUTFILE.split('/')[-1])
	log_echo('Done!')	
	lf.close()

if __name__ == '__main__':
	main(sys.argv[1:])
