# Generate heatmaps from TCGA mRNA expression data
# to compare expression of target gene with other genes in input dataset

import sys, argparse, os, datetime, errno

#def main(argv):

#TARGET_GENE, PICKLE_INFILE, OUT_PATH, PROJECT_NAME = '','','',''
now=str(datetime.datetime.now()).split()
# Get input from user cl arguments
parser = argparse.ArgumentParser(description='CorrelationAnalysis', add_help=True)

parser.add_argument('-i', '--input_df', required=True,
	help="Path to df_log2.pickle file that was "
	"generated from ExtractTargets.py "
	"Ex: '-i BRCA_hippo_df_log2.pickle'")

parser.add_argument('-t', '--target_gene', required=False,
	help="Target gene of interest! Must match name from "
	"the targets file.  Default will choose first"
	"gene in list Ex:  '-t F2RL3'", default='')
	
parser.add_argument('-o', '--output_dir', required=False,
	help="Folder to save output files to?  Ex: '-o ./output'"
	"  Default is current directory" , default='./')

parser.add_argument('-p','--project_name', required=False,
	help="Short name to add on to the output files "
	"ie '-p BRCA_F2RL3_20160606'  Default is target/date_time",
	default=(now[0] + '_' + now[1][:2]+now[1][3:5]))
		## FIXME Add ability to write several targets
		## USE add_argument(' ' action='append', dest=targets,default=[])
parser.add_argument('-z','--z_cutoff', required=False,
	help="Exclude samples with z-scores above and"
	"below a certain threshold.  Default is 5. Set to 0 "
	"for no cutoff.  Set to -1 to disable z-score conversion",
	default=5, type=int)
parser.add_argument('--no_heat', required=False,
	help="Option to disable generation of heatmap",
	action='store_true')
parser.add_argument('--no_corr', required=False,
	help="Option to disable generation of correlation data",
	action='store_true')





import pandas as pd
from scipy import stats
args = parser.parse_args()

def make_out_dir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

OUT_PATH=args.output_dir
make_out_dir(OUT_PATH)
PROJECT_NAME = args.project_name
TARGET_GENE = args.target_gene
PICKLE_INFILE = args.input_df
Z_CUTOFF = args.z_cutoff
print '\n\t***************'
	
	
# Build path to input/output files
HEATMAP_OUTFILE = '%s/%s_%s_heatmap.pdf' % (OUT_PATH, PROJECT_NAME, TARGET_GENE)
CORRELATION_OUTFILE = '%s/%s_%s_correlations.txt' %  (OUT_PATH, PROJECT_NAME, TARGET_GENE)
ANAL_LOG_OUTFILE = '%s/%s_%s_analysis_log.txt' %  (OUT_PATH,PROJECT_NAME, TARGET_GENE)

lf = open(ANAL_LOG_OUTFILE, 'w')
def log_echo(msg):
	print '\t' + msg
	lf.write(msg+'\n')

log_echo('Writing analysis log file to: %s' % (ANAL_LOG_OUTFILE.split('/')[-1]) )
log_echo('TCGA Project: %s' % PROJECT_NAME)
log_echo('Gene of interest: %s' % TARGET_GENE)
df = pd.read_pickle(PICKLE_INFILE)
log_echo('Loaded dataframe from file: %s\n' % PICKLE_INFILE.split('/')[-1])


""" Sample Info """	
# Output sample types (ie Normal tissue, solid tumor, metastatic)

for item in df.SampleType.unique():
	msg='%s: %d' % (item, sum(df.SampleType==item))
	log_echo(msg)
msg='Total: %d\n' % len(df.SampleType)
log_echo(msg)

	## FIXME Add opt to drop subtypes before pearsons calc
	#pre_drop = df.shape[0]
	#df = df[df['TissueType'] != -3]
	#log_msg = 'Dropped %i metastatic cases' % (pre_drop- df.shape[0])
	#log_echo(log_msg)
	
	## FIXME Add option to include SampleType in correlation calculations?
	# Grab all non-categorical data


def ComputeCorrelation():
	""" Correlation calculations """
	# Compute pearson and spearman correlation coefficients for target gene and all target genes
	# Write results to tab-delimited txt file

	df_corr2 = df.select_dtypes(exclude=['object'])
	if TARGET_GENE == '':
		TARGET_GENE = df_corr2.columns[0]
	target_col = df_corr2[TARGET_GENE]
	df_corr = df_corr2.drop(labels=[TARGET_GENE], axis=1)
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

if args.no_corr == False:
	ComputeCorrelation()


""" Calculate z-scores """
# Create new df for z-scores
if Z_CUTOFF >= 0:
	df_z = df.select_dtypes(exclude=['object']).apply(stats.zscore)
elif Z_CUTOFF == -1:  # Option to not plot z-scores
	df_z = df.select_dtypes(exclude=['object'])	

SampleToInt = {'NormalControl':-3, 'SolidTumor':0, 'Metastatic':0}
df_z.insert(len(df_z.columns), 'SampleType', df['SampleType'].apply(lambda x: SampleToInt[x]))
target_col = df_z[TARGET_GENE]
df_z.drop(labels=[TARGET_GENE], axis=1, inplace=True)
df_z.insert(0, TARGET_GENE, target_col)
df_z.sort_values(TARGET_GENE, ascending=False, inplace=True)
	## FIXME Add option to inverse sorting

if Z_CUTOFF > 0: # Filter samlpes if a z-score cutoff was provided
	log_echo('Filtering out cases with (z-score > %i) or (z-score < -%i)' % (Z_CUTOFF, Z_CUTOFF))
	pre_filtered_n = len(df_z)
	log_echo('Number of samples before z-cutoff: %d' % (pre_filtered_n))
	df_z = df_z[~(df_z > Z_CUTOFF).any(axis=1)]
	df_z = df_z[~(df_z < -Z_CUTOFF).any(axis=1)]
	post_filtered_n = len(df_z)
	pct_lost = ( (pre_filtered_n - post_filtered_n)/float(pre_filtered_n)) * 100
	log_echo('Number of samples after z-cutoff: %d' % (post_filtered_n))
	log_echo('Dropped %d samples (%f%% of total dataset)\n' % ((-1*(post_filtered_n - pre_filtered_n)), pct_lost))

def GenerateHeatMap():
	FIG_HEIGHT = 5
	if post_filtered_n > 70:
		FIG_HEIGHT = int(post_filtered_n/45)
	plt.figure(figsize=(25,FIG_HEIGHT))
	sns.heatmap(df_z,yticklabels=False, xticklabels=True)
	#plt.show()
	plt.savefig(HEATMAP_OUTFILE, bbox_inches = 'tight',dpi=300)
	log_echo('Saved heatmap as %s' % HEATMAP_OUTFILE.split('/')[-1])

if args.no_heat == False:	
	import seaborn as sns
	import matplotlib.pyplot as plt
	GenerateHeatMap()

log_echo('Done!')	
lf.close()

