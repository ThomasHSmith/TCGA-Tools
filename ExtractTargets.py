# Author: Thomas H Smith
# Given FPKM-UQ df .pickle file and list of target genes in .txt
# Return a new df pickle containing FPKM-UQ for targets and metadata 

import os, sys, getopt
import pandas as pd
import numpy as np

def main(argv):
	HELP_MSG = 'ExtractTargets.py -i <ValuesDataFrame.pickle> -t <TargetsList.txt> -o <Output dir> -p <Project name>'

	try:
		opts, args = getopt.getopt(argv, 'i:t:o:p:h',['input=', 'targets=', 'output_dir=', 'project_name='])
	# Print help msg if wrong cl args
	except getopt.GetoptError:
		print HELP_MSG
		sys.exit(2)
	for o, arg in opts:
		if o == '-h':
			print HELP_MSG
			sys.exit(2)
		elif o in ('-i','--input'):
			PICKLE_INFILE = arg
		elif o in ('-t','--targets'):
			TARGETS_INFILE = arg
		elif o in ('-o','--output_dir'):
			OUT_PATH = arg
		elif o in ('-p','--project_name'):
			PROJECT_NAME = arg 


	# Gene targets hit list .txt file should be in format: Symbol ENSG
	# ie: F2RL3	  ENSG00000127533
	#     GAPDH	  ENSG00000111640
	#     WWTR1 (TAZ) ENSG00000018408
	# ENSG ID's will be replaced with gene symbols in output

	df = pd.read_pickle(PICKLE_INFILE)
	print 'Loaded dataframe from: %s' % PICKLE_INFILE

	# Retrieve target list with genes and corresponding Ensembl IDs from targets file and store in dict
	targets_dict = {}
	f = open(TARGETS_INFILE, 'r')
	print 'Reading gene targets from: %s' % TARGETS_INFILE
	for line in f:
		words = line.split()
	 	key =' '.join(words[:-1])
	    	ID = words[-1]
	    	targets_dict[key] = ID
	f.close()

	EXCEL_OUTFILE = '%s/%s_targets.xlsx' % (OUT_PATH, PROJECT_NAME)
	PICKLE_OUTFILE = '%s/df_%s_targets.pickle' % (OUT_PATH, PROJECT_NAME)
	EXCEL_LOG2_OUTFILE = '%s/%s_targets_log2.xlsx' % (OUT_PATH, PROJECT_NAME)
	PICKLE_LOG2_OUTFILE = '%s/df_%s_targets_log2.pickle' % (OUT_PATH, PROJECT_NAME)

	# Create new df and populate with target genes
	df_targets = pd.DataFrame()
	print 'Populating new dataframe...'
	for key in targets_dict:
        	try:
			df_targets.insert(0, key, df[targets_dict[key]])
		except KeyError:
			print '%s (%s) not found in Dataset' % (targets_dict[key], key)

	# Create df containing log2-transformed values
	df_targets_log2 = df_targets.apply(lambda x: np.log2(x+1))
	
	# Populate metadata
	df_targets.insert( len(df_targets.columns), 'PtID', df['PtID'])
	df_targets.insert( len(df_targets.columns), 'TumorStage', df['TumorStage'])
	df_targets.insert( len(df_targets.columns), 'SampleType', df['SampleType'])

	df_targets_log2.insert( len(df_targets_log2.columns), 'PtID', df['PtID'])
	df_targets_log2.insert( len(df_targets_log2.columns), 'TumorStage', df['TumorStage'])
	df_targets_log2.insert( len(df_targets_log2.columns), 'SampleType', df['SampleType'])

	filename = EXCEL_OUTFILE
	dir = os.path.dirname(filename)
	if not os.path.exists(dir):
		os.makedirs(dir)
	df_targets.to_excel(EXCEL_OUTFILE)
	df_targets.to_pickle(PICKLE_OUTFILE)
	df_targets_log2.to_excel(EXCEL_LOG2_OUTFILE)
	df_targets_log2.to_pickle(PICKLE_LOG2_OUTFILE)
	print 'Output files were saved to output directory: %s\nDone!' % OUT_PATH

if __name__ == '__main__':
	main(sys.argv[1:])
