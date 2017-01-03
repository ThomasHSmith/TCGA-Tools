# Author: Thomas H Smith
# Given FPKM-UQ df .pickle file and list of target genes in .txt
# Return a new df pickle containing FPKM-UQ for targets and metadata 

import os
import gzip
import pandas as pd
import numpy as np
import sys

# Specify project name - must match Data subdirectory
# containing FPKM-UQ_with_MetaData pickle file
PROJECT = str(sys.argv[1]) 

# Specifc gene hit list .txt file in format: Symbol ENSG
# First gene on list will be considered gene of interest for correlation calculations
# ie: F2RL3	ENSG00000127533
# Should be located in Analysis/TargetLists/
TARGETS_FILENAME = ''

# Directory info - should not need to be changed
BASE_PATH = '/home/tsmith/Data/TCGA'
ANAL_PATH = '/home/tsmith/Analysis'

# Build path to input/output files
PICKLE_INFILE = '%s/%s/df_FPKM-UQ_with_MetaData' % (BASE_PATH, PROJECT)
GENE_TARGETS_INFILE = '%s/TargetLists/%s' % (ANAL_PATH, TARGETS_FILENAME)
OUT_PATH = '%s/TCGA/output/%s' % (ANAL_PATH, PROJECT)

df = pd.read_pickle(PICKLE_INFILE)
print 'Loaded dataframe from %s' % PICKLE_INFILE

# Retrieve target list with genes and corresponding Ensembl IDs from file and store in dict
targets_dict = {}
f = open(GENE_TARGETS_INFILE, 'r')
print 'Reading genes from targets file: %s' % GENE_TARGETS_INFILE
FIRST_LINE = True
for line in f:
	words = line.split()
 	key =' '.join(words[:-1])
    	ID = words[-1]
    	targets_dict[key] = ID
	if FIRST_LINE:
		GENE_OF_INTEREST = key
		FIRST_LINE = False
f.close()

EXCEL_OUTFILE = '%s/%s_%s_targets.xlsx' % (OUT_PATH, PROJECT, GENE_OF_INTEREST)
PICKLE_OUTFILE = '%s/df_%s_%s_targets.pickle' % (OUT_PATH, PROJECT, GENE_OF_INTEREST)
EXCEL_OUTFILE_LOG2 = '%s/%s_%s_targets_log2.xlsx' % (OUT_PATH, PROJECT, GENE_OF_INTEREST)
PICKLE_OUTFILE_LOG2 = '%s/df_%s_%s_targets_log2.pickle' % (OUT_PATH, PROJECT, GENE_OF_INTEREST)

# Create new df and populate with target genes
df_targets = pd.DataFrame()
print 'Populating new dataframe...'
for key in targets_dict:
	if key != GENE_OF_INTEREST:
        	df_targets.insert(0, key, df[targets_dict[key]])
df_targets.insert(0, GENE_OF_INTEREST, df[targets_dict[GENE_OF_INTEREST]])

# Create df containing log2-transformed values
df_targets_log2 = df_targets.apply(lambda x: np.log2(x+1))

df_targets.insert( len(df_targets.columns), 'PtID', df['PtID'])
df_targets.insert( len(df_targets.columns), 'TumorStage', df['TumorStage'])
df_targets.insert( len(df_targets.columns), 'SampleType', df['SampleType'])

df_targets_log2.insert( len(df_targets_log2.columns), 'PtID', df['PtID'])
df_targets_log2.insert( len(df_targets_log2.columns), 'TumorStage', df['TumorStage'])
df_targets_log2.insert( len(df_targets_log2.columns), 'SampleType', df['SampleType'])

df_targets.sort_values(['SampleType', GENE_OF_INTEREST],ascending=False, inplace=True)
df_targets_log2.sort_values(['SampleType', GENE_OF_INTEREST],ascending=False, inplace=True)

df_targets.to_excel(EXCEL_OUTFILE)
df_targets.to_pickle(PICKLE_OUTFILE)
df_targets_log2.to_excel(EXCEL_OUTFILE_LOG2)
df_targets_log2.to_pickle(PICKLE_OUTFILE_LOG2)
