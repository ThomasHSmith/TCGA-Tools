# Parse directory of FPKM-UQ.gz files, extract FPKM-UQ values
# Also extract case info from metadata JSON file
# Write a Pandas dataframe

import json
import pandas as pd
from tqdm import tqdm
import os
import gzip

# Specify metadata JSON file
JSON_METADATA = 'TCGA_HNSC_metadata.json'
# Specify data directory with structure TOP_DATA_DIR/Subdir/CaseDir/FPKM-UQ.txt.gz
TOP_DATA_DIR = 'data_dir'

# Extract relevent info from JSON metadata file
# Sample types: 01=Primary Solid; 06=Mets; 11=SolidNormalTissue
json_data=open(JSON_METADATA)
data = json.load(json_data)

case_dirs, file_names, pt_ids, sample_types, tumor_stages = [], [], [], [], []
for i in range(len(data)):
	case_dirs.append(str(data[i]['file_id']))
	file_names.append(str(data[i]['file_name']))
	submitter_id = str(data[i]['cases'][0]['samples'][0]['submitter_id']).split('-')
	pt_ids.append(submitter_id[2])
	sample_types.append(submitter_id[3][:2])
	try:
        	tumor_stages.append(str(data[i]['cases'][0]['diagnoses'][0]['tumor_stage']))
    	except KeyError:
        	tumor_stages.append('na')

df_meta = pd.DataFrame({'PtID':pt_ids,'CaseDir':case_dirs,'FileName':file_names, 'SampleType':sample_types, 'TumorStage':tumor_stages})

# Helper function to formalize TumorStage field
stages_dict = { 'stage iia':2, 'stage iib': 2, 'stage iiia': 3, 'stage i': 1 , 'stage ia':1, 'stage iiic':3, 'stage iiib':3, 'stage iv':4, 'stage x': 0, 'not reported': 0, 'stage ii': 2, 'stage ib':1, 'stage iii':3, 'na':0}
def conv_stage(x):
	return stages_dict[x]

# Add numerical tumor stage column to metadata df
#df_meta['TumorStageInt'] = df_meta['TumorStage'].apply(conv_stage)

# Identify actual FPKM files (they all end with .gz extension)
# First get list of all files in directory and all subdirectories
top_dir = TOP_DATA_DIR

file_list = []
for root, dirs, files in os.walk(top_dir):
	for name in files:
        	file_list.append(os.path.join(root, name))

# Now generate new list of paths for files only ending in '.gz'
gz_paths = []
for item in file_list:
	extension = item.split('\\')[-1].split('.')[-1]
    	if extension == 'gz':
        	gz_paths.append(item)

# Parse FPKM-UQ files for each case
df = pd.DataFrame()
pbar1 = tqdm(gz_paths, total=len(gz_paths))
for path in pbar1:    
	f = gzip.open(path, 'rb')
    	cols = []
    	vals = []
    	for line in f:
        	words = line.split()
        	ensemblID = words[0].split('.')[0]
        	FPKM_value = float(words[1])
        	cols.append(ensemblID)
        	vals.append(FPKM_value)
	# Append associated metadata
    	cols.append('TumorStage')
    	vals.append(df_meta[df_meta['CaseDir']==path.split('/')[2]]['TumorStage'].get_values()[0])
    	#cols.append('TumorStageInt')
   	#vals.append(int(df_meta[df_meta['CaseDir']==path.split('/')[2]]['TumorStageInt'].get_values()[0]))
   	cols.append('PtID')
    	vals.append(df_meta[df_meta['CaseDir']==path.split('/')[2]]['PtID'].get_values()[0])
    	cols.append('SampleType')
    	vals.append(df_meta[df_meta['CaseDir']==path.split('/')[2]]['SampleType'].get_values()[0])
    	df_temp = pd.DataFrame(data=[vals], columns=cols)
    	df = df.append(df_temp)
    	f.close()
df.to_pickle('df_FPKM-UQ_with_MetaData')
print 'Done! Wrote FPKM-UQ_with_MetaData.pickle'

