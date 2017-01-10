# Parse directory of FPKM-UQ.gz files, extract FPKM-UQ values
# Option to also extract case info from metadata JSON file
# Output as a Pandas dataframe

import pandas as pd
from tqdm import tqdm
import json, sys, getopt, os, gzip

def main(argv):
	HELP_MSG = 'USAGE:\tParseAgg.py -d <Data dir> -m <Metadata file path> -o <Output dir> -p <Project name>\n\t-d:\tSpecify path to top data directory with structure top_dir/Subdir/CaseDir/FPKM-UQ.txt.gz>\n\t-m:\tSpecify path to JSON file containing TCGA metadata\n\t-o:\tSpecify directory to save output pickle\n\t-p\tSpecify a name for the project (ie BRCA), to prepend output file'
	
	
	try:
		opts, args = getopt.getopt(argv, 'd:m:o:p:h',['top_dir=','metadata=','output_dir=','project_name='])
	# Print help msg if wrong cl args
	except getopt.GetoptError:
		print HELP_MSG
		sys.exit(2)
	for o, arg in opts:
		if o == '-h':
			print HELP_MSG
			sys.exit(2)
		elif o in ('-d','--top_dir'):
			TOP_DATA_DIR = arg
			print 'Root data directory path: %s' % TOP_DATA_DIR
		elif o in ('-m','--metadata'):
			JSON_METADATA = arg
			print 'Metadata file path: %s' % JSON_METADATA			
		elif o in ('-o','--output_dir'):
			OUTPUT_DIR = arg.rstrip('/')
			print 'Output directory path: %s' % OUTPUT_DIR	
		elif o in ('-p','--project_name'):
			PROJECT_NAME = arg
			print 'Project name: %s' % PROJECT_NAME
			PICKLE_OUTFILE = '%s/%s_aggregated_FPKM-UQ.pickle' % (OUTPUT_DIR, PROJECT_NAME)

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
	def conv_stage(x):
		stages_dict = { 'stage iia':2, 'stage iib': 2, 'stage iiia': 3, 
				'stage i': 1 , 'stage ia':1, 'stage iiic':3, 
				'stage iiib':3, 'stage iv':4, 'stage x': 0, 
				'not reported': 0, 'stage ii': 2, 'stage ib':1,
				'stage iii':3, 'na':0}
		return stages_dict[x]

	## FIXME Make this an option at command line
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
	print 'Found %d FPKM-UQ files.' % len(gz_paths)

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
	    	vals.append(df_meta[df_meta['CaseDir']==path.split('/')[-2]]['TumorStage'].get_values()[0])
	    	# FIXME
		#cols.append('TumorStageInt')
	   	#vals.append(int(df_meta[df_meta['CaseDir']==path.split('/')[2]]['TumorStageInt'].get_values()[0]))
	   	cols.append('PtID')
	    	vals.append(df_meta[df_meta['CaseDir']==path.split('/')[-2]]['PtID'].get_values()[0])
	    	cols.append('SampleType')
	    	vals.append(df_meta[df_meta['CaseDir']==path.split('/')[-2]]['SampleType'].get_values()[0])
	    	df_temp = pd.DataFrame(data=[vals], columns=cols)
	    	df = df.append(df_temp)
	    	f.close()
	
	filename = PICKLE_OUTFILE
	dir = os.path.dirname(filename)
	if not os.path.exists(dir):
    		os.makedirs(dir)
	df.to_pickle(PICKLE_OUTFILE)
	print 'Done! Wrote aggregated counts to %s.' % PICKLE_OUTFILE

if __name__ == '__main__':
	main(sys.argv[1:])

