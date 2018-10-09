# TCGA-Tools

This is a set of simple scripts that I wrote for the purpose of parsing/aggregating, filtering, and visualizing open-access RNA-Seq gene expression data obtained from the National Cancer Institute Genomic Data Commons [NCI-GDC](https://gdc-portal.nci.nih.gov).  I have specifically used it for the TCGA project datasets.  Parsed datasets can be loaded into [TCGA-Tools-GUI](https://github.com/ThomasHSmith/TCGA-Tools-GUI) for easy analysis and visualization.

## Getting Started


### Prerequisites
Python 2.7 (haven't tested with 3).  The following packages are required: Pandas, Tqdm, numpy, scipy.stats, seaborn, matplotlib.

### Installing
Coming soon

### 1. Data structure/download
FPKM or FPKM-UQ data should be downloaded from NCI-GDC using their data transfer client <https://gdc.nci.nih.gov/access-data/gdc-data-transfer-tool>, make sure to also download the metadata JSON file.  The downloaded data should have the following file structure:

	.
	├── ExtractTargets.py
	├── ParseAgg.py
	├── README.md
	└── TestData
	    ├── BRCA_metadata.json
	    └── data_dir
	        └── data
	            ├── 0a0ba3c5-0443-4c81-83cb-7fb1ddcfd054
	            │   ├── 44ebe859-bc4f-40c1-bdb9-1fa970b680da.FPKM-UQ.txt.gz
	            │   └── logs
	            │       └── 44ebe859-bc4f-40c1-bdb9-1fa970b680da.FPKM-UQ.txt.gz.parcel
	            ├── 0a19beeb-7d12-4ecb-882b-3a9d311bd1b2
	            │   ├── fa9d8ef2-0fd0-4fd8-8266-e4760b77d103.FPKM-UQ.txt.gz
	            │   └── logs
	            │       └── fa9d8ef2-0fd0-4fd8-8266-e4760b77d103.FPKM-UQ.txt.gz.parcel
	            ├── 0a40e0e1-a97d-4643-a8b2-2907ea13812e
	            │   ├── bdd8c340-250b-474a-8802-7653b7884ced.FPKM-UQ.txt.gz
	            │   └── logs
	            │       └── bdd8c340-250b-474a-8802-7653b7884ced.FPKM-UQ.txt.gz.parcel
	            ├── 0aa46e8a-0219-4996-86db-cf021b9c3882
	            │   ├── 4cb8bb61-b1d0-44d3-beb2-c6989b60aa50.FPKM-UQ.txt.gz
	            │   └── logs
	            │       └── 4cb8bb61-b1d0-44d3-beb2-c6989b60aa50.FPKM-UQ.txt.gz.parcel
	            └── 0b5c15df-e0d6-473b-93c9-c2458193647e
	                ├── 2d21d348-762a-424f-8965-e5132ff2359c.FPKM-UQ.txt.gz
	                └── logs
	                    └── 2d21d348-762a-424f-8965-e5132ff2359c.FPKM-UQ.txt.gz.parcel
	
	13 directories, 14 files


Each sub-directory of the 'data' directory represents one case (5 cases are shown in the above example, but some TCGA projects contain several hundred or over a thousand).  Within each case, there is a .FPKM-UQ.txt.gz file, which is the only one of interest to us.  These .txt files contain a list of FPKM-UQ values paired with Ensembl gene IDs, for example:

	ENSG00000134108.11	454530.58045
	ENSG00000263089.1	3083.42506523
	ENSG00000172137.17	32261.3825595
	ENSG00000167700.7	420363.026919
### 2. Parse FPKM-UQ Data with ParseAgg.py
ParseAgg.py will search through a specified top-level directory ('TestData' in the example above), identify and unzip all FPKM-UQ.txt.gz files, extract relevent data, and aggregate into a Pandas DataFrame.  
ParseAgg.py requires the following command-line arguments:

	-d <path to top level data directory>
	-m <path to metadata JSON file>
	-o <path to output directory>
	-p <project name (ie BRCA)>

Thus, to parse the data in the above example:

	$ python ParseAgg.py -d /TestData -m TestData/BRCA_metadata.json -o TestData/Output -p BRCA

The output will be a pickle containing the Pandas DataFrame with the aggregated counts data.  This file size can be several hundred MBs, depending on the number of cases parsed.  Each row of the DataFrame represents an individual case (e.g. patient sample), and columns contain FPKM-UQ values associated with each Ensembl gene ID.  

Three additional columns will also be populated from the JSON metadata file: 'PtID', 'TumorStage', and 'SampleType'.  SampleType will be a two-digit code corresponding to the type of tissue from which the sample was acquired.  A list of SampleType codes is available at <https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes>  Some of the more common codes I've encountered are listed below.

	01	Primary Solid Tumor
	02	Recurrent Solid Tumor
	06	Metastatic
	11	Solid Tissue Normal
	

 

### 3. Extract values for genes of interest with ExtractTargets.py
The .pickle file that ParseAgg.py outputs contains FPKM-UQ values for all 60,486 Ensembl gene IDs from each sample.  ExtractTargets.py will create and output a new DataFrame containing aggregated data for only a specified subset of Ensembl gene IDs.  
ExtractTargets.py requires the following command-line arguments:

	-i <ValuesDataFrame.pickle>
	-t <TargetsList.txt>
	-o <Output directory>
	-p <Project name>

Thus, to parse the data in the above example:

	$ python ExtractTargets.py -i /TestData/Output/BRCA_aggregated_FPKM-UQ.pickle -t /TestData/targets.txt -o /TestData/Output -p BRCA

Where targets.txt is a list containing the Ensembl IDs and associated gene names, for example:

	F2RL3	  ENSG00000127533
	GAPDH	  ENSG00000111640
	WWTR1 (TAZ) ENSG00000018408

ExtractTargets.py will generate DataFrames containing raw count values for selected targets and a separate one containing log2(x+1) transformed values.  Both DataFrames will be saved as .pickle and as Excel spreadsheets (.xlsx) for use in subsequent analysis/visualization.

### 4. CorrelationAnalysis.py
Calculate Pearson's and Spearman correlation coefficients between a target gene of interest and all other genes in the input DataFrame.  CorrelationAnalysis.py only strictly requires an input dataframe .pickle:

		usage: CorrelationAnalysis.py [-h] -i INPUT_DF [-t TARGET_GENE]
                              [-o OUTPUT_DIR] [-p PROJECT_NAME] [-z Z_CUTOFF]
                              [--no_heat] [--no_corr] [--invert]
                              [-f FILE_FORMAT] [--corr_include_ctrl]
                              [--no_sorting]

	CorrelationAnalysis

	optional arguments:
	  -h, --help            show this help message and exit
	  -i INPUT_DF, --input_df INPUT_DF
	                        Path to df_log2.pickle file that was generated from
	                        ExtractTargets.py Ex: '-i BRCA_hippo_df_log2.pickle'
	  -t TARGET_GENE, --target_gene TARGET_GENE
	                        Target gene of interest! Must match name from the
	                        targets file. Default will choose firstgene in list
	                        Ex: '-t F2RL3'
	  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
	                        Folder to save output files to? Ex: '-o ./output'
	                        Default is current directory
	  -p PROJECT_NAME, --project_name PROJECT_NAME
	                        Short name to add on to the output files Ex: '-p
	                        BRCA_F2RL3_20160606' Default is target/date_time
	  -z Z_CUTOFF, --z_cutoff Z_CUTOFF
	                        Exclude samples with z-scores above andbelow a certain
	                        threshold. Default is 5. Set to 0 for no cutoff. Set
	                        to -1 to disable z-score conversion
	  --no_heat             Option to disable generation of heatmap
	  --no_corr             Option to disable generation of correlation data
	  --invert              Option to invert sorting of target gene on heatmap ie
	                        lowest expressers will be on top
	  -f FILE_FORMAT, --file_format FILE_FORMAT
	                        Option to save heatmap as different format (eps, ps,
	                        or pdf) default is pdf
	  --corr_include_ctrl   Option to include control samples in
	                        correlationcalculations. They are discluded by default
	  --no_sorting          Option to disable arranging heatmap with highest
	                        correlated genes from left to right. Sorting is on by
	                        default
	

CorrelationAnalysis.py outputs the Pearson's and Spearman correlation coefficients in a tab-delimited text document.  A heatmap is also generated (from z-scores of the input df, typically log2(FPKM-UQ) values) for a visual representation of trends in expression between your target gene and each of the individual genes in your set.
Future version will allow for more control over how the input data is processed, such as:  
<li>option to drop genes that have zero values over a certain threshold % of the samples  
<li>option to exclude certain sample types (ie drop metastatic samples).
 

## Authors

* **Thomas Smith** - [ThomasHSmith](https://github.com/ThomasHSmith)


## License
This work is licensed under the MIT License - see LICENSE.md for details.
