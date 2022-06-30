# research_dataset

selection_data.ipynb : research in mouse_transcript_v11.h5 and research in GEO dataset

# GCN
## 1. used of pipeline_gcn.py
#### Script will calculate GCNs for provided matrices paths and segment the GCNs at provided thresholds.concensus of the random result


	usage: pipeline_gcn.py [-h] [-pm PATHS_MATRICES] [-o OUTPUT_DIRECTORY] [-pthr PEARSON_THRESHOLD]
                       [--repeat_times REPEAT_TIMES] [--number_of_each_GCN NUMBER_OF_EACH_GCN]
                       [--number_sample_min NUMBER_SAMPLE_MIN]
                       [--number_of_repeat_diffrente NUMBER_OF_REPEAT_DIFFRENTE] [--normal NORMAL]
		  
		  
		  
	eg. python pipeline_gcn.py -pm /home/storage_1/yuping/GCN/pipeline_GCN/dataset_test/paths_matrices_blood.csv -o /home/storage_1/yuping/GCN/pipeline_GCN/test/ -rt 20 -nG 10 -n 16 -norm False
	
(-rt:20 repetitions -nG ramdomly select 10 samples -n : smallest groupe contain 16 samples in the smallest group)




#### optional arguments:

  -pm PATHS_MATRICES, --paths_matrices PATHS_MATRICES
                        path to list of paths to count matrices (format organ comma sex comma age_class comma path)
			
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        path to the output directory
			
  -pthr PEARSON_THRESHOLD, --pearson_threshold PEARSON_THRESHOLD
                        Pearson correlation threshold (default: 0.8)
			
  --repeat_times REPEAT_TIMES, -rt REPEAT_TIMES
                        number of replications (default: 10)
			
  --number_of_each_GCN NUMBER_OF_EACH_GCN, -nG NUMBER_OF_EACH_GCN
                        number of sample for each GCN (default : 5)
			
  --number_sample_min NUMBER_SAMPLE_MIN, -n NUMBER_SAMPLE_MIN
                        number of samples of smallest groupe (default : 16)
			
  --number_of_repeat_diffrente NUMBER_OF_REPEAT_DIFFRENTE, -nd NUMBER_OF_REPEAT_DIFFRENTE
                        number maximun common between each repeat(default : 2)
			
  --normal NORMAL, -norm NORMAL
                        normalize with deseq or not (default = 'True')

#### Description of datatest:

dataset_test/paths_matrices_blood.csv : example for the path file

dataset_test/matrice/: two matrice used for test (only 40 gene in the matrice)

#### Desciption of script:

##### R:

GCNScript_normalization.r
  -- R version 3.6.3
  -- package : DESeq2 ‘1.26.0’


GCNScript4_JT.r
  -- R version 3.6.3
  -- package : WGCNA ‘1.70.3’

##### python : 

concensus_package.py : script for table concensus


  
  
  
