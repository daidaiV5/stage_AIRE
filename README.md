# research_dataset

selection_data.ipynb : research in mouse_transcript_v11.h5 and research in GEO dataset

# GCN
## 1. used of batch_gcn_segmentation.py

	python batch_gcn_segmentation.py -pm paths_matrices_sub_20.csv -pa paths_annot_sub_20.csv -o echantionage -rt 10

modification: 

1.add parameter repeat_time(defaut=1)

2.paths_matrices_sub_20.csv and paths_annot_sub_20.csv : separate by comma

3.-o: the path will used for the next step: same path for input_file

4.GCNScript3_JT.r:allowWGCNAThreads(nThreads=20)

## 2. used of concensus.py

	python concensus.py --input_file /home/storage_1/yuping/GCN/Liver/ --name_of_file filtered_all_cor_0.8.csv --threshold 0.7 --outdir /home/storage_1/yuping/compute_GCN/result_concensus3/

input:

--input_file : path to folder which contain all the result batch_gcn_segmentation.py( path should be the same  in batch_gcn_segmentation.py parameter -o)
  
--name_of_file : which file you want to do te combination(for example if we are intresting in 'liver_class1_male_filtered_all_cor_0.85.csv' we will use filtered_all_cor_0.85.csv)
    
--threshold : threshold of the ratio (dafault = 0.9): so we filtering the result which ratio > 0.9
    
--outdir : the path for output

output :

multilayer_edges.csv and multilayer_nodes.csv for each class of age

## 3.script for class of edges and nodes(abdel)






