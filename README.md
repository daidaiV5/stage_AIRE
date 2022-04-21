# research_dataset

selection_data.ipynb : research in mouse_transcript_v11.h5 and research in GEO dataset

# GCN
## 1. used of batch_gcn_segmentation.py

	python batch_gcn_segmentation.py -pm paths_matrices_sub_20.csv -pa paths_annot_sub_20.csv -o echantionage -rt 10

modification: 

1.add parameter repeat_time(defaut=1)
2.paths_matrices_sub_20.csv and paths_annot_sub_20.csv : separate by comma
3.-o: the path will used for the next step: same path for input_file

## 2. used of concensus.py


		

