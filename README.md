# GCN
## 1. used of batch_gcn_segmentation.py

	python batch_gcn_segmentation.py -pm paths_matrices_sub_20.csv -pa paths_annot_sub_20.csv -o echantionage -rt 10

modification: add parameter repeat_time(defaut=1)

so the path



### package : 
#### pysradb :
version : 0.9.7
https://github.com/saketkc/pysradb


## 2.download the fastq from list of SRA

### Package:
#### SRA tools : 
version : > 2.10(because of the fasterq-dump)
https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump

impotant :change the default path to which SRA files are downloaded :  https://www.biostars.org/p/159950/

	path : echo \'/repository/user/main/public/root= "{new_path}"\' > $HOME/.ncbi/user-settings.mkfg' 

(new_path wil use for the next step: parameter path)
				
#### download_SRA.py : Script for download the fastq 

Automating downloads using Python ,Since there are lots of SRA files associated with our samples,it would take a long time to manually run prefetch and fastq-dump for all the files. 

To automate this process,the fonction will download each SRA file using prefetch and then run fasterq-dump. 
    
    Parameter:
    
		--list_SRA, type=check_file_path, help="list of SRA:SRRXXXXX(.txt)"
    	--fastq, type=str,default='fastq', help="folder name"
		--new_path, type=str,default='/home/storage_1/yuping/raw_data/',help="path to save the SRA fill (echo '/repository/user/main/public/root= new_path ' > $HOME/.ncbi/user-settings.mkfg) before run the script"
    	--e, type=int,default='20', help="threads"
		

