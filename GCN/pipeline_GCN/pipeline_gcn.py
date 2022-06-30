"""
pipeline used for the wgcna and concensus

"""

__author__ = "Yuping DAI"
__copyright__ = "Team AIRE, UMR 7205 ISYEB / Sorbonne Université"
__credits__ = ["Jérôme Teulière", "Eric Bapteste","Yuping DAI"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = ""
__status__ = "Development"

import argparse
import subprocess as sp
import os
import time
import datetime
from itertools import combinations
import pandas as pd
import numpy as np
import random
from random import randrange
import networkx as nx
import shutil
from concensus_package import MultiLayer


def parser():
    """
    Converts command line arguments to variables
    :return: variables
    """
    parser = argparse.ArgumentParser(description="Will calculate GCNs for provided matrices paths and "
                                                 "segment the GCNs at provided thresholds.""concensus of the random result")
    parser.add_argument("-pm", "--paths_matrices", type=str, help="path to list of paths to count matrices "
                                                                  "(format organ comma sex comma age_class comma path)")
    parser.add_argument("-o", "--output_directory", type=str, help="path to the output directory")
    parser.add_argument("-pthr", "--pearson_threshold", type=float, help=" Pearson correlation threshold "
                                                                         "(default: 0.8)", default=0.8)
    parser.add_argument("--repeat_times", "-rt", type=int, help="number of replications (default: 10)", default=10)
    parser.add_argument("--number_of_each_GCN","-nG",type=int,help="number of sample for each GCN (default : 5)", default=5)
    parser.add_argument("--number_sample_min","-n",type=int,help="number of samples of smallest groupe (default : 16)", default=16)
    parser.add_argument("--number_of_repeat_diffrente","-nd",type=int,help="number maximun common between each repeat(default : 2)",default=2)
    parser.add_argument("--normal","-norm",type=str,help="normalize with deseq or not (default = 'True')", default = 'True')
     
    
    args = parser.parse_args()
    paths_matrices = args.paths_matrices
    outdir = args.output_directory
    pthr = args.pearson_threshold
    number_sample=args.number_sample_min
    repeat_times=args.repeat_times
    number_of_each_GCN=args.number_of_each_GCN
    number_of_repeat_diffrente=args.number_of_each_GCN-args.number_of_repeat_diffrente
    normal=args.normal
    return paths_matrices,outdir, pthr, repeat_times,number_of_each_GCN,number_of_repeat_diffrente,number_sample,normal

############################################################################################
####################function of the random##################################################
############################################################################################


def combine(temp_list, n):
    ''' 
    combine : get all possible combinations of a list’s elements
    input:
    temp_list: list for which we will do the random
    n: number of random
    output :
    all possible combinations of list
    
    '''
    temp_list2 = []
    for c in combinations(temp_list, n):#Return n length subsequences of elements from the input iterable.
        temp_list2.append(c)
    return temp_list2

def choose_and_remove( items ):
    ''' 
	random one out of the list,and when chossing one times we will delete in the items

    '''
    index = random.randrange( len(items) )
    
    return items.pop(index)


def maximun_commen(random_result,list_final,nomber_common_maximun):
    '''
    compare the new random result with the old random 
	random result : all samples already present in the random result
	list final : new random result that we want to add in the random result
	nomber_common_maximun : number of the samples between two random 

    '''
    list_commen=[]
    set2 = set(random_result)
    for j in list_final:
        set1 = set(j)
        list_commen.append(len(set1&set2))
    print(list_commen)
    if [i for i in list_commen if i > nomber_common_maximun]: 
        return False
    else:
        return True

def list_to_file(liste_final,path):
    data=pd.DataFrame(liste_final)
    data.to_csv(path,sep=',',index=0,header=0)

############################################################################################
####################read and write files##################################################
############################################################################################

   

def file_to_list(path):
    """
    returns a list of lines (end of lines removed) from a file
    :param path: path to the file to read
    :return: list of strings
    """
    with open(path, "r") as f:
        lines = f.readlines()
        list_lines = [i.rstrip() for i in lines]
    return list_lines


def read_dictionary(path, sep=",", colindexk=(0, 1, 2), colindexv=3, skipfirst=False):
    """
    Read a csv file and creates a dictionary from two columns.
    :param path: path to the csv file
    :param sep: separator, comma by default
    :param colindexk: list of indices of the columns for dictionary keys
    :param colindexv: index of the column for dictionary values
    :param skipfirst: boolean, True to skip first line (column names), False by default
    :return: dictionary
    """
    with open(path, "r") as f:
        lines = f.readlines()
        list_lines = [i.rstrip() for i in lines]
    attr_dict = {}
    if skipfirst:
        for line in list_lines[1:]:
            k = (line.split(sep)[colindexk[0]], line.split(sep)[colindexk[1]], line.split(sep)[colindexk[2]])
            v = line.split(sep)[colindexv]
            attr_dict[k] = v
    else:
        for line in list_lines:
            k = (line.split(sep)[colindexk[0]], line.split(sep)[colindexk[1]], line.split(sep)[colindexk[2]])
            v = line.split(sep)[colindexv]
            attr_dict[k] = v
    return attr_dict


def make_output_file(outfile_name):
    """
    Creates an output file.
    :param outfile_name: name of, or path to, the output file
    :return: None
    """
    open(outfile_name, "w").close()


def make_output_directory(dir_path):
    """
    Creates all nested directories for a given directory path
    :param dir_path: path to the lowest level directory
    (i.e. : "./foo/bar/toto")
    :return: None
    """
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)


def execute(command_line):
    sp.run(command_line, shell=True)

############################################################################################
#############################for concensus##################################################
############################################################################################
def read_layer(file_layer):
    """
    input : file_layer, the name of a gene co-expression network edge file corresponding to an class
    return : layer, a coexpresson's networkx graph,
    """
    layer = nx.Graph()
    with open(file_layer) as infile:
        for line in infile:
            if line.startswith("from"):
                continue
            list_line = line.rstrip().split()
            # if float(list_line[2]) > 0.0: #only positive corellation of coexpression are take in account
            # node1, node2 = tuple(sorted([list_line[0].split(".")[0], list_line[1].split(".")[0]]))
            node1, node2 = tuple(sorted([list_line[0], list_line[1]]))
            cor = list_line[2]
            layer.add_edge(node1, node2, pearson=float(cor))
    return layer  



def read_input_file(input_file,name_of_file,repetitions):
    """
    a fonction for read all file for combination(consensus)
    
    input:
        input_file (str) : a path contain all fold which we want to combine(e.g: ~./20-29/male/wholeblood)
        name_of_file (str): a name of the file (e.g : wholeblood_male_20-29_sub_20_filtered_positive_cor_0.95.csv) 
       
    output : 
        list_layers : a list of networkx graphs in the order provided by the input file
        np.array(list_class): str array containing the number of the documents(number of repetition)
    
    There are thus as many lines as there are number of the documents.

    """
    list_path = []
    list_age_class = []
    files = os.listdir(input_file)
    for i in range(1,repetitions+1):
        i=str(i)
        input_file2=f'{input_file}/{i}/{name_of_file}'
        list_path.append(input_file2)
        list_age_class.append(i)
    list_layers = []
    for path in list_path:
        list_layers.append(read_layer(path))
    return list_layers, np.array(list_age_class)

############################################################################################
#############################for filtering##################################################
############################################################################################


def write_segmented(path_to_gcn, path_to_segmented, path_to_nodelist, absolute=True):
    with open(path_to_gcn, "r") as f, open(path_to_segmented, "a") as g, open(path_to_nodelist, "a") as h:
        all_nodes = set()
        for line in f:
            if "fromNode" in line:
                g.write("fromNode\ttoNode\tcor\n")
            else:
                line = line.split("\t")
                n1 = line[0]
                n2 = line[1]
                cor = float(line[2])
                if absolute:
                    if abs(cor) >= segthr:
                        g.write(f"{n1}\t{n2}\t{round(cor, 2)}\n")
                        all_nodes.add(n1)
                        all_nodes.add(n2)
                else:
                    if cor >= segthr:
                        g.write(f"{n1}\t{n2}\t{round(cor, 2)}\n")
                        all_nodes.add(n1)
                        all_nodes.add(n2)
        for node in all_nodes:
            h.write(f"{node}\n")


############################################################################################
################################ MAIN ######################################################
############################################################################################


# user input
paths_matrices,outdir, pthr, repeat_times,number_of_each_GCN,number_of_repeat_diffrente,number_sample,normal= parser()
outpaths = {}
radicals = {}
dico_new_matrices={}


for repeat_time in range(1,repeat_times+1):
    # hard-coded range of Pearson thresholds to segment the large GCN
    repeat_time=str(repeat_time)

    # create the output directory
    make_output_directory(outdir)
    

    # store paths to matrices and annotations in dictionaries with age-classe, sex, organ as keys
    dico_matrices = read_dictionary(paths_matrices, sep=",")

    for age_class, sex, organ in dico_matrices:
        dico_matrices[age_class,sex,organ]
        dico_new_matrices[age_class,sex,organ,repeat_time]=dico_matrices[age_class,sex,organ]

    # create subfolders for each age-class, to store in an outpaths dictionary:
    for age_class, sex, organ in dico_matrices:
        path = os.path.join(outdir, organ, sex, age_class)
        make_output_directory(path)
        outpaths[(age_class, sex, organ)] = path
        radical =f'{age_class}_{sex}_{organ}'
        radicals[(age_class, sex, organ)] = radical


# use R script to calculate GCNs

segmented = {}
nodelists = {}

### random sampling
list1 = list(range(1,number_sample+1,1))
end_list = [] #end_list is all possible combinations 
end_list.extend(combine(list1, number_of_each_GCN))

list_final=[] 
list_final.append(choose_and_remove(end_list))
i=0
while(i<repeat_times-1):
    random_result=choose_and_remove(end_list)
    x=maximun_commen(random_result,list_final,number_of_repeat_diffrente)
    if(x==True):
        list_final.append(random_result)
        i=i+1 

for (age_class, sex, organ), path_m in dico_matrices.items():
    stamp = int(time.time())
    print(datetime.datetime.fromtimestamp(stamp))
    outpath = outpaths[(age_class, sex, organ)]
    print(outpath)
    path_list=f'{outpath}/list_of_random.txt'
    list_to_file(list_final,path_list)
    radical = radicals[(age_class, sex, organ)]
    print(radical)

    print(normal)
    ####### normalization by DEseq2 ########

    command1=f"Rscript GCNScript_normalization.r --in_mat={path_m} --out={radical} --dir={outpath} --number_sample={number_sample} --norm={normal}"

    
    execute(command1)
    path_m2=f'{outpath}/{radical}_matrice_nor.csv'
    
    ####### WGCNA ########

    command = f"Rscript GCNScript4_JT.r --in_mat={path_m2} --out={radical} --dir={outpath} --pthr={pthr} --repeats={repeat_times} --random_result={path_list}"
    
    execute(command)
    ####### filtering of threshold of the correlation of pearson ########

    thrs = ["0.9"] #### change for filteing with differents thresholds
    for segthr in thrs:
        for repeat_time in range(1,repeat_times+1):
            repeat_time=str(repeat_time)
            segthr = float(segthr)
            print(f"Segmenting GCN at Pearson >= {segthr}")
            segpath1 = os.path.join(outpath, repeat_time,f"{radical}_filtered_all_cor_{segthr}.csv")
            make_output_file(segpath1)
            nodes1 = os.path.join(outpath, repeat_time,f"{radical}_filtered_all_cor_{segthr}_nodes.txt")
            make_output_file(nodes1)
            gcn = os.path.join(outpath,repeat_time,f"{radical}_gcn_edges_cor_{pthr}.csv")
            write_segmented(gcn, segpath1, nodes1, absolute=True)
            print(f"GCN {segthr}: done.")
   
    ####### concensus of the correlation of pearson ########
    liste_of_final_name=f'{outpath}/'
    for segthr in thrs:
        name_of_file=f'{radical}_filtered_all_cor_{segthr}.csv'
        list_layers, array_class = read_input_file(liste_of_final_name,name_of_file,repeat_times)
        liste_outdir=f'{outpath}/concensus_pearson_{segthr}'
        MultiLayer(list_layers, array_class,liste_outdir) 
     
    ####### removing the random result  ########
    for repeat_time in range(1,repeat_times+1):
        dest=f'{outpath}/{repeat_time}/'
        shutil.rmtree(dest, ignore_errors=True)
	


stamp = int(time.time())
print(datetime.datetime.fromtimestamp(stamp))
print("All done!")




