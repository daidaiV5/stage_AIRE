import argparse
import subprocess as sp
import os
import time
import datetime
from itertools import combinations
import pandas as pd
import random
from random import randrange

def parser():
    """
    Converts command line arguments to variables
    :return: variables
    """
    parser = argparse.ArgumentParser(description="Will calculate GCNs for provided matrices paths and "
                                                 "segment the GCNs at provided thresholds.")
    parser.add_argument("-pm", "--paths_matrices", type=str, help="path to list of paths to count matrices "
                                                                  "(format organ tab sex tab age_class tab path)")
    parser.add_argument("-pa", "--paths_annot", type=str, help="path to list of paths to annotation files "
                                                               "(format organ tab sex tab age_class tab path)")
    parser.add_argument("-o", "--output_directory", type=str, help="path to the output directory")
    parser.add_argument("-thr", "--student_threshold", type=float, help="Student P value threshold for correlations"
                                                                        " (default: 0.05)", default=0.05)
    parser.add_argument("-pthr", "--pearson_threshold", type=float, help=" Pearson correlation threshold "
                                                                         "(default: 0.5)", default=0.5)
    parser.add_argument("-min", "--min_reads", type=int, help="minimum number of reads to consider (default: 20)", default=20)
    parser.add_argument("--repeat_times", "-rt", type=int, help="number of replications (default: 1)", default=10)
    parser.add_argument("--number_of_each_GCN","-nG",type=int,help="number of sample for each GCN(default : 5)", default=5)
    parser.add_argument("--number_of_repeat_diffrente","-nd",type=int,help="number maximun common between each repeat(default : 2)",default=2)    
    
    args = parser.parse_args()
    paths_matrices = args.paths_matrices
    paths_annot = args.paths_annot
    outdir = args.output_directory
    thr = 1 - args.student_threshold
    pthr = args.pearson_threshold
    min_reads = args.min_reads
    repeat_times=args.repeat_times
    number_of_each_GCN=args.number_of_each_GCN
    number_of_repeat_diffrente=args.number_of_repeat_diffrente
    return paths_matrices, paths_annot, outdir, thr, pthr, min_reads,repeat_times,number_of_each_GCN,number_of_repeat_diffrente


def combine(temp_list, n):
    ''' 
    combine : get all possible combinations of a list’s elements
    
    temp_list: list for which we will do the random
    n: number of random
    
    '''
    temp_list2 = []
    for c in combinations(temp_list, n):#Return n length subsequences of elements from the input iterable.
        temp_list2.append(c)
    return temp_list2

def choose_and_remove( items ):
    ''' 
    
    '''
    index = random.randrange( len(items) )
    
    return items.pop(index)


def maximun_commen(random_result,list_final,nomber_common_maximun):
    list_commen=[]
    set2 = set(random_result)
    for j in list_final:
        set1 = set(j)
        list_commen.append(len(set1&set2))
    print(list_commen)
    if [i for i in list_commen if i <= nomber_common_maximun]: 
        return True
    else:
        return False

def list_to_file(liste_final,path):
    data=pd.DataFrame(liste_final)
    data.to_csv(path,sep=',',index=0,header=0)
    
    
def compture_sample(annotation):
    '''
    to know how many sample in the matrice 
    
    '''
    df=pd.read_csv(annotation)
    return len(df)    
    

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
paths_matrices, paths_annot, outdir, thr, pthr, min_reads,repeat_times,number_of_each_GCN,number_of_repeat_diffrente = parser()
outpaths = {}
radicals = {}
dico_new_matrices={}

    


for repeat_time in range(1,repeat_times+1):
    # hard-coded range of Pearson thresholds to segment the large GCN
    repeat_time=str(repeat_time)
    thrs = ["0.8","0.85","0.9", "0.95"]

    # create the output directory
    make_output_directory(outdir)
    

    # store paths to matrices and annotations in dictionaries with age-classe, sex, organ as keys
    dico_matrices = read_dictionary(paths_matrices, sep=",")
    dico_annot = read_dictionary(paths_annot, sep=",")

    for age_class, sex, organ in dico_matrices:
        dico_matrices[age_class,sex,organ]
        dico_new_matrices[age_class,sex,organ,repeat_time]=dico_matrices[age_class,sex,organ]

    # create subfolders for each age-class, to store in an outpaths dictionary:
    for age_class, sex, organ in dico_matrices:
        path = os.path.join(outdir, organ, sex, age_class)
        make_output_directory(path)
        outpaths[(age_class, sex, organ)] = path


    # define radicals for file names in a third dictionary
    for (age_class, sex, organ), path in dico_matrices.items():
        radical = path.split("/")[-1].split(".")[0]
        radicals[(age_class, sex, organ)] = radical


# use R script to calculate GCNs

segmented = {}
nodelists = {}

for (age_class, sex, organ), path_m in dico_matrices.items():
    stamp = int(time.time())
    print(datetime.datetime.fromtimestamp(stamp))
    path_a = dico_annot[(age_class, sex, organ)]
    number_sample=compture_sample(path_a)
    list1 = list(range(1,number_sample+1,1))
    end_list = [] #end_list is all possible combinations 
    end_list.extend(combine(list1, number_of_each_GCN))

    list_final=[] # list_final for 
    list_final.append(choose_and_remove(end_list))
    i=0
    while(i<repeat_times-1):
        random_result=choose_and_remove(end_list)
        x=maximun_commen(random_result,list_final,number_of_repeat_diffrente)
        if(x==True):
            list_final.append(random_result)
            i=i+1
        
    outpath = outpaths[(age_class, sex, organ)]
    print(outpath)
    path_list=f'{outpath}/list_of_random.txt'
    list_to_file(list_final,path_list)
    

    radical = radicals[(age_class, sex, organ)]
    
    command = f"Rscript GCNScript4_JT.r --in_mat={path_m} --out={radical} --dir={outpath} --min={min_reads} --thr={thr} --pthr={pthr} --repeats={repeat_time} --random_result={path_list}"
    
    execute(command)

# segment the resulting GCN, extract node lists, and  keep track of paths in a dictionary

    
    for segthr in thrs:
        segthr = float(segthr)
        print(f"Segmenting GCN at Pearson >= {segthr}")
        segpath1 = os.path.join(outpath, f"{radical}_filtered_all_cor_{segthr}.csv")
        segpath2 = os.path.join(outpath, f"{radical}_filtered_positive_cor_{segthr}.csv")
        make_output_file(segpath1)
        make_output_file(segpath2)
        segmented[(age_class, sex, organ, repeat_time ,segthr)] = {"all": segpath1, "pos": segpath2}
        nodes1 = os.path.join(outpath, f"{radical}_filtered_all_cor_{segthr}_nodes.txt")
        nodes2 = os.path.join(outpath, f"{radical}_filtered_positive_cor_{segthr}_nodes.txt")
        make_output_file(nodes1)
        make_output_file(nodes2)
        nodelists[(age_class, sex, organ, repeat_time ,segthr)] = {"all": nodes1, "pos": nodes2}

        gcn = os.path.join(outpath, f"{radical}_gcn_edges_cor_{pthr}.csv")
        write_segmented(gcn, segpath1, nodes1, absolute=True)
        write_segmented(gcn, segpath2, nodes2, absolute=False)
        print(f"GCN {segthr}: done.")

print("All done!")


# export paths to edgefiles in a new file
path_file = f"./path_to_edgefiles.csv"
make_output_file(path_file)
with open(path_file, "a") as f:
    f.write(f"age_class\tsex\torgan\tPearson_thr\tcorrelations\tpath\n")
    for (age_class, sex, organ,repeat_time), dir in outpaths.items():
        radical = radicals[(age_class, sex, organ,repeat_time)]
        filepath = os.path.join(dir, f"{radical}_gcn_edges_cor_{pthr}.csv")
        f.write(f"{age_class}\t{sex}\t{organ}\t{pthr}\tall\t{filepath}\n")
    for (age_class, sex, organ,repeat_time,segthr), paths in segmented.items():
        for cor, path in paths.items():
            f.write(f"{age_class}\t{sex}\t{organ}\t{segthr}\t{cor}\t{path}\n")

# path_file = f"./path_to_nodelists.csv"
# make_output_file(path_file)
# with open(path_file, "a") as f:
#     f.write(f"age_class\tsex\torgan\tPearson_thr\tcorrelations\tpath\n")
#     for (age_class, sex, organ,repeat_time), dir in outpaths.items():
#         radical = radicals[(age_class, sex, organ,repeat_time)]
#         filepath = os.path.join(dir, f"{radical}_gcn_nodes_cor_{pthr}.txt")
#         f.write(f"{age_class}\t{sex}\t{organ}\t{pthr}\tall\t{filepath}\n")
#     for (age_class, sex, organ,repeat_time,segthr), paths in nodelists.items():
#         for cor, path in paths.items():
#             f.write(f"{age_class}\t{sex}\t{organ}\t{segthr}\t{cor}\t{path}\n")


