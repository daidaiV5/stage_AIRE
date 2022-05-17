import argparse
import os
import sys
import networkx as nx
import numpy as np
import pandas as pd
import statistics
import time
import datetime

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



def ratio_of_true(list_boolean,len_of_list):
    """
    input:
        list_boolean: list of booleans, descibes the presence or not of a
        given edge
        len_of_list : lengths of the list (total )
    output:
        value: ratio of edge are present to number of repetition.
    """
    x=list_boolean.count(True)
    value=x/len_of_list
    return value


def count_pearson(list_pearson):
    """
    the fonction uses to calculate the mean of the pearson for edgs presents.
    
    input:
        list_boolean: list of pearson for all repetition, if pearson = 0 that means the edge are absents in this repetition
        
    output:
        means of the pearson 
    """
    list_new_pearson=[]
    for i in list_pearson:
        if i>0:
            list_new_pearson.append(i)
    if len(list_new_pearson)==0:
        return 0
    else:
        return statistics.mean(list_new_pearson)

def create_directory(outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    return outdir

class MultiLayer:
    """
    a multilayer object is a dynamic gene coexpression graph evolving through repetition groupe
     attributes:
         - list_class: a str list, the list of repetition classes
         - multilayer: networkx graph whose edges represent a gene pair coexpressed within at one repetition.
         nodes are genes. nodes and edges have some attributes:
            - layers : an array list of boolean characterizing the presence or the absence as a function of the repetition/layer considered
            - ratio : list of tuples, a tuple corresponds to an interval of True value in the attribute layers
         - list_layers: a networkx graph list whose each element is the coexpression graph of the corresponding class(repetition),
         nodes are genes
         - nb_layer: int, the number of class(repetition)
    """

    """constructor"""

    def __init__(self, list_layers, array_class, threshold,outdir):
        """
        MultiLayer class constuctor
        input:
        list_layers : a list whose each element is the coexpression graph of the corresponding age class,
        nodes are genes
        array_age_class : an str list, the list of age classes
        """
        self.list_layers, self.classes = list_layers, array_class 
        self.nb_layer = len(self.classes)
        self.multilayer = nx.Graph()
        self.threshold = threshold
        self.compute_vectors()
        self.outdir = outdir
        self.compute_ranges_class()
        self.write_multi_edges_csv()

    def compute_vectors(self):
        """
        initializes the multilayer's attributes and computes the boolean vector
        for each edge and node from the attribute list_layers
        """
        for index, layer in enumerate(self.list_layers):
            for node1, node2, data in layer.edges(data=True):
                if self.multilayer.has_edge(node1, node2):
                    self.multilayer[node1][node2]["layers"][index] = True
                    self.multilayer[node1][node2]["pearson"][index] = float(data["pearson"])
                else:
                    self.multilayer.add_edge(node1, node2, layers=[False] * self.nb_layer, pearson=[0] * self.nb_layer)
                    self.multilayer[node1][node2]["layers"][index] = True
                    self.multilayer[node1][node2]["pearson"][index] = float(data["pearson"])
    

    def compute_ranges_class(self):
        """
        computes for each node and edge the atributes range (from layers) and class (from range)
        """
        edge_list=[]
        for node1, node2 in self.multilayer.edges:
            self.multilayer[node1][node2]["ratio"] = ratio_of_true(self.multilayer[node1][node2]["layers"],self.nb_layer)
            if(self.multilayer[node1][node2]["ratio"]<self.threshold):
                edge_list.append((node1, node2))
        self.multilayer.remove_edges_from(edge_list)
        
    def write_multi_edges_csv(self):
        """
        save the edge's multilayer and their attribute into csv format
        input:
            outdir: str, the out directory's name
        """
        subdir =  create_directory(self.outdir)
        with open(os.path.join(subdir, "multilayer_edges.csv"), "w") as outfilecsv:
            line = f'Node1,Node2,mean_pearson,{",".join(list(self.classes))},ratio\n'
            outfilecsv.write(line)
            for node1, node2 in self.multilayer.edges:
                layers = str(list(self.multilayer[node1][node2]['layers']))[1:-1]
                mean_pearson=str(count_pearson(list(self.multilayer[node1][node2]["pearson"])))
                line = f'{node1},{node2},{mean_pearson},{layers},{self.multilayer[node1][node2]["ratio"]}\n'
                outfilecsv.write(line.replace(" ", ""))
#         with open(os.path.join(subdir, "multilayer_nodes.csv"), "w") as outfilecsv:
#             line=f'NodeName\n'
#             outfilecsv.write(line)
#             for node in self.multilayer.nodes:
#                 line=f'{node}\n'
#                 outfilecsv.write(line.replace(" ", ""))
            


def args_check():
    parser = argparse.ArgumentParser(description= """
    This script will use for combination
    
    eg:
    
    python concensus.py --input_file /home/storage_1/yuping/GCN/Liver/ --name_of_file filtered_all_cor_0.8.csv --threshold 0.7 --outdir /home/storage_1/yuping/compute_GCN/result_concensus3/
   
    --input_file : path to folder which contain all the result batch_gcn_segmentation.py( path should be the same  in batch_gcn_segmentation.py parameter -o)
  
    --name_of_file : which file you want to do te combination(for example if we are intresting in 'liver_class1_male_filtered_all_cor_0.85.csv' we will use filtered_all_cor_0.85.csv)
    
    --threshold : threshold of the ratio (dafault = 0.9): so we filtering the result which ratio > 0.9
    
    --outdir : the path for output
    
    

"""
    ,formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument("--input_file", type=str, help="path folder")
    parser.add_argument("--name_of_file", type=str, help="file name 'filtered_all_cor_0.8.csv'")
    parser.add_argument("--threshold", type=float,default=0.9)
    parser.add_argument("--repetition", type=float,default=20)
    parser.add_argument("--outdir", type=str)
    parameters = parser.parse_args()
    return parameters                
            
            
def main():
    parameters = args_check()
    print('------begin----------')
    input_file=parameters.input_file
    for i in os.listdir(input_file):
        liste_of_class_name=f'{input_file}{i}/'
        print(liste_of_class_name)
        for j in os.listdir(liste_of_class_name):
            liste_of_sex_name=f'{input_file}{i}/{j}/'
            for h in os.listdir(liste_of_sex_name):
                stamp = int(time.time())
                print(datetime.datetime.fromtimestamp(stamp))
                liste_of_final_name=f'{input_file}{i}/{j}/{h}'
                name=parameters.name_of_file
                h=h.lower()
                j=j.capitalize()
                name_of_file=f'{h}_{i}_{j}_{name}'
                print(f'{h}_{i}_{j}_{name}')
                print('begin""')
                list_layers, array_class = read_input_file(liste_of_final_name,name_of_file,parameters.repetition)
                liste_outdir=f'{parameters.outdir}/{i}_{j}_{parameters.threshold}'
                MultiLayer(list_layers, array_class, parameters.threshold,liste_outdir)
                print(liste_of_final_name)
    print('------success----------')
    stamp = int(time.time())
    print(datetime.datetime.fromtimestamp(stamp))


if __name__ == "__main__":
    main() 
    
    
    
    
    
