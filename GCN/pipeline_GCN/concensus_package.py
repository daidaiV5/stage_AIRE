import argparse
import os
import sys
import networkx as nx
import numpy as np
import pandas as pd
import statistics
import time
import datetime


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
        if i != 0:
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

    def __init__(self, list_layers, array_class,outdir):
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
        self.threshold = ["0.7","0.8","0.9","1"]
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
        for node1, node2 in self.multilayer.edges:
            self.multilayer[node1][node2]["ratio"] = ratio_of_true(self.multilayer[node1][node2]["layers"],self.nb_layer)

        
    def write_multi_edges_csv(self):
        """
        save the edge's multilayer and their attribute into csv format
        input:
            outdir: str, the out directory's name
        """
        subdir =  create_directory(self.outdir)
        for i in self.threshold:
            name_threshold=f'multilayer_edges_{i}.csv'
            with open(os.path.join(subdir, name_threshold), "w") as outfilecsv:
                line=f'Node1,Node2,mean_pearson\n'
                outfilecsv.write(line)
                for node1, node2 in self.multilayer.edges:
                    i=float(i)
                    if(self.multilayer[node1][node2]["ratio"]>=i):
                        layers = str(list(self.multilayer[node1][node2]['layers']))[1:-1]
                        mean_pearson=str(count_pearson(list(self.multilayer[node1][node2]["pearson"])))
                        line = f'{node1},{node2},{mean_pearson}\n'
                        outfilecsv.write(line.replace(" ", ""))

        with open(os.path.join(subdir, "multilayer_edges_row.csv"), "w") as outfilecsv:
            line = f'Node1,Node2,mean_pearson,{",".join(list(self.classes))},ratio,{",".join(list(self.classes))}\n'
            outfilecsv.write(line)
            for node1, node2 in self.multilayer.edges:
                layers = str(list(self.multilayer[node1][node2]['layers']))[1:-1]
                mean_pearson=str(count_pearson(list(self.multilayer[node1][node2]["pearson"])))
                line = f'{node1},{node2},{mean_pearson},{layers},{self.multilayer[node1][node2]["ratio"]},{self.multilayer[node1][node2]["pearson"]}\n'
                outfilecsv.write(line.replace(" ", ""))
                
#         with open(os.path.join(subdir, "multilayer_nodes.csv"), "w") as outfilecsv:
#             line=f'NodeName\n'
#             outfilecsv.write(line)
#             for node in self.multilayer.nodes:
#                 line=f'{node}\n'
#                 outfilecsv.write(line.replace(" ", ""))
            




    
    
    
