# -*- coding: utf-8 -*-
"""
Created on Mon May 25 02:08:53 2020

@author: Heeyoun
"""

import os
import matplotlib.pyplot as plt
import seaborn as sns; sns.set(color_codes=True)
import pandas as pd
#import math

def search(path):
    res = []
    for root, dirs, files in os.walk(path):
        rootpath = os.path.join(os.path.abspath(path), root)
        
        for file in files:
            filepath = os.path.join(rootpath, file)
            
            res.append(filepath)
    return res


def violinplot(data_df, title):
    fig = plt.figure(figsize=(6,6), dpi=600)
    
    
    
    sns.violinplot(x = 'GO Class', y = 'Cscore', data = data_df, color = '0.8', order = ['Biological Process', 'Cellular Component', 'Molecular Function'])
    sns.swarmplot(x = 'GO Class', y = 'Cscore', data = data_df, size = 3, order = ['Biological Process', 'Cellular Component', 'Molecular Function'])    

    #plt.xlabel(GO_class)
    plt.ylabel('Cscore')
    plt.title(title)
    #plt.legend(legends, loc="lower right")
    plt.savefig(title + "_violinplot.png")
    
    plt.clf()
    plt.cla()
    plt.close(fig)

def pathfinder(Stream_list, root, GO_path_dic):
    #path_ = []
    stems = []
    
    for stream in Stream_list:
        streams = stream.split('-')

        if streams[1] == root:
            #path_.append(stream)
            stems.append(streams[0])
    GO_path_dic [root] = stems
    
    for stem in stems:
        pathfinder (Stream_list, stem, GO_path_dic)
      


def pathfinder2 (GO_end, GO_end_path):

    for key in GO_path_dic.keys():
        daughter_list = GO_path_dic [key]
        if GO_end in daughter_list:
            GO_end_path.append(key)
            pathfinder2 (key, GO_end_path)
    
    

def svg_file_parsing(svg_file_lines, class_):
    
    GO_stream_line = ''
    GO_streams = []
    for line in svg_file_lines:
        if "<!--" in line and "-->"in line:
            GO_stream_line = line[:-1]
        elif 'class="edge"' in line:
            GO_stream_lines = GO_stream_line.split(';')
            GO0__ = GO_stream_lines[0].split('&')
            GO0_ = GO0__[0].split(' ')
            GO0 = GO0_[1]
            
            GO1_ = GO_stream_lines[2].split(' ')
            GO1 = GO1_[0]
            GO_streams.append(GO0 + '-' + GO1)
            
    root_GO = GO_root[class_]

    pathfinder(GO_streams, root_GO, GO_path_dic)

            
            
GO_node_cut = 3

GO_root = {'MF':'GO:0003674', 'BP':'GO:0008150', 'CC':'GO:0005575'}

nowdir = os.getcwd()
file_list = search(nowdir)
sample_proteins_dic = {}
all_proteins = []
all_GO_dic = {}
pathdata = []


for file_ in file_list:
    if "GOsearchresult_final" in file_ and ".csv" in file_ :
        csv_dic = {}
        file_path = file_.split('\\')
        gene_name = file_path[-3]
        file_names = file_path[-1].split('.')
        file_names_ = file_names[0].split('_')
        
        csv_file = open(file_)
        csv_file_lines = csv_file.readlines()
        for line in csv_file_lines:
            lines = line.split('\t')
        
        csv_file_df = pd.read_table(file_,header = None)
        csv_file_df.columns = ['GO_index', 'class', 'score', 'GO_names']
        csv_dic = csv_file_df.set_index('GO_index').T.to_dict('list')
        all_GO_dic [gene_name + '_' + file_names_[-1]] = csv_dic
        
        svg_file = open(file_[:-4] + '.svg')
        svg_file_lines = svg_file.readlines()
        GO_path_dic = {}
        svg_file_parsing (svg_file_lines, file_names_[2])
        

        path_list = []
        root = GO_root [file_names_[2]]
        
        GO_end_list = []
        for key in GO_path_dic.keys():
            daughter_list = GO_path_dic [key]
            if len(daughter_list) == 0:
                GO_end_list.append(key)
        endpath_list_all = []
        
        for GO_end in GO_end_list:
            endpath_list = []
            GO_end_path = [GO_end]
            pathfinder2 (GO_end, GO_end_path)
            #print GO_end_path
            index_list = list(filter(lambda x: GO_end_path[x] == root, range(len(GO_end_path))))
            before_index = 0
            
            for index_ in index_list:
                if index_ == len(GO_end_path):
                    pass
                else:
                    

                    
                    if before_index == 0:
                        #print index_, GO_end_path[before_index:index_+1]
                        endpath_list.append(GO_end_path[before_index:index_+1])
                        pass
                    else:
                        last_path = endpath_list[-1]
                        #print 'last_path', last_path
                        GO_found_list = GO_path_dic [GO_end_path[before_index]]
                        for GO_found in GO_found_list:
                            if GO_found in last_path:
                                new_index = last_path.index(GO_found)
                                if new_index == 0:
                                    new_list = [last_path [new_index]] + GO_end_path[before_index:index_+1]
                                else:
                                    new_list = last_path [:new_index + 1] + GO_end_path[before_index:index_+1]
                                #print GO_end_path[before_index], GO_found, new_index, ';'.join(new_list)
                                endpath_list.append(new_list)
                                
                            else:
                                pass

                    before_index = index_ + 1
            
            GO_path_element_list = []
            for endpath_ in endpath_list:
                GO_path_element_list = GO_path_element_list + endpath_[:- GO_node_cut]
                
 
            
            GO_path_element_list = list(set(GO_path_element_list))
            if len(GO_path_element_list) == 0:
                pass
            else:
                GO_end_series = csv_file_df[csv_file_df ['GO_index'] == GO_end]
                
                if file_names_[2] == 'BP':
                    GO_class = 'Biological Process'
                elif file_names_[2] == 'CC':
                    GO_class = 'Cellular Component'
                elif file_names_[2] == 'MF':
                    GO_class = 'Molecular Function'
                
                
                if GO_end_series ['score'].values == []:
                    pass
                else:
                    for val in GO_end_series ['score'].values:
                        pathdata.append([gene_name, GO_class, GO_end, val, ''.join(GO_end_series ['GO_names'].values), len(endpath_list), len(GO_path_element_list), ';'.join(GO_path_element_list)])

pathdata_df = pd.DataFrame(pathdata, columns = ['Genes', 'GO Class', 'End of GO Path', 'Cscore', 'GO Description', 'No. of GO Path', 'No. of GO elements', 'GO elements'])

pathdata_df.to_excel('GO_path_table.xlsx')





GO_info_dic = {}
data = []
protein_lists = []
for key in all_GO_dic.keys():
    keys = key.split('_')
    protein_lists.append(keys[0])
    if keys [-1] == 'BP':
        GO_class = 'Biological Process'
    elif keys [-1] == 'CC':
        GO_class = 'Cellular Component'
    elif keys [-1] == 'MF':
        GO_class = 'Molecular Function'

    GO_dic = all_GO_dic [key]

    
    for GO_ in GO_dic.keys():

        GO_list = GO_dic[GO_]

        data.append([keys[0]+'_'+keys[1], GO_class, GO_, GO_list[1], GO_list[2]])
        GO_info_dic [GO_] = GO_list[2]

data_df = pd.DataFrame (data, columns = ['Protein_Acc', 'GO Class', 'GO_ID', 'Cscore', 'GO_description'])

data_df.to_excel('all_GO_lists_gold_silver.xlsx')
protein_lists = list(set(protein_lists))





GO_class = ['Biological Process', 'Cellular Component', 'Molecular Function']

        
violinplot(data_df, 'All GO') 


        
violinplot(pathdata_df, 'End of GO path') 
           
            
        
        
        
            
        