# -*- coding: utf-8 -*-
"""
Created on Mon May 25 02:08:53 2020

@author: Heeyoun
"""

import os
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import plot_precision_recall_curve


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set(color_codes=True)
import pandas as pd
from statannot import add_stat_annotation
import pickle
#import math
import ssl
import json

import requests, sys



def search(path):
    res = []
    for root, dirs, files in os.walk(path):
        rootpath = os.path.join(os.path.abspath(path), root)
        
        for file in files:
            filepath = os.path.join(rootpath, file)
            
            res.append(filepath)
    return res

def ROC_curve(fpr, tpr, roc_auc, grades, title):
    fig = plt.figure(figsize=(6,6), dpi=600)
    #colors = ['r', 'g', 'b']
    i = 0
    legends = []
    for grade in grades:
        legends.append(str(round(roc_auc[i],2)) + ' ' + grade)
        lw = 0.5
        plt.plot(fpr[i], tpr[i], lw=lw)
        i = i + 1
    
    plt.plot([0, 1], [0, 1], color='navy', lw=lw *2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title)
    plt.legend(legends, loc="lower right")
    plt.savefig(title + " ROC Curve.png")
    
    plt.clf()
    plt.cla()
    plt.close(fig)

def violinplot(data_df, grade, title):
    fig = plt.figure(figsize=(6,6), dpi=600)
    
    order = ['1', '0']
    
    ax = sns.violinplot(x = grade, y = 'Cscore', data = data_df, color = '0.8')
    ax = sns.swarmplot(x = grade, y = 'Cscore', data = data_df, size = 3)    
    
    ax, test_results = add_stat_annotation(ax, data = data_df, x = grade, y = 'Cscore', order = order, box_pairs = [('1', '0')], test = 'Mann-Whitney', text_format = 'star', loc='outside', verbose = 2)
    
    # text_format='full', loc='inside', comparisons_correction=None,
    #                 line_offset_to_box=0.2, line_offset=0.1, line_height=0.05, text_offset=8
    
    plt.xlabel(grade)
    plt.ylabel('Cscore')
    plt.title(title)
    #plt.legend(legends, loc="lower right")
    plt.savefig(title + '_' + grade + "_violinplot.png")
    
    plt.clf()
    plt.cla()
    plt.close(fig)

def violinplot2(data_df, grade, title):
    fig = plt.figure(figsize=(6,6), dpi=600)
    ax = sns.violinplot(x = grade, y = 'No. of GO elements', data = data_df, color = '0.8')
    ax = sns.swarmplot(x = grade, y = 'No. of GO elements', data = data_df, size = 3)    
    order = ['1', '0']
    ax, test_results = add_stat_annotation(ax, data = data_df, x = grade, y = 'Cscore', order = order, box_pairs = [('1', '0')], test = 'Mann-Whitney', text_format = 'star', loc='outside', verbose = 2)
    
    plt.xlabel(grade)
    plt.ylabel('No. of GO elements')
    plt.title(title)
    #plt.legend(legends, loc="lower right")
    plt.savefig(title + '_' + grade + "_No. of GO elements_violinplot.png")
    
    plt.clf()
    plt.cla()
    plt.close(fig)

def FDR_plot(data_df, title):
    fig = plt.figure(figsize=(6,6), dpi=600)
    sns.lineplot(x = 'Cscore', y = 'fdr', data = data_df)
    #sns.swarmplot(x = grade, y = 'No. of GO elements', data = data_df, size = 3)    

    plt.xlabel('Cscore')
    plt.ylabel('FDR')
    plt.title(title)
    #plt.legend(legends, loc="lower right")
    plt.savefig(title + "_fdrplot.png")
    
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

def ROC_curve2(fpr, tpr, roc_auc, GO_class, grade, title):
    fig = plt.figure(figsize=(6,6), dpi=600)
    #colors = ['r', 'g', 'b']
    i = 0
    legends = []
    GOnames = ['BP', 'CC', 'MF']
    color=['blue', 'orange', 'purple']
    for GO_ in GO_class:
        legends.append(str(round(roc_auc[i],2)) + ' ' + GOnames[i])
        lw = 0.5
        plt.plot(fpr[i], tpr[i], color=color[i], lw=2)
        i = i + 1
    
    plt.plot([0, 1], [0, 1], color='navy', lw=lw *2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(grade)
    plt.legend(legends, loc="lower right")
    plt.savefig(title + '_' + grade + " ROC Curve.png")
    
    plt.clf()
    plt.cla()
    plt.close(fig)
    
def PR_curve (dddata_df, classes, title):    
    
    plt.figure(figsize=(7, 8))
    f_scores = np.linspace(0.2, 0.8, num=4)
    lines = []
    labels = []
    for f_score in f_scores:
        x = np.linspace(0.01, 1)
        y = f_score * x / (2 * x - f_score)
        l, = plt.plot(x[y >= 0], y[y >= 0], color='gray', alpha=0.2)
        plt.annotate('f1={0:0.1f}'.format(f_score), xy=(0.9, y[45] + 0.02))
    
    color=['blue', 'orange']
    i = 0

    for grade in classes:
        
        average_precision = average_precision_score(dddata_df[grade], dddata_df['Cscore'])
        print('Average precision-recall score: {0:0.2f}'.format(average_precision))
        labels.append(grade + ' Averaged PR score: ' + str(round(average_precision,3)))
        precision, recall, _ = precision_recall_curve(dddata_df[grade], dddata_df['Cscore'])
        #disp.ax_.set_title('2-class Precision-Recall curve: '
        #                   'AP={0:0.2f}'.format(average_precision))
        l, = plt.plot(recall, precision, color=color[i], lw=2)
        lines.append(l)
        i = i + 1
        
    fig = plt.gcf()
    fig.subplots_adjust(bottom=0.25)
    #l.legend(classes, loc="lower right")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    
    plt.title('Precision-Recall curve of ' + title )
    plt.legend(lines, labels, loc="lower right")        
    
    #plt.title(title)
    #plt.legend(legends, loc="lower right")
    plt.savefig(title + "_PR_curve.png", dpi = 600)
    
    plt.clf()
    plt.cla()
    plt.close(fig)

def PR_curve2 (ddata_df, GO_classes, grade, title):    
    GO_names = ['BP', 'CC', 'MF']    
    plt.figure(figsize=(7, 8))
    f_scores = np.linspace(0.2, 0.8, num=4)
    lines = []
    labels = []
    for f_score in f_scores:
        x = np.linspace(0.01, 1)
        y = f_score * x / (2 * x - f_score)
        l, = plt.plot(x[y >= 0], y[y >= 0], color='gray', alpha=0.2)
        plt.annotate('f1={0:0.1f}'.format(f_score), xy=(0.9, y[45] + 0.02))
    
    color=['blue', 'orange', 'purple']
    i = 0

    for GO_ in GO_classes:
        dddata_df = ddata_df[ddata_df["GO_Class"] == GO_]
        average_precision = average_precision_score(dddata_df[grade], dddata_df['Cscore'])
        print('Average precision-recall score: {0:0.2f}'.format(average_precision))
        labels.append(GO_names[i] + ' Averaged PR score: ' + str(round(average_precision,3)))
        precision, recall, _ = precision_recall_curve(dddata_df[grade], dddata_df['Cscore'])
        #disp.ax_.set_title('2-class Precision-Recall curve: '
        #                   'AP={0:0.2f}'.format(average_precision))
        l, = plt.plot(recall, precision, color=color[i], lw=2)
        lines.append(l)
        i = i + 1
        
    fig = plt.gcf()
    fig.subplots_adjust(bottom=0.25)
    #l.legend(classes, loc="lower right")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    
    plt.title('Precision-Recall curve of ' + title + '_' + grade)
    plt.legend(lines, labels, loc="lower right")        
    
    #plt.title(title)
    #plt.legend(legends, loc="lower right")
    plt.savefig(title + '_' + grade + "_PR_curve.png", dpi = 600)
    
    plt.clf()
    plt.cla()
    plt.close(fig)

           
GO_node_cut = 3

GO_root = {'MF':'GO:0003674', 'BP':'GO:0008150', 'CC':'GO:0005575'}

nowdir = os.getcwd()
file_list = search(nowdir)
sample_proteins_dic = {}
all_proteins = []
all_GO_dic = {}
pathdata = []

runlist_df = pd.read_excel('ITASSER_runlist.xlsx', index_col = 0)

runlist_dic = runlist_df.T.to_dict('list')

ancestor_dic = {}

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
        

        runlist_dic2 = runlist_dic [gene_name]
        try:
            gold_list = runlist_dic2 [1].split(';')

            ancestor_gold_list = []
        #     for gold_ in gold_list:
        #         goldies = gold_.split(':')
        #         requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO%3A" + goldies[1] + "/ancestors?relations=is_a%2Cpart_of%2Coccurs_in%2Cregulates"
        #         r = requests.get(requestURL, headers={ "Accept" : "application/json"}, verify = False)
                
        #         if not r.ok:
        #           r.raise_for_status()
        #           sys.exit()
                
        #         responseBody = json.loads(r.text)
        #         result_dic = responseBody.get('results')[0]
                
        #         try:
        #             ancestor_gold_list = ancestor_gold_list + result_dic ['ancestors']
        #             ancestor_dic [gold_] = result_dic ['ancestors']
        #         except:
        #             pass            
            
            
        #     #glod_list = gold_list + ancestor_dic [gene_name][0]
        except:
             glod_list = []
        # ancestor_gold_list = list(set(ancestor_gold_list))
        
        
        
        try:
            silver_list = runlist_dic2 [3].split(';')

            ancestor_silver_list = []
        #     for gold_ in silver_list:
        #         goldies = gold_.split(':')
        #         requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO%3A" + goldies[1] + "/ancestors?relations=is_a%2Cpart_of%2Coccurs_in%2Cregulates"
        #         r = requests.get(requestURL, headers={ "Accept" : "application/json"}, verify = False)
                
        #         if not r.ok:
        #           r.raise_for_status()
        #           sys.exit()
                
        #         responseBody = json.loads(r.text)
        #         result_dic = responseBody.get('results')[0]
                
        #         try:
        #             ancestor_silver_list = ancestor_silver_list + result_dic ['ancestors']
        #             ancestor_dic [gold_] = result_dic ['ancestors']
        #         except:
        #             pass            
    
            
            
            
        except:
             silver_list = []

        # ancestor_silver_list = list(set(ancestor_silver_list))

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
                            #print index_, GO_end_path [:new_index] + GO_end_path[before_index:index_+1]
#                        GO_end_path2 = [GO_end_path[before_index]]
#                        pathfinder2 (GO_end_path[before_index], GO_end_path2)
#                        print GO_end_path2
                    before_index = index_ + 1
            
            GO_path_element_list = []
            for endpath_ in endpath_list:
                GO_path_element_list = GO_path_element_list + endpath_[:- GO_node_cut]
                
            gold = False
            silver = False
            #panther = 0
#========================================================================            
            gold_list = gold_list + ancestor_gold_list
            silver_list = silver_list + ancestor_silver_list            
#=====================================================================            
            if GO_end in gold_list:
                gold = True
                silver = True
                #panther = True
            else:
                gold = False
                if GO_end in silver_list:
                    
                    silver = True
                    #panther = True
                else:
                    silver = False
                    
                    # try:
                    #     panther_list = Gene_panther_dic [gene_name]
                    
                    
                    #     if GO_end in panther_list:
                    #         panther = 1
                    #     else:
                    #         panther = 0
                    # except:
                    #     panther = 0    
            
            GO_path_element_list = list(set(GO_path_element_list))
            
            if len(GO_path_element_list) == 0:
                pass
            else:
                
                GO_end_series = csv_file_df[csv_file_df ['GO_index'] == GO_end]
                if GO_end_series ['score'].values == []:
                    pass
                else:
                    for val in GO_end_series ['score'].values:
                        pathdata.append([gene_name, file_names_[2], GO_end, val, ''.join(GO_end_series ['GO_names'].values), len(endpath_list), len(GO_path_element_list), ';'.join(GO_path_element_list), gold, silver])

pathdata_df = pd.DataFrame(pathdata, columns = ['Genes', 'GO_Class', 'End of GO Path', 'Cscore', 'GO Description', 'No. of GO Path', 'No. of GO elements', 'GO elements', 'Gold grade', 'Gold_Silver grade'])

pathdata_df.to_excel('GO_path_table.xlsx')





GO_info_dic = {}
data = []
protein_lists = []
ancestor_dic = {}
key_list = all_GO_dic.keys()

go_data = []
all_gold = []
all_gold_ancestor = []
all_silver = []
all_silver_ancestor = []

for key in key_list:
    keys = key.split('_')
    protein_lists.append(keys[0])
    if keys [-1] == 'BP':
        GO_class = 'Biological Process'
    elif keys [-1] == 'CC':
        GO_class = 'Cellular Component'
    elif keys [-1] == 'MF':
        GO_class = 'Molecular Function'
    
    runlist_dic2 = []
    GO_dic = all_GO_dic [key]
    try:
        runlist_dic2 = runlist_dic [keys[0]+'_'+keys[1]]
    except:
        print(keys[0]+'_'+keys[1] + " is not found.")
        pass
    try:
        gold_list = runlist_dic2 [1].split(';')
        #print('gold:', keys[0]+'_'+keys[1], len(gold_list), len(ancestor_dic [keys[0]+'_'+keys[1]][0]))
        #glod_list = gold_list + ancestor_dic [keys[0]+'_'+keys[1]][0]
        ancestor_gold_list = []
        for gold_ in gold_list:
            goldies = gold_.split(':')
            requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO%3A" + goldies[1] + "/ancestors?relations=is_a%2Cpart_of%2Coccurs_in%2Cregulates"
            r = requests.get(requestURL, headers={ "Accept" : "application/json"}, verify = False)
            
            if not r.ok:
              r.raise_for_status()
              sys.exit()
            
            responseBody = json.loads(r.text)
            result_dic = responseBody.get('results')[0]
            
            try:
                ancestor_gold_list = ancestor_gold_list + result_dic ['ancestors']
                ancestor_dic [gold_] = result_dic ['ancestors']
            except:
                pass            
        
    except:
        gold_list = []
        
    ancestor_gold_list = list(set(ancestor_gold_list))
    
    
    

    
        
        
    try:
        silver_list = runlist_dic2 [3].split(';')
        #print('silver:', keys[0]+'_'+keys[1], len(silver_list), len(ancestor_dic [keys[0]+'_'+keys[1]][0]))
        #silver_list = silver_list + ancestor_dic [keys[0]+'_'+keys[1]][1]
        ancestor_silver_list = []
        for gold_ in silver_list:
             goldies = gold_.split(':')
             requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO%3A" + goldies[1] + "/ancestors?relations=is_a%2Cpart_of%2Coccurs_in%2Cregulates"
             r = requests.get(requestURL, headers={ "Accept" : "application/json"}, verify = False)
            
             if not r.ok:
               r.raise_for_status()
               sys.exit()
            
             responseBody = json.loads(r.text)
             result_dic = responseBody.get('results')[0]
           
             try:
                 ancestor_silver_list = ancestor_silver_list + result_dic ['ancestors']
                 ancestor_dic [gold_] = result_dic ['ancestors']
             except:
                 pass            
    except:
        silver_list = []

    ancestor_silver_list = list(set(ancestor_silver_list))
    
    go_data.append([key, len(gold_list), len(ancestor_gold_list), len(silver_list), len(ancestor_silver_list)])
    
    all_gold = all_gold + gold_list
    all_gold_ancestor =  all_gold_ancestor + ancestor_gold_list
    all_silver = all_silver + silver_list
    all_silver_ancestor =  all_silver_ancestor + ancestor_silver_list
    
    
    gold_list = gold_list + ancestor_gold_list    
    silver_list = silver_list + ancestor_silver_list




    
    for GO_ in GO_dic.keys():
        
        
        
        
        gold = False
        silver = False
        #panther = 0
        GO_list = GO_dic[GO_]
        if GO_ in gold_list:
            gold = True
            silver = True
            #panther = 1
        else:
            gold = False
            if GO_ in silver_list:
                
                silver = True

        try:        
            ancestor = ancestor_dic [GO_]
        except:
            ancestor = []
        data.append([keys[0]+'_'+keys[1], GO_class, GO_, GO_list[1], GO_list[2], gold, silver])
        GO_info_dic [GO_] = GO_list[2]

data_df = pd.DataFrame (data, columns = ['Protein_Acc', 'GO_Class', 'GO_ID', 'Cscore', 'GO_description', 'Gold grade', 'Gold_Silver grade'])

data_df.to_excel('all_GO_lists_gold_silver.xlsx')
protein_lists = list(set(protein_lists))


data_gogo_df = pd.DataFrame(go_data, columns = ['Protein_Acc', 'no. of Gold', 'no. of Gold Ancestors', 'no. of Silver', 'no. of Silver Ancestors'])
data_gogo_df.to_excel('GO_Ancestor_no.xlsx')


print(len(all_gold), len(all_gold_ancestor), len(all_silver), len(all_silver_ancestor))

#sns.violinplot(x = 'Gold grade', y = 'Score', data = data_df, color = '0.8')
#sns.swarmplot(x = 'Gold grade', y = 'Score', data = data_df, size = 2)
#plt.show()
#
#sns.violinplot(x = 'Gold/Silver grade', y = 'Score', data = data_df, color = '0.8')
#sns.swarmplot(x = 'Gold/Silver grade', y = 'Score', data = data_df, size = 2)
#plt.show()
#
#sns.violinplot(x = 'Gold/Silver/Panther', y = 'Score', data = data_df, color = '0.8')
#sns.swarmplot(x = 'Gold/Silver/Panther', y = 'Score', data = data_df, size = 2)
#plt.show()



for grade in ['Gold grade', 'Gold_Silver grade']:
    #data_go_df = data_df[data_df['GO_Class'] == GO_class]
    fpr_list = []
    tpr_list = []
    roc_auc_list = []
    
    for GO_class in ['Biological Process', 'Cellular Component', 'Molecular Function']:
        data_go_df = data_df[data_df['GO_Class'] == GO_class]
        fpr, tpr, _ = roc_curve(data_go_df[grade], data_go_df['Cscore'])
        roc_auc = auc(fpr, tpr)
    
        fpr_list.append(fpr)
        tpr_list.append(tpr)
        roc_auc_list.append(roc_auc)

        mask = data_go_df.applymap(type) != bool
        d = {True: '1', False: '0'}

        data_go_df2 = data_go_df.where(mask, data_go_df.replace(d))
        print(GO_class, grade)        
        violinplot(data_go_df2, grade, GO_class) 
    # ROC_curve(fpr_list, tpr_list, roc_auc_list, ['Gold grade', 'Gold_Silver grade'], GO_class)
    
    # PR_curve (data_go_df, ['Gold grade', 'Gold_Silver grade'], GO_class)

    ROC_curve2(fpr_list, tpr_list, roc_auc_list, ['Biological Process', 'Cellular Component', 'Molecular Function'], grade, 'Unfiltered')
    
    PR_curve2 (data_df, ['Biological Process', 'Cellular Component', 'Molecular Function'], grade, 'Unfiltered')



           
        
#for GO_class in ['BP', 'CC', 'MF']:
for grade in ['Gold grade', 'Gold_Silver grade']:
    #data_go_df3 = pathdata_df[pathdata_df['GO Class'] == GO_class]
    
    fpr_list = []
    tpr_list = []
    roc_auc_list = []
    
  
    #for grade in ['Gold grade', 'Gold_Silver grade']:
    for GO_class in ['BP', 'CC', 'MF']:
        data_go_df3 = pathdata_df[pathdata_df['GO_Class'] == GO_class]
        fpr, tpr, _ = roc_curve(data_go_df3[grade], data_go_df3['Cscore'])
        roc_auc = auc(fpr, tpr)
    
        fpr_list.append(fpr)
        tpr_list.append(tpr)
        roc_auc_list.append(roc_auc)
        
        mask = data_go_df3.applymap(type) != bool
        d = {1: '1', 0: '0'}

        data_go_df4 = data_go_df3.where(mask, data_go_df3.replace(d))
        
        print(GO_class, grade)
        violinplot(data_go_df4, grade, GO_class)
        
#    ROC_curve(fpr_list, tpr_list, roc_auc_list, ['Gold grade', 'Gold_Silver grade'], GO_class + 'Cscore')            
    
#----------------------------------------------------------Precision Recall    
#    PR_curve (data_go_df3, ['Gold grade', 'Gold_Silver grade'], GO_class)    

    
    ROC_curve2(fpr_list, tpr_list, roc_auc_list, ['BP', 'CC', 'MF'], grade, 'Filtered')            
    
#----------------------------------------------------------Precision Recall    
    PR_curve2 (pathdata_df, ['BP', 'CC', 'MF'], grade, 'Filtered') 

        
        
            
        