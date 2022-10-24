#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 15:21:06 2021

@author: zfd297
"""


import mte_main as mt
import pandas as pd
from sklearn.manifold import TSNE
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from idtxl import idtxl_io
#import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(rc={"lines.linewidth": 2}, palette  = "deep", style = "ticks")
from sklearn.metrics import precision_recall_curve, roc_curve, auc
from itertools import product, permutations, combinations, combinations_with_replacement
from tqdm import tqdm



def vis(df, p):
    tsne = TSNE(n_components=2,perplexity=p).fit_transform(df.T.values)
    tdf = pd.DataFrame(tsne, columns=['t-SNE 1', 't-SNE 2'],index=df.columns)
    fig, ax = plt.subplots(1,1,figsize=(5,5))
    ax.scatter(tsne[:,0], tsne[:,1], c = [float(col.split('_')[1]) for col in df.columns])
    #ax.axis('off')
    plt.yticks([])
    plt.xticks([])
    plt.tight_layout()
    plt.savefig("figure.pdf",transparent=True)
    
    return tsne


def convert_sif_to_rankedgrn(sif_name, fix_weight=False):
    
    
    load_sif = np.loadtxt(sif_name,dtype=object)
    if fix_weight == True:
        data_dic = {
            'Gene1': load_sif[:,0],
            'Gene2': load_sif[:,2],
            'EdgeWeight': [1 for i in range(load_sif.shape[0])]
            }
    else:

        data_dic = {
            'Gene1': load_sif[:,0],
            'Gene2': load_sif[:,2],
            'EdgeWeight': load_sif[:,1]
            }
    df = pd.DataFrame(data_dic)
    df.to_csv(sif_name[:-7]+"rankedEdges.csv",index=False,sep="\t")

def run_mte_1(exp_data,fdr_value,max_lag_sources=5,min_lag_sources=2):

    
    re1 = mt.run_mte_single_1(exp_data,fdr_value, max_lag_sources,min_lag_sources)
    
    adj_matrix = re1[1].get_adjacency_matrix(weights="max_te_lag", fdr=False)
    edge_list = adj_matrix.get_edge_list()
    gene_name = exp_data.index
    temp_dict = dict(zip([i for i in range(len(gene_name))],list(gene_name)))
    
    final_edge_list = []
    for i in edge_list:
        edge_tmp = (temp_dict[i[0]], temp_dict[i[1]],i[2])
        final_edge_list.append(edge_tmp)
    
    
    return re1




def get_te_value_max(idxtl_re,exp_data):
    
    connections = []
    ad_mat= idxtl_re.get_adjacency_matrix(weights="max_te_lag", fdr=False)
    edge_list = ad_mat.get_edge_list()
    edge_list_str = [str(i[0])+"|"+str(i[1]) for i in edge_list]
    target_list = np.array([i[1] for i in edge_list])
    target_list = np.unique(target_list)
    for i in target_list:
        tmp_dict = {}
        tmp_target = i
        print(tmp_target)
        select_sources = idxtl_re.get_single_target(tmp_target,fdr=False).selected_vars_sources
        tes = idxtl_re.get_single_target(tmp_target,fdr=False).te
        sources_list = np.array([j[0] for j in select_sources])
        
        unique_source = np.unique(sources_list)
        unique_source_index = np.argsort(unique_source)
        unique_source_sort = unique_source[unique_source_index]
        for k in range(len(unique_source_sort)):
            tmp_edge = (unique_source_sort[k],tmp_target,tes[k])
        
            connections.append(tmp_edge)
    
    connections_f = []
    for j in connections:
        if str(j[0])+"|"+str(j[1]) in edge_list_str:
            connections_f.append((j[0],j[1],j[2]))
        
    gene_name = exp_data.index
    temp_dict = dict(zip([i for i in range(len(gene_name))],list(gene_name)))
    
    final_edge_list = []
    for i in connections_f:
        edge_tmp = (temp_dict[i[0]], temp_dict[i[1]],i[2])
        final_edge_list.append(edge_tmp)
        
        
    return final_edge_list





def compute_roc_single(predDF, trueEdgesDF):
    

    
    possibleEdges = list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                             r = 2))
    TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}
    PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}
    
    for key in TrueEdgeDict.keys():
        if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == key.split('|')[0]) &
               (trueEdgesDF['Gene2'] == key.split('|')[1])])>0:
                TrueEdgeDict[key] = 1
                
    predEdgeDF = predDF     
    for key in PredEdgeDict.keys():
        subDF = predEdgeDF.loc[(predEdgeDF['Gene1'] == key.split('|')[0]) &
                           (predEdgeDF['Gene2'] == key.split('|')[1])]
        if len(subDF)>0:
            PredEdgeDict[key] = np.abs(subDF.EdgeWeight.values[0])
    
    key  = 'g2|g1'
    subDF = predEdgeDF.loc[(predEdgeDF['Gene1'] == key.split('|')[0]) &
                           (predEdgeDF['Gene2'] == key.split('|')[1])]
    
    
    
    outDF = pd.DataFrame([TrueEdgeDict,PredEdgeDict]).T
    outDF.columns = ['TrueEdges','PredEdges']
    
    
    fpr, tpr, thresholds = roc_curve(y_true=outDF['TrueEdges'],
                                         y_score=outDF['PredEdges'], pos_label=1)

    
    return fpr, tpr, auc(fpr, tpr),thresholds