#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 11:35:41 2021

@author: zfd297
"""

import mte_main as mt
import pandas as pd
from sklearn.manifold import TSNE
import matplotlib
import matplotlib.pyplot as plt

def vis(df, p):
    tsne = TSNE(n_components=2,perplexity=p).fit_transform(df.T.values)
    tdf = pd.DataFrame(tsne, columns=['t-SNE 1', 't-SNE 2'],index=df.columns)
    fig, ax = plt.subplots(1,1,figsize=(5,5))
    ax.scatter(tsne[:,0], tsne[:,1], c = [float(col.split('_')[1]) for col in df.columns])
    ax.axis('off')
    plt.tight_layout()

df=pd.read_csv("/Users/zfd297/workspace/verify_mte/0815_newdata/dyn-
               
               furcating.txt/ExpressionData.csv",index_col=0)
vis(df, 400)
                                                                               
       
                                                                                         


exp_data = pd.read_csv("/Users/zfd297/workspace/software/BoolODE-0.1/dyn-consecutive-bifurcating/ExpressionData.csv",index_col=0)
pseudo_time_csv = pd.read_csv("/Users/zfd297/workspace/software/BoolODE-0.1/dyn-consecutive-bifurcating/PseudoTime.csv",index_col=0)


re1 = mt.run_mte_single(exp_data, pseudo_time_csv,fdr_value=False , max_lag_sources=5,min_lag_sources=2,save=True)
adj_matrix = re1[1].get_adjacency_matrix(weights="max_te_lag", fdr=False)
edge_list = adj_matrix.get_edge_list()
gene_name = exp_data.index
temp_dict = dict(zip([i for i in range(len(gene_name))],list(gene_name)))

final_edge_list = []
for i in edge_list:
    edge_tmp = (temp_dict[i[0]], temp_dict[i[1]],i[2])
    final_edge_list.append(edge_tmp)
        
f = open('./conse_bifurcating_mte.sif','w')
for k in final_edge_list:
    write_str = str(k[0]) +" "+ str(i[2]) + " "+str(k[1]) + "\n"
    f.write(write_str)
f.close()
        

    
exp_data = pd.read_csv("/Users/zfd297/workspace/verify_mte/dyn_linear_exp/dyn-linearExpressionData.csv",index_col=0)
pseudo_time_csv = pd.read_csv("/Users/zfd297/workspace/verify_mte/dyn_linear_exp/dyn-linearPseudoTime.csv",index_col=0)


def robust_run(exp_data,pseudo_time_csv, max_l, min_l,save_path):


    re1 = mt.run_mte_single(exp_data, pseudo_time_csv,fdr_value=False , max_lag_sources=max_l,min_lag_sources=min_l,save=False)
    adj_matrix = re1[1].get_adjacency_matrix(weights="max_te_lag", fdr=False)
    edge_list = adj_matrix.get_edge_list()
    gene_name = exp_data.index
    temp_dict = dict(zip([i for i in range(len(gene_name))],list(gene_name)))
    
    final_edge_list = []
    for i in edge_list:
        edge_tmp = (temp_dict[i[0]], temp_dict[i[1]],i[2])
        final_edge_list.append(edge_tmp)
            
    f = open(save_path,'w')
    for k in final_edge_list:
        write_str = str(k[0]) +" "+ str(k[2]) + " "+str(k[1]) + "\n"
        f.write(write_str)
    f.close()


exp_data = pd.read_csv("/Users/zfd297/workspace/verify_mte/dyn_linear_exp/dyn-linearExpressionData.csv",index_col=0)
pseudo_time_csv = pd.read_csv("/Users/zfd297/workspace/verify_mte/dyn_linear_exp/dyn-linearPseudoTime.csv",index_col=0)

parameter_list = [(5,2),(10,2),(10,5),(15,2),(15,5),(15,10),(20,2),(20,5),
                  (20,10),(20,15),(25,2),(25,5), (25,10),(25,15),(25,20), (30,2),
                  (30,5),(30,10),(30,15),(30,20),(30,25)]

for k in parameter_list:
    for i in range(3):
        save_path = "./robust_linear/"+str(k[0])+"_"+str(k[1])+"_round_" +str(i)+".sif"
        robust_run(exp_data,pseudo_time_csv, k[0], k[1], save_path)





exp_data = pd.read_csv("/Users/zfd297/workspace/verify_mte/cycle/cycleExpressionData.csv",index_col=0)
pseudo_time_csv = pd.read_csv("/Users/zfd297/workspace/verify_mte/cycle/cyclePseudoTime.csv",index_col=0)

parameter_list = [(5,2),(10,2),(10,5),(15,2),(15,5),(15,10),(20,2),(20,5),
                  (20,10),(20,15),(25,2),(25,5), (25,10),(25,15),(25,20), (30,2),
                  (30,5),(30,10),(30,15),(30,20),(30,25)]
for k in parameter_list:
    for i in range(3):
        save_path = "./robust_cycle/"+str(k[0])+"_"+str(k[1])+"_round_" +str(i)+".sif"
        robust_run(exp_data,pseudo_time_csv,k[0], k[1], save_path)



exp_data = pd.read_csv("/Users/zfd297/workspace/verify_mte/dyn_reduce/dyn-converging_reduceExpressionData.csv",index_col=0)
pseudo_time_csv = pd.read_csv("/Users/zfd297/workspace/verify_mte/dyn_reduce/dyn-converging_reducePseudoTime.csv",index_col=0)

parameter_list = [(5,2),(10,2),(10,5),(15,2),(15,5),(15,10),(20,2),(20,5),
                  (20,10),(20,15),(25,2),(25,5), (25,10),(25,15),(25,20), (30,2),
                  (30,5),(30,10),(30,15),(30,20),(30,25)]
for k in parameter_list:
    for i in range(3):
        save_path = "./robust_bifurcating/"+str(k[0])+"_"+str(k[1])+"_round_" +str(i)+".sif"
        robust_run(exp_data,pseudo_time_csv,k[0], k[1], save_path)






######

exp_file = pd.read_csv("./exp_data_pancera/ExpressionData.csv",index_col=0)
ps_time = pd.read_csv("./exp_data_pancera/PseudoTime.csv",index_col=0)["PseudoTime"]
cell_sort_index = np.argsort(ps_time)
cell_sort_index = np.array(cell_sort_index)
exp_file_sort = exp_file.iloc[:,cell_sort_index]

mte_file = np.loadtxt("./pancrea/PIDC_cutoff_34.sif",dtype=object)

pos_neg = mt.find_negtive(mte_file, exp_file)

pos = mte_file[pos_neg[0],:]
neg = mte_file[pos_neg[1],:]

method = "PIDC"
f = open('./tmp_posi_neg_sif/'+ method+"_.sif",'w')
for k in pos:
    write_str = str(k[0]) +" "+ "+" + " "+str(k[2]) + "\n"
    f.write(write_str)
for k in neg:
    write_str = str(k[0]) +" "+ "-" + " "+str(k[2]) + "\n"
    f.write(write_str)
f.close()





def get_true_grn(filename):
    
    true_sif = []
    df = pd.read_csv(filename)
    targets = np.array(df["Gene2"])
    regulators = np.array(df["Gene1"])
    rules = np.array(df["Type"])
    
    f = open('./tmp_bool_true.sif','w')
    
    for i in range(len(targets)):
        
        if rules[i] == "+":
            
            tmp = (regulators[i], targets[i])
            true_sif.append(tmp)
            f.write(regulators[i] + " " + "1" + " " + targets[i]+"\n")
        else:
            tmp = (regulators[i], targets[i])
            true_sif.append(tmp)
            f.write(regulators[i] + " " + "0" + " " + targets[i]+"\n")
            
         
    f.close()
        
    return true_sif
    
bool_grn = get_true_grn("/Users/zfd297/workspace/verify_mte/cycle/cyclerefNetwork.csv")

import numpy as np
import os

def convert_csv_to_sif(sample_name):


    filePath = '/Users/zfd297/workspace/software/Beeline/outputs/example/' + sample_name + "/"
    for i in os.listdir(filePath)[1:]:
        
        read_csv_path = filePath+i+"/"+"rankedEdges.csv"
        df = pd.read_csv(read_csv_path,sep="\t")
        targets = np.array(df["Gene2"])
        regulators = np.array(df["Gene1"])
        rules = np.array(df["EdgeWeight"])
        
        write_path = "/Users/zfd297/workspace/verify_mte/" + sample_name +"/" +i +".sif"
        f = open(write_path,'w')
        
        for i in range(len(regulators)):
            f.write(regulators[i] + " " + str(rules[i]) + " " + targets[i]+"\n")
            
        f.close()
        
        
convert_csv_to_sif("dyn_linear")

######run junil data

exp_data = pd.read_csv("/Users/zfd297/workspace/verify_mte/junil_data_mte/sampleDEGlogRatio05autophagyGluconeoAmpk_4wishbone.csv",index_col=0)
pseudo_time_csv = pd.read_csv("/Users/zfd297/workspace/verify_mte/junil_data_mte/dyn-converging_reducePseudoTime.csv",index_col=0)

cell_select_brach12 = np.loadtxt("/Users/zfd297/workspace/verify_mte/junil_data_mte/cell_select_branch1_2.txt")
cell_select_brach12_index = np.where(cell_select_brach12==1)[0]
branch12_pstime = np.loadtxt("/Users/zfd297/workspace/verify_mte/junil_data_mte/pseudotime_branch1_2.txt")
branch12_pstime = branch12_pstime[cell_select_brach12_index]
branch12_pstime_sort_index = np.argsort(branch12_pstime)

exp_branch12 = exp_data.iloc[cell_select_brach12_index,:]
exp_branch12_sort = exp_branch12.iloc[branch12_pstime_sort_index,:]

tf_genes = np.loadtxt("/Users/zfd297/workspace/verify_mte/junil_data_mte/TFgene_names2",dtype=object)

exp_branch12_sort_t = exp_branch12_sort.T
genes = exp_branch12_sort.index

exp_branch12_sort_tf_genes = exp_branch12_sort_t.loc[tf_genes,:]

re1 = mt.run_mte_single_1(exp_branch12_sort_tf_genes, fdr_value=False , max_lag_sources=5,min_lag_sources=2,save=True)


import umap

reducer = umap.UMAP()
embedding = reducer.fit_transform(exp_data.values)
embedding.shape

embedding_branch12 = embedding[cell_select_brach12_index,:]
embedding_branch12_sort = embedding_branch12[branch12_pstime_sort_index,:]

plt.scatter(embedding_branch12_sort[:,0], embedding_branch12_sort[:,1],c=branch12_pstime[branch12_pstime_sort_index)
                                                                                         
                                                                                         
                                                                                         
                                                                                         
cell_select_brach13 = np.loadtxt("/Users/zfd297/workspace/verify_mte/junil_data_mte/cell_select_branch1_3.txt")
cell_select_brach13_index = np.where(cell_select_brach13==1)[0]
branch13_pstime = np.loadtxt("/Users/zfd297/workspace/verify_mte/junil_data_mte/pseudotime_branch1_3.txt")
branch13_pstime = branch13_pstime[cell_select_brach13_index]
branch13_pstime_sort_index = np.argsort(branch13_pstime)

exp_branch13 = exp_data.iloc[cell_select_brach13_index,:]
exp_branch13_sort = exp_branch13.iloc[branch13_pstime_sort_index,:]

embed = pd.read_csv("/Users/zfd297/workspace/verify_mte/Ndufa8_q001/cell_trajectory2D.txt", sep=" ",index_col=0)
embed_mat = embed.values[:,1:]
embed_mat = embed_mat.astype(float)


ps_time = np.loadtxt("/Users/zfd297/workspace/verify_mte/Ndufa8_q001/trajectory_GGGAATGAGGCTAGCA-1.txt")
ps_time = np.array(ps_time)
ps_time_sort_index = np.argsort(ps_time)
ps_time_sort = np.array(ps_time[ps_time_sort_index])

exp_data_sort = exp_data.iloc[ps_time_sort_index,:]

sample_index = [i for i in range(0, len(ps_time_sort), 3)]
sample_index = np.array(sample_index)

sample_pstime = ps_time[sample_index]

sample_exp = exp_data_sort.iloc[sample_index,:]

sample_exp_tf = sample_exp.loc[:,tf_genes]
sample_exp_tf_t = sample_exp_tf.T

sample_exp_tf_t.to_csv("6999_sample_exp.csv")


exp_data_sort_t = exp_data_sort.T
exp_data_sort_t_tf = exp_data_sort_t.loc[tf_genes, :]
exp_data_sort_t_tf.to_csv("exp_data_sort.csv")

sample_exp_tf_t_1 = sample_exp_tf_t.iloc[1:4,:]


re1 = mt.run_mte_single_1(sample_exp_tf_t_1,fdr_value=False , max_lag_sources=5,min_lag_sources=2,save=False)


with open('./the_first_pickle.pickle','rb') as p:
    p.write()   #将列表t保存起来

from idtxl import idtxl_io
idtxl_io.save_pickle(ll,"nace")

sample_6999 = idtxl_io.load_pickle("6999_samle_idtxl")
sample_2100 = idtxl_io.load_pickle("sample_exp_idtxl")


adj_matrix = sample_6999.get_adjacency_matrix(weights="max_te_lag", fdr=0.01)
edge_list = adj_matrix.get_edge_list()
gene_name = exp_data.index
temp_dict = dict(zip([i for i in range(len(gene_name))],list(gene_name)))


adj_matrix = sample_6999.get_adjacency_matrix(weights="max_p_lag", fdr=False)

adata = sca.read("/Users/zfd297/workspace/verify_mte/dyngen_data/100tf.h5ad")

exp_mat = adata.X.A.T

df = pd.DataFrame(exp_mat, index = adata.var.index, columns = adata.obs.index)

sim_time = np.array(adata.obs["sim_time"])
sim_time_index = np.argsort(sim_time)

df_sort = df.iloc[:, sim_time_index]

df_sort.to_csv("./other_data/100tf.csv")

plt.scatter(adata.obsm["dimred"][:,0], adata.obsm["dimred"][:,1],c=sim_time)