import pandas as pd
from idtxl.multivariate_te import MultivariateTE
from idtxl.data import Data
from idtxl.visualise_graph import plot_network
import matplotlib.pyplot as plt
import numpy as np
from idtxl.multivariate_mi import MultivariateMI


def run_mte_single(exp_data, pseudo_time_csv,fdr_value,max_lag_sources,min_lag_sources,save=False):
    
    exp_mat = exp_data.values
    pseudo_time = pseudo_time_csv["Time"]
    pseudo_time_sort_index = np.argsort(pseudo_time)
    pseudo_time_sort = pseudo_time[pseudo_time_sort_index]
    
    pseudo_time_genes = np.array(pseudo_time_sort.index)
    
    exp_data_sort = exp_data.loc[:,pseudo_time_genes]
    
    
    exp_sort_mat = exp_data_sort.values
    
    
    data_2 = Data(exp_data_sort, dim_order='ps',normalise=True) 
    
    network_analysis = MultivariateTE()
    settings = {'cmi_estimator': 'JidtGaussianCMI',
                'n_perm_max_stat': 200,
        'n_perm_min_stat': 200,
        'n_perm_omnibus': 500,
        'n_perm_max_seq': 500,
        'max_lag_sources': max_lag_sources,
        'min_lag_sources': min_lag_sources}
    
    
    
    # c) Run analysis
    results = network_analysis.analyse_network(settings=settings, data=data_2)
    
    results.print_edge_list(weights='max_te_lag', fdr=False)
    plot_network(results=results, weights='max_te_lag', fdr=False)
    plt.show()
    
    
    adj_matrix = results.get_adjacency_matrix(weights="max_te_lag", fdr=fdr_value)
    edge_list = adj_matrix.get_edge_list()
    
    
    new_edge_list = []
    for i in edge_list:
        
        new_edge_list.append([i[0], i[1]])
        
    
    gene_name = exp_data_sort.index
    temp_dict = dict(zip([i for i in range(len(gene_name))],list(gene_name)))
    
    final_edge_list = []
    for i in new_edge_list:
        edge_tmp = (temp_dict[i[0]], temp_dict[i[1]])
        final_edge_list.append(edge_tmp)
    
    
    if save == True:
        f = open('./tmp_mte.sif','w')
        for k in final_edge_list:
            write_str = str(k[0]) +" "+ "1" + " "+str(k[1]) + "\n"
            f.write(write_str)
        f.close()

        
    return final_edge_list, results


def run_mte_single_1(exp_data,fdr_value,max_lag_sources,min_lag_sources,save=False):
    
    exp_mat = exp_data.values
    exp_data_sort = exp_data
    exp_sort_mat = exp_data_sort.values
    
    
    data_2 = Data(exp_data_sort, dim_order='ps') 
    
    network_analysis = MultivariateTE()
    settings = {'cmi_estimator': 'JidtGaussianCMI',
                'n_perm_max_stat': 200,
        'n_perm_min_stat': 200,
        'n_perm_omnibus': 500,
        'n_perm_max_seq': 500,
        'max_lag_sources': max_lag_sources,
        'min_lag_sources': min_lag_sources}
    
    
    
    # c) Run analysis
    results = network_analysis.analyse_network(settings=settings, data=data_2)
    
    results.print_edge_list(weights='max_te_lag', fdr=False)
    plot_network(results=results, weights='max_te_lag', fdr=False)
    plt.show()
    
    
    adj_matrix = results.get_adjacency_matrix(weights="max_te_lag", fdr=fdr_value)
    edge_list = adj_matrix.get_edge_list()
    
    
    new_edge_list = []
    for i in edge_list:
        
        new_edge_list.append([i[0], i[1]])
        
    
    gene_name = exp_data_sort.index
    temp_dict = dict(zip([i for i in range(len(gene_name))],list(gene_name)))
    
    final_edge_list = []
    for i in new_edge_list:
        edge_tmp = (temp_dict[i[0]], temp_dict[i[1]])
        final_edge_list.append(edge_tmp)
    
    
    if save == True:
        f = open('./tmp_mte.sif','w')
        for k in final_edge_list:
            write_str = str(k[0]) +" "+ "1" + " "+str(k[1]) + "\n"
            f.write(write_str)
        f.close()

        
    return final_edge_list, results



####find negtive
# mte = np.loadtxt("/Users/zfd297/workspace/software/BoolODE-0.1/bifurcating_mte_1.sif",dtype=object)
# re = []
# for i in mte:
#     tmp1 = pd.Series(exp_data_sort.loc[i[0]])
#     tmp2 = pd.Series(exp_data_sort.loc[i[2]])
    
#     corr=tmp1.corr(tmp2)

#     re.append(corr)
    
    
# re = np.array(re)
# index2 = np.where(re<0)[0]

# cc = mte[index2,:]



#true_grn = np.loadtxt("/Users/zfd297/workspace/software/BoolODE-0.1/data/long-linear.txt",delimiter="\t",dtype=object)


import numpy as np
import pandas as pd
import re


def convert_boolodegrn_to_sif(filename):
    
    sif_ = []
    
    grn = np.loadtxt(filename, delimiter="\t",dtype=object)
    
    activator_list = []
    repressor_list = []
    for i in grn[1:]:
        
        tmp_target = i[0]
        tmp_regulators = i[1]
        
        if "not" in tmp_regulators:
            
            tmp_regulators = tmp_regulators.split("not")
            
            tmp_activators = tmp_regulators[0]
            tmp_repressors = tmp_regulators[1]
            
            pattern = re.compile(r'g\d.') 
            
            activator = pattern.findall(tmp_activators)
            repressor = pattern.findall(tmp_repressors)
            
            activator_list.append(activator)
            repressor_list.append(repressor)
            
            for k in activator:
                tmp_relation = k+" "+tmp_target
                sif_.append(tmp_relation)
            for k in repressor:
               tmp_relation = k+" "+tmp_target
               sif_.append(tmp_relation)

        else:
            
            pattern = re.compile(r'g\d.') 
            activator = pattern.findall(tmp_regulators)
            activator_list.append(activator)
            for k in activator:
                tmp_relation = k+" "+tmp_target
                sif_.append(tmp_relation)
            
    

        
    return sif_


def find_negtive(mte, exp_final):
    
    re = []
    for i in mte:
        tmp1 = pd.Series(exp_final.loc[i[0]])
        tmp2 = pd.Series(exp_final.loc[i[2]])
        
        corr=tmp1.corr(tmp2)
        
        
        re.append(corr)
        
    
    re = np.array(re)
    index_neg = np.where(re<0)[0]
    index_pos = np.where(re>0)[0]
    return index_pos,index_neg,re




# true_grn = convert_boolodegrn_to_sif("/Users/zfd297/workspace/software/BoolODE-0.1/data/dyn-linear-long.txt")

# true_grn_connect = []
# for i in true_grn:
#     tmp = i.split(" ")
#     tmp_append = tmp[0]+tmp[-1]
#     true_grn_connect.append(tmp_append)

# count=0
# for i in final_edge_list:
    
#     tmp = i[0]  +i[1]
#     if tmp in true_grn_connect:
#         print(tmp)
#         count+=1
        
        
####plot heatmap
# plt.figure(None,(20,10))
# plt.imshow(exp_sort_mat, cmap='hot', interpolation='nearest')
# plt.show()

# import seaborn as sns

# ax = sns.heatmap(exp_data_sort)
# plt.show()


# from sklearn.metrics import r2_score, mean_absolute_error, median_absolute_error



# rolling_mean = exp_data_sort.rolling(window=1).mean()
# ax = sns.heatmap(rolling_mean)
# plt.show()

# plt.plot(rolling_mean.loc["g8",:])
# plt.show()