import novosparc
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy import stats
import os
import scanpy as sc
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.metrics import accuracy_score
import random

dataset=sc.read_h5ad('dataset.h5ad')
atlas=sc.read_h5ad('atlas.h5ad')
atlas_genes = atlas.var.index.tolist()
gene_names = dataset.var.index.tolist()
num_cells, num_genes = dataset.shape

atlas_dir = '../../novosparc/novosparc/datasets/bdtnp/'
target_space_path = os.path.join(atlas_dir, 'geometry.txt')
locations = pd.read_csv(target_space_path, sep=' ')
num_locations = 3039
locations_apriori = locations[:num_locations][['xcoord', 'ycoord', 'zcoord']].values

result_matrix=pd.DataFrame(columns=['atlas gene','spearman score'])
random_matrix=pd.DataFrame(columns=['atlas gene','iter','spearman score'])
for gene in atlas_genes:
	print(gene)

	full_list=atlas_genes.copy()
	full_list.remove(gene)
	atlas_reduced=atlas[:,full_list]

	tissue = novosparc.cm.Tissue(dataset=dataset, locations=locations_apriori)

	num_neighbors_s = 3
	num_neighbors_t = 5
	markers = list(set(full_list).intersection(gene_names))
	atlas_matrix = atlas_reduced.to_df()[markers].values
	markers_idx = pd.DataFrame({'markers_idx': np.arange(num_genes)}, index=gene_names)
	markers_to_use = np.concatenate(markers_idx.loc[markers].values)
	tissue.setup_reconstruction(atlas_matrix=atlas_matrix, 
	                            markers_to_use=markers_to_use, 
	                            num_neighbors_s=num_neighbors_s, 
	                            num_neighbors_t=num_neighbors_t)

	alpha_linear = 0.35
	epsilon = 5e-3
	tissue.reconstruct(alpha_linear=alpha_linear, epsilon=epsilon)
	sdge = tissue.sdge 
	dataset_reconst = sc.AnnData(pd.DataFrame(sdge.T, columns=gene_names))
	dataset_reconst.obsm['spatial'] = locations_apriori

	y_test=np.array(sc.get.obs_df(dataset_reconst,keys=gene))
	y_test= (y_test-np.min(y_test))/(np.max(y_test)-np.min(y_test))
	y_truth=np.array(sc.get.obs_df(atlas,keys=gene))
	y_truth= (y_truth-np.min(y_truth))/(np.max(y_truth)-np.min(y_truth))

	spearman_score=stats.spearmanr(y_test,y_truth)
	result_matrix.loc[len(result_matrix)]=[gene,spearman_score]

	for i in range(100):
		y_test_random=y_test
		random.Random(i).shuffle(y_test_random) # Random(i) permet d'avoir la même seed à chaque fois qu'on fera retourner le script
		spearman_score=stats.spearmanr(y_test_random,y_truth)
		random_matrix.loc[len(random_matrix)]=[gene,i+1,spearman_score]

result_matrix.to_csv("atlas_spearman_correlation.csv")
random_matrix.to_csv("random_spearman_correlation.csv")

s=r'\b(\d\.\d+)\b'
p=r'\b(\d\.\d+e?-?\d+|0\.0)\b'
result_matrix['actual score']=result_matrix['spearman score'].str.extract(s,expand=False).astype(float)
result_matrix['actual pvalue']=result_matrix['spearman score'].str.findall(p).str[-1].astype(float)

p=r'\b(\-?\d\.\d+e?-?\d+|0\.0)\b'
random_matrix['actual score']=random_matrix['spearman score'].str.findall(p).str[0]
random_matrix['actual score']=random_matrix['actual score'].astype(float)
random_matrix['actual pvalue']=random_matrix['spearman score'].str.findall(p).str[-1]
random_matrix['actual pvalue']=random_matrix['actual pvalue'].astype(float)

mydict={"non random":result_matrix["actual score"],"random":random_matrix['actual score']}
fig, ax = plt.subplots()
ax.boxplot(mydict.values(),widths=0.5)
ax.set_xticklabels(mydict.keys(),fontsize=15)
ax.set_ylabel('spearman correlation', fontsize=15)
plt.savefig('boxplot_spearman_corr.png')

print(stats.kstest(spearman_table["actual score"],spearman_random["actual score"]))

result_matrix.nlargest(5,columns="actual score",keep="first").to_csv("highest_spearman_score.csv")
result_matrix.nlowest(5,columns="actual score",keep="first").to_csv("lowest_spearman_score.csv")
