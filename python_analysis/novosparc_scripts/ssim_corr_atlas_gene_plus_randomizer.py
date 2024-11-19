"""
This script aims at performing a leave-one-out cross-validation.
The idea is to remove a gene in the atlas, performing the reconstruction, and observed how it is reconstructed compared to its initial pattern in the atlas with a spearman correlation test.
With this script, each of the 84 genes of the atlas are tested.
To perform a statistical comparison, we compare the spearman score of the reconstruction against a "random" reconstruction where the values are kept the same but ranodmized across the embryo
This comparison is saved as a boxplot at the end of the script.
"""

# Loading packages
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

# Defining ssim score
def ssim_1d_manual(array1, array2):
    # Constants to stabilize the division
    C1 = (0.01) ** 2  # Assuming data range of 255 for 8-bit images; adjust as needed
    C2 = (0.03) ** 2
    
    # Mean of each array
    mu_x = array1.mean()
    mu_y = array2.mean()
    
    # Variance and covariance
    sigma_x = array1.var()
    sigma_y = array2.var()
    sigma_xy = np.cov(array1, array2)[0, 1]
    
    # Compute SSIM based on the formula
    ssim_score = ((2 * mu_x * mu_y + C1) * (2 * sigma_xy + C2)) / ((mu_x**2 + mu_y**2 + C1) * (sigma_x + sigma_y + C2))
    return ssim_score

# Loading dataset
dataset=sc.read_h5ad('dataset.h5ad')
atlas=sc.read_h5ad('atlas.h5ad')
atlas_genes = atlas.var.index.tolist()
gene_names = dataset.var.index.tolist()
num_cells, num_genes = dataset.shape

# Loading original atlas
atlas_dir = '../../data/novosparc_atlas/'
target_space_path = os.path.join(atlas_dir, 'geometry.txt')
locations = pd.read_csv(target_space_path, sep=' ')
num_locations = 3039
locations_apriori = locations[:num_locations][['xcoord', 'ycoord', 'zcoord']].values

# Performing the leave-one-out cross-validation
result_matrix=pd.DataFrame(columns=['atlas gene','ssim score'])
random_matrix=pd.DataFrame(columns=['atlas gene','iter','ssim score'])
for gene in atlas_genes:
	print(gene)

    # removing one gene from the atlas
	full_list=atlas_genes.copy()
	full_list.remove(gene)
	atlas_reduced=atlas[:,full_list]

    # Perform the reconstruction without it
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

    # normalizing the values observed after the reconstruction and in the atlas
	A=atlas[:,gene].X
	A=(A-np.min(A))/(np.max(A)-np.min(A))
	A=A.reshape(1,-1)
	B=dataset_reconst[:,gene].X
	B=(B-np.min(B))/(np.max(B)-np.min(B))
	B=B.reshape(1,-1)

	# Perform the spearman correlation between the reconstruction and the atlas
	ssim_score=ssim_1d_manual(A,B)
	result_matrix.loc[len(result_matrix)]=[gene,ssim_score]

    # Perform a "random" reconstruction 100 times and save the spearman score in another table. 
	for i in range(100):
		B_random=B
		random.Random(i).shuffle(B_random[0]) # Random(i) permet d'avoir la même seed à chaque fois qu'on fera retourner le script
		ssim_score=ssim_1d_manual(B_random,A)
		random_matrix.loc[len(random_matrix)]=[gene,i+1,ssim_score]

result_matrix.to_csv("atlas_ssim_correlation.csv")
random_matrix.to_csv("random_ssim_correlation.csv")

# Make a boxplot comparing the corss-validation to a random reconstructions
mydict={"non random":result_matrix["ssim score"],"random":random_matrix['ssim score']}
fig, ax = plt.subplots()
ax.boxplot(mydict.values(),widths=0.5)
ax.set_xticklabels(mydict.keys(),fontsize=15)
ax.set_ylabel('ssim correlation', fontsize=15)
plt.savefig('boxplot_ssim_corr.png')

print(stats.kstest(spearman_table["ssim score"],spearman_random["ssim score"]))

# Save the 5 genes with the highest and lowest spearman scores
result_matrix.nlargest(5,columns="ssim score",keep="first").to_csv("highest_ssim_score.csv")
result_matrix.nlowest(5,columns="ssim score",keep="first").to_csv("lowest_ssim_score.csv")
