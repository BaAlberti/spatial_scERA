# novoSpaRc analyses
This folder contains all the scripts where the scRNA-seq dataset was projected onto a virtual embryo with novoSpaRc.

#### novosparc_initial_analysis.ipynb
This script is the core of the novoSpaRc analysis. In it you will find the steps to generate the reconstruction from the scRNA-seq data.
The remaining of the script is calling the **home_functions.py** file in order to perform 3D visualisation of enhancer spatial activity, cluster projections, and gene expressions.
#### ssim_corr_atlas_gene_plus_randomization.py
This script was generated to perform the leave-one-out cross-validation. For each gene in the atlas, a loop will first remove this gene from the atlas, perform the reconstruction, and then compare with (structural similarity index measure) ssim score the predicted gene activity against its initial activity in the atlas. <br>
It will also perform a 100 randomization of the predictions in order to generate a "random distribution" for each gene in the atlas.
#### leave_one_out_cross_validation.ipynb
From the output of the last python script, this one was used to visualy inspect the reconstruction of the genes with the highest and lowest spearman scores. 
### reconstruction_with_novosparc_tutorial_dataset.ipynb
Notebook performing the same pipeline as in the novosparc_initial_analysis script with Karaiskos dataset as scRNA-seq input. This notebook was made to observe if the strippy pattern observed in our reconstruction was comming from our data or the atlas itself.
