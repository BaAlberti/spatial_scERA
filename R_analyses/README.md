# R analyses
This folder contains all the scripts where the scRNA-seq and PCR datasets are analysed with Seurat.

#### seurat_initial_analysis.R
This script is the core of the Seurat analysis of the scRNA-seq data only. In it you will find the major steps of scRNA-seq analysis conducted both on the full clustering and reduced clusteriing versions.
It also contains all the lines of codes used to generates the table resuming different information for each enhancer and each cluster
#### PCR addition.R
This script is starting from the output of the _trimming_PCR_ script. It will associate back the cells containing enhancers found from the PCR to the scRNA-seq dataset. 
It will generates the same tables as **seurat_initial_analysis.R** but either for PCR cells only or for common cells between PCR and scRNA-seq analyses.
### image_reproduction.R
This script is used to generate all the UMAPs with different informations regarding genes or enhancers.

### karaiskos_DGE.R
This last script is the one used to analyse Karaiskos dataset and obtain the list of 60 differentially expressed genes to look for enhancers at the beginning of the study.
