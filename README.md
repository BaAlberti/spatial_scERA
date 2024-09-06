# spatial_scERA

### Goal
We developped the spatial_scERA method to predict the spatial activity of _Drosophila melanogaster_ enhancers during embryogenesis, thanks to scRNA-seq and novoSpaRc tools.
Our first experiment caracterised the spatial activity of 25 enhancers in a stage 6 embryo.
This git repository contains all the scripts used to carry out the analyses and to generats the figures found in the paper: <link to the paper>

### Architecture
The repository is composed of several folders:
- **envs**: contains the _conda yml_ file and _Renv files_ in order to load and use the same packages/libraries used in the analyses. 
- **bash**: This folder contains the scripts used to generate the preprocessed data (mapping scRNA-seq datasets with cellranger and trimming PCR datasets with cutadapt).  
- **data**: All the intermediate files of the scripts will be redirected to the folders inside. It is also the location where you should download the processed files from ArrayExpress in order to run the scripts.
- **R_analyses**: contains all the scripts made in R. More info on the scripts themselves can be found in the dedicated README
- **python_analyses**: This folder is separataed in two.
  - The first folder **trimming_PCR** has the necessary scripts to list all the enhancer present in the cells from the PCR data.
  - The second folder **novosparc_scripts** has all the scripts where novoSpaRc is involved. The **home made functions** coded to improve novoSpaRc are in the _home_functions.py_ file. More details on each scripts can be found in the dedicated README

### Global pipeline
1. Load the processed data at E-MTAB-14447 (scRNA-seq), and E-MTAB-14445 (targeted PCR amplification) + the custome genome with the doi: 10.5281/zenodo.13711749 if you want to run cellranger
2. Run the seurat_initial_analysis.R to analyse the scRNA-seq datasets
3. Run the trimming_PCR script to generate the cell/enhancer tables from PCR data
4. Run the novosparc_initial_analysis notebook to generate the predictions of each enhancer spatial activity.

### Note
The other scripts are usefull to generate images or perform supplementary analysis. You can learn more in the dedicated README and in the commentaries in the scripts.


  
