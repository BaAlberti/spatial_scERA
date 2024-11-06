# Description of the different folders and datatypes that you can find or load in this data folder

This folder is gathering every data generated along the analysis, or is storing data to run them. <br>
For several folders you need to download files, the links are provided in this README. As all the scripts are getting their inputs/outputs from here, it is encouraged to follow the structure of this data folder. <br>

### custom_genome
**Download at**: Lien vers Zenodo avec le génome de référence. <br>
This folder contains the reference genome to use with CellRanger if you want to map the reads from the fastq files you can download at E-MTAB-14447
### cutadapt_output
Stores the cutadapt processed fastqs files from the preprocessing script **cutadapt.sh**. This folder contains the trimmed fastqs from PCR sequencing. They are used later in the **clean_trimming_script_novopaper.ipynb** script.
### novosparc_atlas
Stores the atlas used with novoSpaRc for this project.
### novosparc_dataset
**Download at**: https://github.com/rajewsky-lab/novosparc/tree/master/novosparc/datasets/drosophila_scRNAseq <br>
This folder contains the dataset used in the tutorial of novosparc.
### PCR_sequencing_output
**Download at**: Lien vers arrayexpress (E-MTAB-14445) <br>
This folder contains the three sequenced library of the PCR targeted amplification.
### scRNA-seq_matrix
**Download at**: lien vers arrayexpress (E-MTAB-14447) <br>
This folder contains the 3 different matrices generated after CellRanger mapping. Each matrix is divided in 3 files: the **matrix.mtx** is the count table; **barcode.tsv** contains the cell barcodes (rows of the matrices); and **features.tsv** contains the gene names (columns of the matrix).
### scRNA-seq_sequencing_output
**Download at**: lien vers arrayexpress (E-MTAB-14447) <br>
This folder contains the 3 sequencing results (in fastq format) from the scRNA-seq libraries. 
### seurat_objects
Stores the generated seurat objects.
### seurat_outputs
Stores the generated files comming from the seurat analysis (cell names, matrix, list of DEG, cells per clusters, cells per enhancers). They are usefull for later analysis with novoSpaRc. 