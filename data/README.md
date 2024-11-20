# Description of the different folders and datatypes that you can find or load in this data folder

This folder is gathering every data generated along the analysis, or is storing data to run them. <br>
For several folders you need to download files, the links are provided in this README. As all the scripts are getting their inputs/outputs from here, it is encouraged to follow the structure of this data folder. <br>

### custom_genome
**Download at**: [zenodo][zenodo]. <br>
This folder contains the reference genome to use with CellRanger if you want to map the reads from the fastq files you can download at E-MTAB-14447
### cutadapt_output
Stores the cutadapt processed fastqs files from the preprocessing script **cutadapt.sh**. This folder contains the trimmed fastqs from PCR sequencing. They are used later in the **clean_trimming_script_novopaper.ipynb** script.
### novosparc_atlas
Stores the atlas used with novoSpaRc for this project.
### novosparc_dataset
**Download at**: [novosparc dataset][novosparc] <br>
This folder contains the dataset used in the tutorial of novosparc. It is used only for a figure comparing the stripe level in the reconstruction between our and this tutorial datasets.
### PCR_sequencing_output
**Download at**:  [PCR data][E-MTAB-14445] <br>
This folder contains the three sequenced library of the PCR targeted amplification.
### scRNA-seq_matrix
**Download at**: [scRNA-seq data][E-MTAB-14447] <br>
This folder contains the 3 different matrices generated after CellRanger mapping. Each matrix is divided in 3 files: the **matrix.mtx** is the count table; **barcode.tsv** contains the cell barcodes (rows of the matrices); and **features.tsv** contains the gene names (columns of the matrix).
### scRNA-seq_sequencing_output
**Download at**: [scRNA-seq data][E-MTAB-14447] <br>
This folder contains the 3 sequencing results (in fastq format) from the scRNA-seq libraries. 
### seurat_objects
Stores the generated seurat objects.
### seurat_outputs
Stores the generated files comming from the seurat analysis (cell names, matrix, list of DEG, cells per clusters, cells per enhancers). They are usefull for later analysis with novoSpaRc. 

[zenodo]: https://zenodo.org/records/14006160
[novosparc]: https://github.com/rajewsky-lab/novosparc/tree/master/novosparc/datasets/drosophila_scRNAseq
[E-MTAB-14445]: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14445?key=d55c8c87-d977-4157-979a-cf7a160d96de
[E-MTAB-14447]: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14447?key=7ffc7424-efbd-4a63-ab75-07c486994e0b