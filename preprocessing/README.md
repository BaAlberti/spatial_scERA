# Description of the two shell scripts 

## cellranger_mapping_command.sh
This script is making the mapping step for each of the 3 datasets. <br>
For reproducibility, you should work with the version 7.2.0 of cellranger. If you are not trying to replicate the results any version of CellRanger will work. <br>
The reference genome must be downloaded on zenodo () and stored in data/custom_genome ; The fastq files should be downloaded on ArrayExpress () and stored in data/scRNA-seq_sequencing_output. <br>
You can change the localcore parameters depending on the caracteristics of your computer/cluster

## cutadapt.sh
This script is trimming the PCR reads. Removed reads are the one which does not contain the cell+UMI and enhancer barcodes + the three common sequences of the reporter construct. <br>
The script works from the fastq reads downloaded on ArrayExpress () and stored in data/PCR_sequencing_output. <br>
The outputs will be stored in the data/cutadapt_output folder.