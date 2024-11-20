# Description of the two shell scripts 

## cellranger_mapping_command.sh
This script is making the mapping step for each of the 3 datasets. <br>
For reproducibility, you should work with the version 7.2.0 of cellranger. If you are not trying to replicate the results any version of CellRanger will work. <br>
The reference genome must be downloaded on [zenodo][zenodo] and stored in data/custom_genome ; The fastq files should be downloaded on [ArrayExpress][E-MTAB-14447] and stored in data/scRNA-seq_sequencing_output. <br>
You can change the localcore parameters depending on the caracteristics of your computer/cluster

## cutadapt.sh
This script is trimming the PCR reads. Removed reads are the one which does not contain the cell+UMI and enhancer barcodes + the three common sequences of the reporter construct. <br>
The script works from the fastq reads downloaded on [ArrayExpress][E-MTAB-14445] and stored in data/PCR_sequencing_output. <br>
The outputs will be stored in the data/cutadapt_output folder.

[zenodo]: https://zenodo.org/records/14006160
[E-MTAB-14447]: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14447?key=7ffc7424-efbd-4a63-ab75-07c486994e0b
[E-MTAB-14445]: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14445?key=d55c8c87-d977-4157-979a-cf7a160d96de