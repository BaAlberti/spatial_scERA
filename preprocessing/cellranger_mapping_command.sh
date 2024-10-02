cellranger-7.2.0/cellranger count --id=mapping_6_1_scRNAseq_25E \
                   --transcriptome=../data/custom_genome/dm6.46_27_enhancers_30aug/index/ \
                   --fastqs=../data/Lib6_1 \
                   --sample=Lib6_1-SCI7T085-SCI5T085_22FH2GLT3 \
                   --r1-length=28 \
                   --localcores=20

cellranger-7.2.0/cellranger count --id=mapping_6_2_scRNAseq_25E \
                   --transcriptome=../data/custom_genome/dm6.46_27_enhancers_30aug/index/ \
                   --fastqs=../data/Lib6_2 \
                   --sample=Lib6_2-SCI7T002-SCI5T002_22FH2GLT3 \
                   --r1-length=28 \
                   --localcores=20

cellranger-7.2.0/cellranger count --id=mapping_6_3_scRNAseq_25E \
                   --transcriptome=../data/custom_genome/dm6.46_27_enhancers_30aug/index/ \
                   --fastqs=../data/reads/Lib6_3 \
                   --sample=Lib6_3-SCI7T014-SCI5T014_22FH2GLT3 \
                   --r1-length=28 \
                   --localcores=20