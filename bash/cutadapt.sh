cutadapt -g 'GACGTCatcgtcctgcagg' -G 'CTACACGACGCTCTTCCGATCT' \
--untrimmed-output trimming_libS6.3/CO_untrimmed_5p.fastq \
--untrimmed-paired-output trimming_libS6.3/CO_untrimmed_3p.fastq \
--pair-filter=any \
-O 10 \
-o trimming_libS6.3/out1.bigR1.fastq -p trimming_libS6.3/out1.bigR2.fastq \
LibS6.3/LibS6_EKDL230017004-1A_22FH2GLT3_L2_1.fq LibS6.3/LibS6_EKDL230017004-1A_22FH2GLT3_L2_2.fq > trimming_libS6.3/stat_step1.txt

cutadapt -g 'CTACACGACGCTCTTCCGATCT' -G 'GACGTCatcgtcctgcagg' \
--untrimmed-output trimming_libS6.3/RC_untrimmed_5p.fastq \
--untrimmed-paired-output trimming_libS6.3/RC_untrimmed_3p.fastq \
--minimum-length 100 \
--too-short-output trimming_libS6.3/step2_too_short_5p.fastq \
--too-short-paired-output trimming_libS6.3/step2_too_short_3p.fastq \
--pair-filter=any \
-O 10 \
-o trimming_libS6.3/out2.bigR2_RC.fastq -p trimming_libS6.3/out2.bigR1_RC.fastq \
trimming_libS6.3/CO_untrimmed_5p.fastq trimming_libS6.3/CO_untrimmed_3p.fastq > trimming_libS6.3/stat_step2.txt

seqkit grep -s -P -m 2 -j 10 -R 29:60 -p TTTTTTTTTTTTTTTTTTTT trimming_libS6.3/RC_untrimmed_5p.fastq > trimming_libS6.3/polyT_certain_reads_missoriented_step2.fastq # Je récupère les reads qui n'ont pas marché au step 2 car ils n'ont plus le R1 (pas infaillible mais pas mieux)
seqkit seq trimming_libS6.3/polyT_certain_reads_missoriented_step2.fastq -n -i > trimming_libS6.3/list_ID.txt # extraits les ID des reads qui n'ont pas le R1
seqkit grep -f trimming_libS6.3/list_ID.txt LibS6.3/LibS6_EKDL230017004-1A_22FH2GLT3_L2_2.fq -o trimming_libS6.3/R2_RC_wo_R1.fastq


cutadapt -g 'GACGTCatcgtcctgcagg' \
--untrimmed-output trimming_libS6.3/CO_untrimmed_5p_step3.fastq \
--untrimmed-paired-output trimming_libS6.3/CO_untrimmed_3p_step3.fastq \
--pair-filter=any \
-O 10 \
-o trimming_libS6.3/out3.bigR1.fastq -p trimming_libS6.3/out3.bigR2.fastq \
trimming_libS6.3/R2_RC_wo_R1.fastq trimming_libS6.3/polyT_certain_reads_missoriented_step2.fastq > trimming_libS6.3/stat_step3.txt

cat trimming_libS6.3/out2.bigR1_RC.fastq trimming_libS6.3/out3.bigR1.fastq >> trimming_libS6.3/out1.bigR1.fastq
cat trimming_libS6.3/out2.bigR2_RC.fastq trimming_libS6.3/out3.bigR2.fastq >> trimming_libS6.3/out1.bigR2.fastq

cutadapt -a 'GCTGCCGCTTCGAGCAGACATGCATATG' \
--discard-untrimmed \
-O 15 \
-o trimming_libS6.3/out4.bigR1.fastq -p trimming_libS6.3/out4.bigR2.fastq \
trimming_libS6.3/out1.bigR1.fastq trimming_libS6.3/out1.bigR2.fastq > trimming_libS6.3/stat_step4.txt

cutadapt -l 28 -o trimming_libS6.3/C_BC.fastq trimming_libS6.3/out4.bigR2.fastq
mv trimming_libS6.3/out4.bigR1.fastq trimming_libS6.3/E_BC.fastq

###

cutadapt -g 'GACGTCatcgtcctgcagg' -G 'CTACACGACGCTCTTCCGATCT' \
--untrimmed-output trimming_libS6.1/CO_untrimmed_5p.fastq \
--untrimmed-paired-output trimming_libS6.1/CO_untrimmed_3p.fastq \
--pair-filter=any \
-O 10 \
-o trimming_libS6.1/out1.bigR1.fastq -p trimming_libS6.1/out1.bigR2.fastq \
LibS6.1/LibS61_MKDL240001181-1A_22HWWWLT3_L7_1.fq.gz LibS6.1/LibS61_MKDL240001181-1A_22HWWWLT3_L7_2.fq.gz > trimming_libS6.1/stat_step1.txt

cutadapt -g 'CTACACGACGCTCTTCCGATCT' -G 'GACGTCatcgtcctgcagg' \
--untrimmed-output trimming_libS6.1/RC_untrimmed_5p.fastq \
--untrimmed-paired-output trimming_libS6.1/RC_untrimmed_3p.fastq \
--pair-filter=any \
-O 10 \
-o trimming_libS6.1/out2.bigR2_RC.fastq -p trimming_libS6.1/out2.bigR1_RC.fastq \
trimming_libS6.1/CO_untrimmed_5p.fastq trimming_libS6.1/CO_untrimmed_3p.fastq > trimming_libS6.1/stat_step2.txt

seqkit grep -s -P -m 2 -j 10 -R 29:60 -p TTTTTTTTTTTTTTTTTTTT trimming_libS6.1/RC_untrimmed_5p.fastq > trimming_libS6.1/polyT_certain_reads_missoriented_step2.fastq # Je récupère les reads qui n'ont pas marché au step 2 car ils n'ont plus le R1 (pas infaillible mais pas mieux)
seqkit seq trimming_libS6.1/polyT_certain_reads_missoriented_step2.fastq -n -i > trimming_libS6.1/list_ID.txt # extraits les ID des reads qui n'ont pas le R1
seqkit grep -f trimming_libS6.1/list_ID.txt LibS6.1/LibS61_MKDL240001181-1A_22HWWWLT3_L7_2.fq.gz -o trimming_libS6.1/R2_RC_wo_R1.fastq


cutadapt -g 'GACGTCatcgtcctgcagg' \
--untrimmed-output trimming_libS6.1/CO_untrimmed_5p_step3.fastq \
--untrimmed-paired-output trimming_libS6.1/CO_untrimmed_3p_step3.fastq \
--pair-filter=any \
-O 10 \
-o trimming_libS6.1/out3.bigR1.fastq -p trimming_libS6.1/out3.bigR2.fastq \
trimming_libS6.1/R2_RC_wo_R1.fastq trimming_libS6.1/polyT_certain_reads_missoriented_step2.fastq > trimming_libS6.1/stat_step3.txt

cat trimming_libS6.1/out2.bigR1_RC.fastq trimming_libS6.1/out3.bigR1.fastq >> trimming_libS6.1/out1.bigR1.fastq
cat trimming_libS6.1/out2.bigR2_RC.fastq trimming_libS6.1/out3.bigR2.fastq >> trimming_libS6.1/out1.bigR2.fastq

cutadapt -a 'GCTGCCGCTTCGAGCAGACATGCATATG' \
--discard-untrimmed \
-O 15 \
-o trimming_libS6.1/out4.bigR1.fastq -p trimming_libS6.1/out4.bigR2.fastq \
trimming_libS6.1/out1.bigR1.fastq trimming_libS6.1/out1.bigR2.fastq > trimming_libS6.1/stat_step4.txt

cutadapt -l 28 -o trimming_libS6.1/C_BC.fastq trimming_libS6.1/out4.bigR2.fastq
mv trimming_libS6.1/out4.bigR1.fastq trimming_libS6.1/E_BC.fastq

#####

cutadapt -g 'GACGTCatcgtcctgcagg' -G 'CTACACGACGCTCTTCCGATCT' \
--untrimmed-output trimming_libS6.2/CO_untrimmed_5p.fastq \
--untrimmed-paired-output trimming_libS6.2/CO_untrimmed_3p.fastq \
--pair-filter=any \
-O 10 \
-o trimming_libS6.2/out1.bigR1.fastq -p trimming_libS6.2/out1.bigR2.fastq \
LibS6.2/LibS62_MKDL240002948-1A_22MG3MLT3_L7_1.fq.gz LibS6.2/LibS62_MKDL240002948-1A_22MG3MLT3_L7_2.fq.gz > trimming_libS6.2/stat_step1.txt

cutadapt -g 'CTACACGACGCTCTTCCGATCT' -G 'GACGTCatcgtcctgcagg' \
--untrimmed-output trimming_libS6.2/RC_untrimmed_5p.fastq \
--untrimmed-paired-output trimming_libS6.2/RC_untrimmed_3p.fastq \
--minimum-length 100 \
--too-short-output trimming_libS6.2/step2_too_short_5p.fastq \
--too-short-paired-output trimming_libS6.2/step2_too_short_3p.fastq \
--pair-filter=any \
-O 10 \
-o trimming_libS6.2/out2.bigR2_RC.fastq -p trimming_libS6.2/out2.bigR1_RC.fastq \
trimming_libS6.2/CO_untrimmed_5p.fastq trimming_libS6.2/CO_untrimmed_3p.fastq > trimming_libS6.2/stat_step2.txt

seqkit grep -s -P -m 2 -j 10 -R 29:60 -p TTTTTTTTTTTTTTTTTTTT trimming_libS6.2/RC_untrimmed_5p.fastq > trimming_libS6.2/polyT_certain_reads_missoriented_step2.fastq # Je récupère les reads qui n'ont pas marché au step 2 car ils n'ont plus le R1 (pas infaillible mais pas mieux)
seqkit seq trimming_libS6.2/polyT_certain_reads_missoriented_step2.fastq -n -i > trimming_libS6.2/list_ID.txt # extraits les ID des reads qui n'ont pas le R1
seqkit grep -f trimming_libS6.2/list_ID.txt LibS6.2/LibS62_MKDL240002948-1A_22MG3MLT3_L7_2.fq.gz -o trimming_libS6.2/R2_RC_wo_R1.fastq


cutadapt -g 'GACGTCatcgtcctgcagg' \
--untrimmed-output trimming_libS6.2/CO_untrimmed_5p_step3.fastq \
--untrimmed-paired-output trimming_libS6.2/CO_untrimmed_3p_step3.fastq \
--pair-filter=any \
-O 10 \
-o trimming_libS6.2/out3.bigR1.fastq -p trimming_libS6.2/out3.bigR2.fastq \
trimming_libS6.2/R2_RC_wo_R1.fastq trimming_libS6.2/polyT_certain_reads_missoriented_step2.fastq > trimming_libS6.2/stat_step3.txt

cat trimming_libS6.2/out2.bigR1_RC.fastq trimming_libS6.2/out3.bigR1.fastq >> trimming_libS6.2/out1.bigR1.fastq
cat trimming_libS6.2/out2.bigR2_RC.fastq trimming_libS6.2/out3.bigR2.fastq >> trimming_libS6.2/out1.bigR2.fastq

cutadapt -a 'GCTGCCGCTTCGAGCAGACATGCATATG' \
--discard-untrimmed \
-O 15 \
-o trimming_libS6.2/out4.bigR1.fastq -p trimming_libS6.2/out4.bigR2.fastq \
trimming_libS6.2/out1.bigR1.fastq trimming_libS6.2/out1.bigR2.fastq > trimming_libS6.2/stat_step4.txt

cutadapt -l 28 -o trimming_libS6.2/C_BC.fastq trimming_libS6.2/out4.bigR2.fastq
mv trimming_libS6.2/out4.bigR1.fastq trimming_libS6.2/E_BC.fastq