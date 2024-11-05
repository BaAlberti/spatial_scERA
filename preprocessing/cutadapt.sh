cutadapt -g 'GACGTCatcgtcctgcagg' -G 'CTACACGACGCTCTTCCGATCT' \
--untrimmed-output ../data/cutadapt_output/trimming_libS6.1/CO_untrimmed_5p.fastq \
--untrimmed-paired-output ../data/cutadapt_output/trimming_libS6.1/CO_untrimmed_3p.fastq \
--pair-filter=any \
-O 10 \
-o ../data/cutadapt_output/trimming_libS6.1/out1.bigR1.fastq -p ../data/cutadapt_output/trimming_libS6.1/out1.bigR2.fastq \
../data/PCR_sequencing_output/LibS6.1/LibS61_MKDL240001181-1A_22HWWWLT3_L7_1.fq.gz ../data/PCR_sequencing_output/LibS6.1/LibS61_MKDL240001181-1A_22HWWWLT3_L7_2.fq.gz > ../data/cutadapt_output/trimming_libS6.1/stat_step1.txt

cutadapt -g 'CTACACGACGCTCTTCCGATCT' -G 'GACGTCatcgtcctgcagg' \
--untrimmed-output ../data/cutadapt_output/trimming_libS6.1/RC_untrimmed_5p.fastq \
--untrimmed-paired-output ../data/cutadapt_output/trimming_libS6.1/RC_untrimmed_3p.fastq \
--pair-filter=any \
-O 10 \
-o ../data/cutadapt_output/trimming_libS6.1/out2.bigR2_RC.fastq -p ../data/cutadapt_output/trimming_libS6.1/out2.bigR1_RC.fastq \
../data/cutadapt_output/trimming_libS6.1/CO_untrimmed_5p.fastq ../data/cutadapt_output/trimming_libS6.1/CO_untrimmed_3p.fastq > ../data/cutadapt_output/trimming_libS6.1/stat_step2.txt

seqkit grep -s -P -m 2 -j 10 -R 29:60 -p TTTTTTTTTTTTTTTTTTTT ../data/cutadapt_output/trimming_libS6.1/RC_untrimmed_5p.fastq > ../data/cutadapt_output/trimming_libS6.1/polyT_certain_reads_missoriented_step2.fastq # Je récupère les reads qui n'ont pas marché au step 2 car ils n'ont plus le R1 (pas infaillible mais pas mieux)
seqkit seq ../data/cutadapt_output/trimming_libS6.1/polyT_certain_reads_missoriented_step2.fastq -n -i > ../data/cutadapt_output/trimming_libS6.1/list_ID.txt # extraits les ID des reads qui n'ont pas le R1
seqkit grep -f ../data/cutadapt_output/trimming_libS6.1/list_ID.txt ../data/PCR_sequencing_output/LibS6.1/LibS61_MKDL240001181-1A_22HWWWLT3_L7_2.fq.gz -o ../data/cutadapt_output/trimming_libS6.1/R2_RC_wo_R1.fastq


cutadapt -g 'GACGTCatcgtcctgcagg' \
--untrimmed-output ../data/cutadapt_output/trimming_libS6.1/CO_untrimmed_5p_step3.fastq \
--untrimmed-paired-output ../data/cutadapt_output/trimming_libS6.1/CO_untrimmed_3p_step3.fastq \
--pair-filter=any \
-O 10 \
-o ../data/cutadapt_output/trimming_libS6.1/out3.bigR1.fastq -p ../data/cutadapt_output/trimming_libS6.1/out3.bigR2.fastq \
../data/cutadapt_output/trimming_libS6.1/R2_RC_wo_R1.fastq ../data/cutadapt_output/trimming_libS6.1/polyT_certain_reads_missoriented_step2.fastq > ../data/cutadapt_output/trimming_libS6.1/stat_step3.txt

cat ../data/cutadapt_output/trimming_libS6.1/out2.bigR1_RC.fastq ../data/cutadapt_output/trimming_libS6.1/out3.bigR1.fastq >> ../data/cutadapt_output/trimming_libS6.1/out1.bigR1.fastq
cat ../data/cutadapt_output/trimming_libS6.1/out2.bigR2_RC.fastq ../data/cutadapt_output/trimming_libS6.1/out3.bigR2.fastq >> ../data/cutadapt_output/trimming_libS6.1/out1.bigR2.fastq

cutadapt -a 'GCTGCCGCTTCGAGCAGACATGCATATG' \
--discard-untrimmed \
-O 15 \
-o ../data/cutadapt_output/trimming_libS6.1/out4.bigR1.fastq -p ../data/cutadapt_output/trimming_libS6.1/out4.bigR2.fastq \
../data/cutadapt_output/trimming_libS6.1/out1.bigR1.fastq ../data/cutadapt_output/trimming_libS6.1/out1.bigR2.fastq > ../data/cutadapt_output/trimming_libS6.1/stat_step4.txt

cutadapt -l 28 -o ../data/cutadapt_output/trimming_libS6.1/C_BC.fastq ../data/cutadapt_output/trimming_libS6.1/out4.bigR2.fastq
mv ../data/cutadapt_output/trimming_libS6.1/out4.bigR1.fastq ../data/cutadapt_output/trimming_libS6.1/E_BC.fastq

###

cutadapt -g 'GACGTCatcgtcctgcagg' -G 'CTACACGACGCTCTTCCGATCT' \
--untrimmed-output ../data/cutadapt_output/trimming_libS6.2/CO_untrimmed_5p.fastq \
--untrimmed-paired-output ../data/cutadapt_output/trimming_libS6.2/CO_untrimmed_3p.fastq \
--pair-filter=any \
-O 10 \
-o ../data/cutadapt_output/trimming_libS6.2/out1.bigR1.fastq -p ../data/cutadapt_output/trimming_libS6.2/out1.bigR2.fastq \
../data/PCR_sequencing_output/LibS6.2/LibS62_MKDL240002948-1A_22MG3MLT3_L7_1.fq.gz ../data/PCR_sequencing_output/LibS6.2/LibS62_MKDL240002948-1A_22MG3MLT3_L7_2.fq.gz > ../data/cutadapt_output/trimming_libS6.2/stat_step1.txt

cutadapt -g 'CTACACGACGCTCTTCCGATCT' -G 'GACGTCatcgtcctgcagg' \
--untrimmed-output ../data/cutadapt_output/trimming_libS6.2/RC_untrimmed_5p.fastq \
--untrimmed-paired-output ../data/cutadapt_output/trimming_libS6.2/RC_untrimmed_3p.fastq \
--minimum-length 100 \
--too-short-output ../data/cutadapt_output/trimming_libS6.2/step2_too_short_5p.fastq \
--too-short-paired-output ../data/cutadapt_output/trimming_libS6.2/step2_too_short_3p.fastq \
--pair-filter=any \
-O 10 \
-o ../data/cutadapt_output/trimming_libS6.2/out2.bigR2_RC.fastq -p ../data/cutadapt_output/trimming_libS6.2/out2.bigR1_RC.fastq \
../data/cutadapt_output/trimming_libS6.2/CO_untrimmed_5p.fastq ../data/cutadapt_output/trimming_libS6.2/CO_untrimmed_3p.fastq > ../data/cutadapt_output/trimming_libS6.2/stat_step2.txt

seqkit grep -s -P -m 2 -j 10 -R 29:60 -p TTTTTTTTTTTTTTTTTTTT ../data/cutadapt_output/trimming_libS6.2/RC_untrimmed_5p.fastq > ../data/cutadapt_output/trimming_libS6.2/polyT_certain_reads_missoriented_step2.fastq # Je récupère les reads qui n'ont pas marché au step 2 car ils n'ont plus le R1 (pas infaillible mais pas mieux)
seqkit seq ../data/cutadapt_output/trimming_libS6.2/polyT_certain_reads_missoriented_step2.fastq -n -i > ../data/cutadapt_output/trimming_libS6.2/list_ID.txt # extraits les ID des reads qui n'ont pas le R1
seqkit grep -f ../data/cutadapt_output/trimming_libS6.2/list_ID.txt ../data/PCR_sequencing_output/LibS6.2/LibS62_MKDL240002948-1A_22MG3MLT3_L7_2.fq.gz -o ../data/cutadapt_output/trimming_libS6.2/R2_RC_wo_R1.fastq


cutadapt -g 'GACGTCatcgtcctgcagg' \
--untrimmed-output ../data/cutadapt_output/trimming_libS6.2/CO_untrimmed_5p_step3.fastq \
--untrimmed-paired-output ../data/cutadapt_output/trimming_libS6.2/CO_untrimmed_3p_step3.fastq \
--pair-filter=any \
-O 10 \
-o ../data/cutadapt_output/trimming_libS6.2/out3.bigR1.fastq -p ../data/cutadapt_output/trimming_libS6.2/out3.bigR2.fastq \
../data/cutadapt_output/trimming_libS6.2/R2_RC_wo_R1.fastq ../data/cutadapt_output/trimming_libS6.2/polyT_certain_reads_missoriented_step2.fastq > ../data/cutadapt_output/trimming_libS6.2/stat_step3.txt

cat ../data/cutadapt_output/trimming_libS6.2/out2.bigR1_RC.fastq ../data/cutadapt_output/trimming_libS6.2/out3.bigR1.fastq >> ../data/cutadapt_output/trimming_libS6.2/out1.bigR1.fastq
cat ../data/cutadapt_output/trimming_libS6.2/out2.bigR2_RC.fastq ../data/cutadapt_output/trimming_libS6.2/out3.bigR2.fastq >> ../data/cutadapt_output/trimming_libS6.2/out1.bigR2.fastq

cutadapt -a 'GCTGCCGCTTCGAGCAGACATGCATATG' \
--discard-untrimmed \
-O 15 \
-o ../data/cutadapt_output/trimming_libS6.2/out4.bigR1.fastq -p ../data/cutadapt_output/trimming_libS6.2/out4.bigR2.fastq \
../data/cutadapt_output/trimming_libS6.2/out1.bigR1.fastq ../data/cutadapt_output/trimming_libS6.2/out1.bigR2.fastq > ../data/cutadapt_output/trimming_libS6.2/stat_step4.txt

cutadapt -l 28 -o ../data/cutadapt_output/trimming_libS6.2/C_BC.fastq ../data/cutadapt_output/trimming_libS6.2/out4.bigR2.fastq
mv ../data/cutadapt_output/trimming_libS6.2/out4.bigR1.fastq ../data/cutadapt_output/trimming_libS6.2/E_BC.fastq

###

cutadapt -g 'GACGTCatcgtcctgcagg' -G 'CTACACGACGCTCTTCCGATCT' \
--untrimmed-output ../data/cutadapt_output/trimming_libS6.3/CO_untrimmed_5p.fastq \
--untrimmed-paired-output ../data/cutadapt_output/trimming_libS6.3/CO_untrimmed_3p.fastq \
--pair-filter=any \
-O 10 \
-o ../data/cutadapt_output/trimming_libS6.3/out1.bigR1.fastq -p ../data/cutadapt_output/trimming_libS6.3/out1.bigR2.fastq \
../data/PCR_sequencing_output/LibS6.3/LibS6_EKDL230017004-1A_22FH2GLT3_L2_1.fq ../data/PCR_sequencing_output/LibS6.3/LibS6_EKDL230017004-1A_22FH2GLT3_L2_2.fq > ../data/cutadapt_output/trimming_libS6.3/stat_step1.txt

cutadapt -g 'CTACACGACGCTCTTCCGATCT' -G 'GACGTCatcgtcctgcagg' \
--untrimmed-output ../data/cutadapt_output/trimming_libS6.3/RC_untrimmed_5p.fastq \
--untrimmed-paired-output ../data/cutadapt_output/trimming_libS6.3/RC_untrimmed_3p.fastq \
--minimum-length 100 \
--too-short-output ../data/cutadapt_output/trimming_libS6.3/step2_too_short_5p.fastq \
--too-short-paired-output ../data/cutadapt_output/trimming_libS6.3/step2_too_short_3p.fastq \
--pair-filter=any \
-O 10 \
-o ../data/cutadapt_output/trimming_libS6.3/out2.bigR2_RC.fastq -p ../data/cutadapt_output/trimming_libS6.3/out2.bigR1_RC.fastq \
../data/cutadapt_output/trimming_libS6.3/CO_untrimmed_5p.fastq ../data/cutadapt_output/trimming_libS6.3/CO_untrimmed_3p.fastq > ../data/cutadapt_output/trimming_libS6.3/stat_step2.txt

seqkit grep -s -P -m 2 -j 10 -R 29:60 -p TTTTTTTTTTTTTTTTTTTT ../data/cutadapt_output/trimming_libS6.3/RC_untrimmed_5p.fastq > ../data/cutadapt_output/trimming_libS6.3/polyT_certain_reads_missoriented_step2.fastq # Je récupère les reads qui n'ont pas marché au step 2 car ils n'ont plus le R1 (pas infaillible mais pas mieux)
seqkit seq ../data/cutadapt_output/trimming_libS6.3/polyT_certain_reads_missoriented_step2.fastq -n -i > ../data/cutadapt_output/trimming_libS6.3/list_ID.txt # extraits les ID des reads qui n'ont pas le R1
seqkit grep -f ../data/cutadapt_output/trimming_libS6.3/list_ID.txt ../data/PCR_sequencing_output/LibS6.3/LibS6_EKDL230017004-1A_22FH2GLT3_L2_2.fq -o ../data/cutadapt_output/trimming_libS6.3/R2_RC_wo_R1.fastq


cutadapt -g 'GACGTCatcgtcctgcagg' \
--untrimmed-output ../data/cutadapt_output/trimming_libS6.3/CO_untrimmed_5p_step3.fastq \
--untrimmed-paired-output ../data/cutadapt_output/trimming_libS6.3/CO_untrimmed_3p_step3.fastq \
--pair-filter=any \
-O 10 \
-o ../data/cutadapt_output/trimming_libS6.3/out3.bigR1.fastq -p ../data/cutadapt_output/trimming_libS6.3/out3.bigR2.fastq \
../data/cutadapt_output/trimming_libS6.3/R2_RC_wo_R1.fastq ../data/cutadapt_output/trimming_libS6.3/polyT_certain_reads_missoriented_step2.fastq > ../data/cutadapt_output/trimming_libS6.3/stat_step3.txt

cat ../data/cutadapt_output/trimming_libS6.3/out2.bigR1_RC.fastq ../data/cutadapt_output/trimming_libS6.3/out3.bigR1.fastq >> ../data/cutadapt_output/trimming_libS6.3/out1.bigR1.fastq
cat ../data/cutadapt_output/trimming_libS6.3/out2.bigR2_RC.fastq ../data/cutadapt_output/trimming_libS6.3/out3.bigR2.fastq >> ../data/cutadapt_output/trimming_libS6.3/out1.bigR2.fastq

cutadapt -a 'GCTGCCGCTTCGAGCAGACATGCATATG' \
--discard-untrimmed \
-O 15 \
-o ../data/cutadapt_output/trimming_libS6.3/out4.bigR1.fastq -p ../data/cutadapt_output/trimming_libS6.3/out4.bigR2.fastq \
../data/cutadapt_output/trimming_libS6.3/out1.bigR1.fastq ../data/cutadapt_output/trimming_libS6.3/out1.bigR2.fastq > ../data/cutadapt_output/trimming_libS6.3/stat_step4.txt

cutadapt -l 28 -o ../data/cutadapt_output/trimming_libS6.3/C_BC.fastq ../data/cutadapt_output/trimming_libS6.3/out4.bigR2.fastq
mv ../data/cutadapt_output/trimming_libS6.3/out4.bigR1.fastq ../data/cutadapt_output/trimming_libS6.3/E_BC.fastq