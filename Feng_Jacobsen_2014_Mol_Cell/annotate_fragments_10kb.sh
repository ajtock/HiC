#!/bin/bash

paste -d "\t" TAIR10_fragments_10kb.txt \
<(awk '{print $3}' /home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/HiC_comparison/log2_REC8_HA_Rep1_ChIP_genome_norm_coverage_10kb.txt) \
<(awk '{print $3}' /home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/WT/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/HiC_comparison/log2_REC8_HA_Rep2_ChIP_genome_norm_coverage_10kb.txt) \
<(awk '{print $3}' /home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/HiC_comparison/log2_REC8_MYC_Rep1_ChIP_genome_norm_coverage_10kb.txt) \
<(awk '{print $3}' /home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/kss/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/HiC_comparison/log2_kss_REC8_HA_Rep1_ChIP_genome_norm_coverage_10kb.txt) \
<(awk '{print $3}' /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/genomeProfiles/HiC_comparison/log2_WT_SPO11oligo_RPI1_genome_norm_coverage_10kb.txt) \
<(awk '{print $3}' /projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/log2ChIPinput/genomeProfiles/HiC_comparison/log2_suvh456_SPO11oligo_RPI34_genome_norm_coverage_10kb.txt) \
<(awk '{print $3}' /projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/log2ChIPinput/genomeProfiles/HiC_comparison/log2_suvh456_SPO11oligo_RPI35_genome_norm_coverage_10kb.txt) \
<(awk '{print $3}' /projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/genomeProfiles/HiC_comparison/log2_WT_H3K9me2_ChIP_genome_norm_coverage_10kb.txt) \
<(awk '{print $3}' /projects/ajt200/BAM_masters/H3K9me2/kss/coverage/log2ChIPinput/genomeProfiles/HiC_comparison/log2_kss_H3K9me2_ChIP_genome_norm_coverage_10kb.txt) \
> TAIR10_fragments_10kb_ann_tmp.txt
echo -e "fragmentNumber\tchrom\tstart\tend\tden_wt_REC8_HA_Rep1\tden_wt_REC8_HA_Rep2\tden_wt_REC8_MYC_Rep1\tden_kss_REC8_HA_Rep1\tden_wt_SPO11oligos_Rep1\tden_kss_SPO11oligos_Rep1\tden_kss_SPO11oligos_Rep2\tden_wt_H3K9me2\tden_kss_H3K9me2" \
| cat - <(tail -n +2 TAIR10_fragments_10kb_ann_tmp.txt) > TAIR10_fragments_10kb_ann.txt
rm TAIR10_fragments_10kb_ann_tmp.txt
