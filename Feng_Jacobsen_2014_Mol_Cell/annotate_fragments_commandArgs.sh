#!/bin/bash

i=$1

paste -d "\t" TAIR10_fragments_${i}.txt \
<(awk '{print $3}' /home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/HiC_comparison/log2_REC8_HA_Rep1_ChIP_genome_norm_coverage_${i}.txt) \
<(awk '{print $3}' /home/ajt200/analysis/180622_Chris_lambing_ChIP_REC8_HA_Col_kss/kss/coverage/common_input_MYC_Rep2/log2ChIPinput/genomeProfiles/HiC_comparison/log2_kss_REC8_HA_Rep1_ChIP_genome_norm_coverage_${i}.txt) \
> TAIR10_fragments_${i}_ann_tmp.txt
echo -e "fragmentNumber\tchrom\tstart\tend\tden_wt_REC8_HA_Rep1\tden_kss_REC8_HA_Rep1" \
| cat - <(tail -n +2 TAIR10_fragments_${i}_ann_tmp.txt) > TAIR10_fragments_${i}_ann.txt
rm TAIR10_fragments_${i}_ann_tmp.txt
