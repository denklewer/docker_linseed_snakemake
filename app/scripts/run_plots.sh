#!/bin/bash
input=$1
tool=$2
while IFS= read -r dataset
do
  echo ${dataset}
  RScript "/Users/aladyeva.e/Dropbox (ArtyomovLab)/Deconvolution/Docker/docker_linseed_snakemake/app/scripts/plot_data_to_SC_${tool}.R" \
--dataset "/Users/aladyeva.e/Downloads/HLCA_SCN" \
--results "/Volumes/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/data/results/${dataset}/08_2022" \
--save "/Volumes/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/data/reports/${dataset}/HLCA" \
--min_ct 3 --max_ct 10
done < "$input"