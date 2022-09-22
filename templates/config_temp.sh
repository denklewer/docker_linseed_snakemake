#!/bin/sh

CT=$1
PROJECT=$2

cat  << EOF
dataset: "${PROJECT}"
analysis_name: "${PROJECT}_ct${CT}"
num_inits: 5
top_mad: 100000
filter_genes: 1500
filter_samples: 0
scale_iterations: 20
cell_types: ${CT}
blocks_pipeline: "config/blocks.csv"
count:
  time: 6000 # minutes
  mem_ram: 32 # GB, max RAM
  threads: 8
  email: "aladyeva.e@wustl.edu"
  nodes: -1
  docker: "dockerreg01.accounts.ad.wustl.edu/artyomov_lab/docker_linseed_snakemake:cpp"
EOF
