## Run command
```{bash}
/app/scripts/run_linseedv2.py --min_ct 8 --max_ct 8 \
--snakemake_path "/scratch1/fs1/martyomov/aladyevae/Deconvolution/snakemake/HNSC_TPM_F_LIMITX_MAD_2_GENES_1000_X_Sub/09_2022" \
--data_path "/storage1/fs1/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/data/datasets/HNSC_TPM_F/09_2022" \
--inits_path "/storage1/fs1/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/data/inits/HNSC_TPM_F_LIMITX_MAD_2_GENES_1000_X_Sub/09_2022" \
--results_path "/storage1/fs1/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/data/results/HNSC_TPM_F_LIMITX_MAD_2_GENES_1000_X_Sub/09_2022" \
--reports_path "/storage1/fs1/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/data/reports/HNSC_TPM_F_LIMITX_MAD_2_GENES_1000_X_Sub/09_2022" \
--blocks_path "/storage1/fs1/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/snakemake/templates/LIMIT_X/blocks.csv" \
--snakefile_path "/app/templates/Snakefile" --num_inits 10 \
--min_mad 2 --filter_genes 1000 --dataset "HNSC_TPM_F" --analysis_name "HNSC_TPM_F_LIMITX_MAD_2_GENES_1000_X_Sub" -b --init_strategy "SelectX"
```

## track snakemake logs
```{bash}
touch -f $P_LOG && tail -f $P_LOG
touch -f ./logs/pipeline.log && tail -f ./logs/pipeline.log
```
