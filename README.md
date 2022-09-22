#BEFORE START

```
PROJECT=HNSC_HK_CC_MALE
DATA_DIR=/storage1/fs1/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/data
PROCESS_DIR=/scratch1/fs1/martyomov/aladyevae/Deconvolution/snakemake
TEMPLATES_DIR=/storage1/fs1/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/snakemake/templates/
MIN_CT=6
MAX_CT=12
```

```
mkdir -p $DATA_DIR/datasets/$PROJECT/07_2022

for i in $(seq $MIN_CT $MAX_CT)
do 
mkdir -p $DATA_DIR/inits/$PROJECT/07_2022/ct$i
mkdir -p $DATA_DIR/results/$PROJECT/07_2022/ct$i
mkdir -p $DATA_DIR/config/$PROJECT/07_2022/ct$i

cp /storage1/fs1/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/snakemake/templates/blocks.csv $DATA_DIR/config/$PROJECT/07_2022/ct$i/blocks.csv
/storage1/fs1/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/snakemake/templates/config_temp.sh $i $PROJECT > $DATA_DIR/config/$PROJECT/07_2022/ct$i/config.yaml

done
```

1. Prepare dataset (need to be in RDS file)
2. Prepare configs - config.yaml and blocks.csv


```
for i in $(seq $MIN_CT $MAX_CT)
do 
mkdir -p $PROCESS_DIR/$PROJECT/07_2022/ct$i/resources
mkdir -p $PROCESS_DIR/$PROJECT/07_2022/ct$i/config

cp $TEMPLATES_DIR/Snakefile $PROCESS_DIR/$PROJECT/07_2022/ct$i/

ln -s $DATA_DIR/datasets/$PROJECT/07_2022 $PROCESS_DIR/$PROJECT/07_2022/ct$i/resources/datasets
cp $DATA_DIR/config/$PROJECT/07_2022/ct$i/* $PROCESS_DIR/$PROJECT/07_2022/ct$i/config/
ln -s $DATA_DIR/inits/$PROJECT/07_2022/ct$i $PROCESS_DIR/$PROJECT/07_2022/ct$i/resources/inits
ln -s $DATA_DIR/results/$PROJECT/07_2022/ct$i $PROCESS_DIR/$PROJECT/07_2022/ct$i/results

done
```

```
for i in $(seq $MIN_CT $MAX_CT)
do 
cd $PROCESS_DIR/$PROJECT/07_2022/ct$i

P_WD=`pwd`
mkdir -p "$P_WD/tmp"
echo "__LSF_JOB_CUSTOM_TMPDIR__=$P_WD/tmp" > lsf_docker_env_file.env
chmod a+r lsf_docker_env_file.env

LSF_DOCKER_ENV_FILE=$P_WD/lsf_docker_env_file.env
export SMK_DOCKER_IMG="dockerreg01.accounts.ad.wustl.edu/artyomov_lab/docker_linseed_snakemake:cpp"
export P_LOG=$P_WD/logs/pipeline.log

L_CORES=4; LSF_DOCKER_ENV_FILE=$P_WD/lsf_docker_env_file.env; mkdir -p logs; bsub -cwd $HOME -n $L_CORES -G compute-martyomov -q general -oo $P_LOG -R 'span[hosts=1]' -a "docker($SMK_DOCKER_IMG)" /usr/bin/script -fqe /dev/null  -c "source /etc/bash.bashrc; cd $P_WD; export TMPDIR=$P_WD/tmp; snakemake --profile lsf_demo --local-cores $L_CORES --jobs 50 -pr --conda-frontend conda --restart-times 0"

done
```

# track snakemake logs
touch -f $P_LOG && tail -f $P_LOG
touch -f ./logs/pipeline.log && tail -f ./logs/pipeline.log
