#!/opt/conda/bin/python

""" 
/app/scripts/run_linseedv2.py --min_ct 6 --max_ct 10 \
--snakemake_path "/home/aladyevae/projects/snakemake/deconvolution/HNSC_HK_CC_MALE_MAD10_GENES_2000/08_2022" \
--data_path "/home/aladyevae/projects/snakemake/deconvolution/data/datasets/HNSC_HK_CC_MALE/07_2022" \
--inits_path "/home/aladyevae/projects/snakemake/deconvolution/data/inits/HNSC_HK_CC_MALE_MAD10_GENES_2000/08_2022" \
--results_path "/home/aladyevae/projects/snakemake/deconvolution/data/results/HNSC_HK_CC_MALE_MAD10_GENES_2000/08_2022" \
--reports_path "/home/aladyevae/projects/snakemake/deconvolution/data/reports/HNSC_HK_CC_MALE_MAD10_GENES_2000/08_2022" \
--blocks_path "/home/aladyevae/projects/snakemake/deconvolution/templates/blocks.csv" \
--snakefile_path "/app/templates/Snakefile" --num_inits 10 \
--min_mad 10 --filter_genes 2000 --dataset "HNSC_HK_CC_MALE" --analysis_name "HNSC_HK_CC_MALE_MAD10_GENES_2000" --dt "20220802_220422" -l

/app/scripts/run_linseedv2.py --min_ct 6 --max_ct 10 \
--snakemake_path "/scratch1/fs1/martyomov/aladyevae/Deconvolution/snakemake/HNSC_HK_CC_MALE_MAD10_GENES_2000/08_2022" \
--data_path "/storage1/fs1/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/data/datasets/HNSC_HK_CC_MALE/07_2022" \
--inits_path "/storage1/fs1/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/data/inits/HNSC_HK_CC_MALE_MAD10_GENES_2000/08_2022" \
--results_path "/storage1/fs1/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/data/results/HNSC_HK_CC_MALE_MAD10_GENES_2000/08_2022" \
--reports_path "/storage1/fs1/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/data/reports/HNSC_HK_CC_MALE_MAD10_GENES_2000/08_2022" \
--blocks_path "/storage1/fs1/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/snakemake/templates/MAD10_GENES_2000/blocks.csv" \
--snakefile_path "/app/templates/Snakefile" --num_inits 10 \
--min_mad 10 --filter_genes 2000 --dataset "HNSC_HK_CC_MALE" --analysis_name "HNSC_HK_CC_MALE_MAD10_GENES_2000" -b --dt "20220809_161254"
"""

import os
from optparse import OptionParser
from subprocess import Popen
import datetime as dt
from pathlib import Path
import yaml
import shutil

parser = OptionParser()
parser.add_option("-d", "--dataset")
parser.add_option("-a", "--analysis_name")
parser.add_option("--data_path")
parser.add_option("--snakemake_path")
parser.add_option("--inits_path")
parser.add_option("--results_path")
parser.add_option("--reports_path")
parser.add_option("--blocks_path")
parser.add_option("--snakefile_path", default = "/app/templates/Snakefile")
parser.add_option("--dt")
parser.add_option("--init_strategy")
parser.add_option("--num_inits", type=int, default=5)
parser.add_option("--top_mad", type=int)
parser.add_option("--top_median", type=int)
parser.add_option("--min_mad", type=float)
parser.add_option("--max_mad", type=float)
parser.add_option("--min_median", type=float)
parser.add_option("--max_median", type=float)
parser.add_option("--thresh", type=int)
parser.add_option("--r_tilda", type=float)
parser.add_option("--scale_iterations", type=int, default=20)
parser.add_option("--no_cutoff_samples", action="store_true", dest="no_cutoff_samples")
parser.add_option("--no_cutoff_genes", action="store_true", dest="no_cutoff_genes")
parser.add_option("--no_filter_by_plane", action="store_true", dest="no_filter_by_plane")
parser.add_option("--min_ct", type=int)
parser.add_option("--max_ct", type=int)
parser.add_option("--preprocessing_cell_types", type=int, default=20)
parser.add_option("-l", action="store_true", dest="local", default=True)
parser.add_option("-b", action="store_false", dest="local")
parser.add_option("--apply_filters", action="store_true", dest="apply_filters", default=False)
parser.add_option("--docker_image", default = "dockerreg01.accounts.ad.wustl.edu/artyomov_lab/docker_linseed_snakemake")
parser.add_option("--docker_tag", default = "cpp")
parser.add_option("-e","--email", default = "aladyeva.e@wustl.edu")
parser.add_option("--restarts", type=int, default = 1)


(options, args) = parser.parse_args()

DT_STAMP = options.dt
if DT_STAMP is None:
    DT_STAMP = dt.datetime.now().strftime('%Y%m%d_%H%M%S')
    
Path(options.snakemake_path).mkdir(parents=True, exist_ok=True)
Path(options.reports_path).mkdir(parents=True, exist_ok=True)

with open(os.path.join(options.reports_path,"{0}.html".format(DT_STAMP)),"w+") as report:
    report.writelines("""<html>
    <body>""")
    report.writelines("Dataset: {0}<br>".format(options.dataset))
    report.writelines("Analysis: {0}<br>".format(options.analysis_name))
    report.writelines("Report ID: {0}<br><br>".format(DT_STAMP))

    for ct in range(options.min_ct,options.max_ct+1):
        WORK_DIR = os.path.join(options.snakemake_path,"ct{}".format(ct))
        Path(os.path.join(WORK_DIR,"resources")).mkdir(parents=True, exist_ok=True)
        Path(os.path.join(WORK_DIR,"config")).mkdir(parents=True, exist_ok=True)

        if not os.path.islink(os.path.join(WORK_DIR,"resources","datasets")):
            os.symlink(options.data_path, os.path.join(WORK_DIR,"resources","datasets"))

        Path(options.inits_path).mkdir(parents=True, exist_ok=True)
        if not os.path.islink(os.path.join(WORK_DIR,"resources","inits")):
            os.symlink(options.inits_path, os.path.join(WORK_DIR,"resources","inits"))

        Path(os.path.join(options.results_path,"ct{}".format(ct))).mkdir(parents=True, exist_ok=True)
        if not os.path.islink(os.path.join(WORK_DIR,"results")):
            os.symlink(os.path.join(options.results_path,"ct{}".format(ct)), os.path.join(WORK_DIR,"results"))

        Path(os.path.join(options.reports_path,"ct{}".format(ct))).mkdir(parents=True, exist_ok=True)
        if not os.path.islink(os.path.join(WORK_DIR,"reports")):
            os.symlink(os.path.join(options.reports_path,"ct{}".format(ct)), os.path.join(WORK_DIR,"reports"))

        config_dict = {}
        config_dict['num_inits']=options.num_inits
        if not options.top_mad is None:
            config_dict['top_mad']=options.top_mad
        if not options.top_median is None:
            config_dict['top_median']=options.top_median
        if not options.min_mad is None:
            config_dict['min_mad']=options.min_mad
        if not options.max_mad is None:
            config_dict['max_mad']=options.max_mad
        if not options.min_median is None:
            config_dict['min_median']=options.min_median
        if not options.max_median is None:
            config_dict['max_median']=options.max_median
        if not options.thresh is None:
            config_dict['thresh']=options.thresh
        if not options.r_tilda is None:
            config_dict['r_tilda']=options.r_tilda
        if options.no_cutoff_samples:
            config_dict['cutoff_samples']=not options.no_cutoff_samples
        if options.no_cutoff_genes:
            config_dict['cutoff_genes']=not options.no_cutoff_genes
        if options.no_filter_by_plane:
            config_dict['filter_by_plane']=not options.no_filter_by_plane
        config_dict['scale_iterations']=options.scale_iterations
        config_dict['init_strategy']=options.init_strategy
        config_dict['cell_types']=ct
        config_dict['preprocessing_cell_types']=options.preprocessing_cell_types
        config_dict['dt']=DT_STAMP
        config_dict['dataset']=options.dataset
        config_dict['analysis_name']="{0}_ct{1}".format(options.analysis_name,ct)
        config_dict['blocks_pipeline']="config/blocks.csv"
        config_dict['count']={'time':6000,'mem_ram':32,'threads':8,'email':options.email,
                              'nodes':-1,'docker':"{0}:{1}".format(options.docker_image,options.docker_tag)}
        with open(os.path.join(WORK_DIR,"config",'config.yaml'), 'w+') as ff:
            yaml.dump(config_dict, ff, allow_unicode=True, default_flow_style=False)

        shutil.copyfile(options.blocks_path,os.path.join(WORK_DIR,"config","blocks.csv"))
        shutil.copyfile(options.snakefile_path,os.path.join(WORK_DIR,"Snakefile"))
        os.chdir(WORK_DIR)
        print(os.getcwd())
        if options.local:
            cmd = ["snakemake","-c4","--rerun-incomplete"]
            print(cmd)
            p = Popen(cmd, shell=False, stdin=None, stdout=None, stderr=None, close_fds=True)
        else:
            snk_cmd = "snakemake --profile lsf_demo --local-cores $L_CORES --jobs 50 -pr --conda-frontend conda --restart-times {0} --rerun-incomplete".format(options.restarts)
            if options.apply_filters:
                snk_cmd += " -f apply_filters"
            cmd = """ P_WD=`pwd`; mkdir -p "$P_WD/tmp";
                echo "__LSF_JOB_CUSTOM_TMPDIR__=$P_WD/tmp" > lsf_docker_env_file.env;
                chmod a+r lsf_docker_env_file.env; export SMK_DOCKER_IMG="{1}:{2}";
                export P_LOG=$P_WD/logs/pipeline.log; export L_CORES=4; export LSF_DOCKER_ENV_FILE=$P_WD/lsf_docker_env_file.env; mkdir -p logs;
                bsub -cwd $HOME -n $L_CORES -G compute-martyomov -q general -oo $P_LOG -R 'span[hosts=1]' -a "docker($SMK_DOCKER_IMG)" /usr/bin/script -fqe /dev/null  -c "source /etc/bash.bashrc; cd $P_WD; export TMPDIR=$P_WD/tmp; export __LSF_JOB_CUSTOM_TMPDIR__=$P_WD/tmp; {0}" """.format(snk_cmd,
                options.docker_image,
                options.docker_tag)
            print(cmd)
            p = Popen(cmd, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
            
        report.writelines("<a href='./ct{1}/{0}.html'>{1} cell types</a><br>".format(DT_STAMP,ct))

    report.writelines("""</body>
    </html>""")