FROM dockerreg01.accounts.ad.wustl.edu/artyomov_lab/linseed_v2:nnls_d_cpp

RUN R -e "install.packages('rmarkdown')"
RUN R -e "install.packages('uwot')"
RUN R -e "install.packages('ggpubr')"

RUN  apt-get update \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# Install snakemake and other packages
COPY environment.yaml /app/environment.yaml
RUN /opt/conda/bin/conda update  --yes -n base -c defaults conda setuptools
RUN /opt/conda/bin/conda env update -n base --file /app/environment.yaml
RUN /opt/conda/bin/conda clean   --yes --all
RUN rm /app/environment.yaml

# Solve locale issues when running bash.
#   /bin/bash: warning: setlocale: LC_ALL: cannot change locale (en_US.UTF-8)
#
# It breaks conda version check in snakemake:
RUN apt-get clean && apt-get update && apt-get install -y locales && \
    echo "LC_ALL=en_US.UTF-8" >> /etc/environment  && \
    echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen  && \
    echo "LANG=en_US.UTF-8" > /etc/locale.conf  && \
    locale-gen en_US.UTF-8

COPY app /app
RUN chmod +x /app/scripts/run_linseedv2.py