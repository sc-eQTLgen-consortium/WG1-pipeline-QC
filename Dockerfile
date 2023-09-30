
################## BASE IMAGE ######################

FROM ubuntu:22.04

################## METADATA ######################

LABEL base_image="ubuntu:22.04"
LABEL version="1.0.0"
LABEL software="WG1 Pipeline"
LABEL about.summary="WG1 sceQTLGen Consortium Imputation, Demulitplexing, and Doublet Detection Pipeline"
LABEL about.documentation="https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC"
LABEL about.tags="Genomics"

################## MAINTAINER ######################

MAINTAINER Drew Neavin <d.neavin@garvan.org.au>, Martijn Vochteloo <m.vochteloo@umcg.nl>

################## INSTALLATION ######################

ADD . /tmp/repo
WORKDIR /tmp/repo
ENV PATH=/opt:/usr/games:/opt/conda/envs/py36/bin:/opt/conda/bin:/opt/minimap2-2.26:/opt/bedtools2/bin:/opt/.cargo/bin:/opt/souporcell:/opt/souporcell/troublet/target/release:/opt/vartrix-1.1.22:/opt/freebayes-1.3.7:/opt/freebayes-1.3.7/scripts:/opt/popscle/bin:/opt/DoubletDetection:/opt/Eagle_v2.4.1:/opt/Minimac4:/opt/GenotypeHarmonizer-1.4.27:/opt/picard-3.1.0/build/libs:$PATH
ENV SHELL=/bin/bash
ENV LC_ALL=C
ENV LANG=C.UTF-8
ENV TZ=Europe

# Needed to prevent asking for geagrpahic location when installing things.
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update -y
# libc-bin libc6 libsystemd0 libudev1
RUN apt-get upgrade -y

# binutils, binutils-common, binutils-x86-64-linux-gnu, build-essential, bzip2, cpp
# cpp-11, dirmngr, dpkg-dev, fakeroot, fontconfig-config, fonts-dejavu-core, g++
# g++-11, gcc, gcc-11, gcc-11-base, gnupg, gnupg-l10n, gnupg-utils, gpg, gpg-agent
# gpg-wks-client, gpg-wks-server, gpgconf, gpgsm, libalgorithm-diff-perl
# libalgorithm-diff-xs-perl, libalgorithm-merge-perl, libasan6, libassuan0
# libatomic1, libbinutils, libbrotli1, libbsd0, libc-dev-bin, libc-devtools
# libc6-dev, libcc1-0, libcrypt-dev, libctf-nobfd0, libctf0, libdeflate0
# libdpkg-perl, libexpat1, libfakeroot, libfile-fcntllock-perl, libfontconfig1
# libfreetype6, libgcc-11-dev, libgd3, libgdbm-compat4, libgdbm6, libgomp1, libisl23
# libitm1, libjbig0, libjpeg-turbo8, libjpeg8, libksba8, libldap-2.5-0
# libldap-common, liblocale-gettext-perl, liblsan0, libmd0, libmpc3, libmpfr6
# libnpth0, libnsl-dev, libperl5.34, libpng16-16, libquadmath0, libreadline8
# libsasl2-2, libsasl2-modules, libsasl2-modules-db, libsqlite3-0
# libstdc++-11-dev, libtiff5, libtirpc-dev, libtsan0, libubsan1, libwebp7, libx11-6
# libx11-data, libxau6, libxcb1, libxdmcp6, libxpm4, linux-libc-dev
# lto-disabled-list, make, manpages, manpages-dev, netbase, patch, perl
# perl-modules-5.34, pinentry-curses, readline-common, rpcsvc-proto, ucf, xz-utils
RUN apt-get install -y build-essential # includes g++, gcc, make, libc6-dev, and dpkg-dev
RUN apt-get install -y wget

##################################
############# PYTHON #############
##################################

# Install Python. Also required for bedtools2.
RUN apt-get install -y python3
#RUN cd /opt \
#    && wget https://www.python.org/ftp/python/3.6.15/Python-3.6.15.tgz \
#    && tar -xzf Python-3.6.15.tgz \
#    && rm Python-3.6.15.tgz \
#    && cd Python-3.6.15 \
#    && ./configure \
#    && make \
#    && make install

# Install miniconda for the virtual environment
RUN cd /opt \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-py39_23.5.2-0-Linux-x86_64.sh \
    && /bin/bash Miniconda3-py39_23.5.2-0-Linux-x86_64.sh -b -p /opt/conda \
    && rm Miniconda3-py39_23.5.2-0-Linux-x86_64.sh

# Create and activate virtual environment
# TODO this does not work.
# RUN eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
# RUN export PATH=/opt/conda/bin/:$PATH
# ca-certificates-2023.08.22, libuuid-1.41.5, openssl-3.0.11, pip-23.2.1, python-3.10.13, setuptools-68.0.0, wheel-0.41.2
RUN conda create -n py310 python=3.10
# RUN conda init zsh
# RUN eval "$(conda shell.bash hook)"
# RUN conda activate py310

# Downgrade setuptools for pyvcf.
RUN /opt/conda/envs/py310/bin/pip install --force-reinstall -v "setuptools==58.0.1"

# MarkupSafe-2.1.3 appdirs-1.4.4 attrs-23.1.0 certifi-2023.7.22 charset-normalizer-3.2.0 configargparse-1.7
# connection-pool-0.0.3 datrie-0.8.2 docutils-0.20.1 dpath-2.1.6 fastjsonschema-2.18.0 gitdb-4.0.10
# gitpython-3.1.37 humanfriendly-10.0 idna-3.4 jinja2-3.1.2 jsonschema-4.19.1 jsonschema-specifications-2023.7.1
# jupyter-core-5.3.2 nbformat-5.9.2 packaging-23.1 plac-1.4.0 platformdirs-3.10.0 psutil-5.9.5 pulp-2.7.0
# pyyaml-6.0.1 referencing-0.30.2 requests-2.31.0 reretry-0.11.8 rpds-py-0.10.3 smart-open-6.4.0 smmap-5.0.1
# snakemake-7.32.4 stopit-1.1.2 tabulate-0.9.0 throttler-1.2.2 toposort-1.10 traitlets-5.10.1 urllib3-2.0.5
# wrapt-1.15.0 yte-1.5.1
RUN /opt/conda/envs/py310/bin/pip install snakemake==7.32.4

# Requirements for souporcell.
# None
RUN /opt/conda/envs/py310/bin/pip install numpy==1.26.0
# absl-py-2.0.0 astunparse-1.6.3 cachetools-5.3.1 flatbuffers-23.5.26 gast-0.5.4 google-auth-2.23.2
# google-auth-oauthlib-1.0.0 google-pasta-0.2.0 grpcio-1.58.0 h5py-3.9.0 keras-2.14.0 libclang-16.0.6
# markdown-3.4.4 ml-dtypes-0.2.0 oauthlib-3.2.2 opt-einsum-3.3.0 protobuf-4.24.3 pyasn1-0.5.0
# pyasn1-modules-0.3.0 requests-oauthlib-1.3.1 rsa-4.9 six-1.16.0 tensorboard-2.14.1
# tensorboard-data-server-0.7.1 tensorflow-2.14.0 tensorflow-estimator-2.14.0 tensorflow-io-gcs-filesystem-0.34.0
# termcolor-2.3.0 typing-extensions-4.8.0 werkzeug-2.3.7 wrapt-1.14.1
RUN /opt/conda/envs/py310/bin/pip install tensorflow==2.14.0
# None
RUN /opt/conda/envs/py310/bin/pip install pyvcf==0.6.8
# aiohttp-3.8.5 aiosignal-1.3.1 async-timeout-4.0.3 clikit-0.6.2 crashtest-0.3.1 frozenlist-1.4.0
# httpstan-4.10.1 marshmallow-3.20.1 multidict-6.0.4 pastel-0.2.1 pylev-1.4.0 pysimdjson-5.0.2
# pystan-3.7.0 webargs-8.3.0 yarl-1.9.2
RUN /opt/conda/envs/py310/bin/pip install pystan==3.7.0
# importlib-metadata-6.8.0 pyfaidx-0.7.2.2 zipp-3.17.0
RUN /opt/conda/envs/py310/bin/pip install pyfaidx==0.7.2.2
# None
RUN /opt/conda/envs/py310/bin/pip install scipy==1.11.3

# Requirements for Scrublet.
# pandas-2.1.1 python-dateutil-2.8.2 pytz-2023.3.post1 tzdata-2023.3
RUN /opt/conda/envs/py310/bin/pip install pandas==2.1.1
# contourpy-1.1.1 cycler-0.12.0 fonttools-4.43.0 kiwisolver-1.4.5 matplotlib-3.8.0 pillow-10.0.1 pyparsing-3.1.1
RUN /opt/conda/envs/py310/bin/pip install matplotlib==3.8.0
# PyWavelets-1.4.1 annoy-1.17.3 cython-3.0.2 imageio-2.31.4 joblib-1.3.2 lazy_loader-0.3 llvmlite-0.41.0
# networkx-3.1 numba-0.58.0 numpy-1.25.2 pynndescent-0.5.10 scikit-image-0.21.0 scikit-learn-1.3.1
# scrublet-0.2.3 tbb-2021.10.0 threadpoolctl-3.2.0 tifffile-2023.9.26 tqdm-4.66.1 umap-learn-0.5.4
RUN /opt/conda/envs/py310/bin/pip install scrublet==0.2.3

# Requirements for CrossMap.
# CrossMap-0.6.6 bx-python-0.10.0 pyBigWig-0.3.22 pysam-0.21.0
RUN /opt/conda/envs/py310/bin/pip install CrossMap==0.6.6

# Requirements for DoubletDetection.
# anndata-0.9.2 asttokens-2.4.0 backcall-0.2.0 comm-0.1.4 decorator-5.1.1 doubletdetection-4.2
# exceptiongroup-1.1.3 executing-1.2.0 igraph-0.10.8 ipython-8.16.0 ipywidgets-8.1.1 jedi-0.19.0
# jupyterlab-widgets-3.0.9 leidenalg-0.10.1 louvain-0.8.1 matplotlib-inline-0.1.6 natsort-8.4.0
# parso-0.8.3 patsy-0.5.3 pexpect-4.8.0 phenograph-1.5.7 pickleshare-0.7.5 prompt-toolkit-3.0.39
# ptyprocess-0.7.0 pure-eval-0.2.2 pygments-2.16.1 scanpy-1.9.5 seaborn-0.12.2 session-info-1.0.0
# stack-data-0.6.2 statsmodels-0.14.0 stdlib_list-0.9.0 texttable-1.6.7 wcwidth-0.2.7 widgetsnbextension-4.0.9
RUN /opt/conda/envs/py310/bin/pip install doubletdetection==4.2

RUN conda clean --all

#############################
############# R #############
#############################

RUN apt-get install -y checkinstall # cmake
# cmake cmake-data dh-elpa-helper emacsen-common libarchive13 libcurl4 libicu70 libjsoncpp25 librhash0 libuv1 libxml2
RUN apt-get install -y cmake # ggpubr - nloptr.

# install two helper packages we need: dirmngr (installed by build-essential) and software-properties-common
RUN apt-get install -y --no-install-recommends software-properties-common
# add the signing key (by Michael Rutter) for these repos
# To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN wget -qO- "https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc"
RUN tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/"
# install R and its dependencies.
RUN apt-get install -y --no-install-recommends r-base
RUN apt-get install -y --no-install-recommends r-base-dev

# This needs to be in front of installing the R package tidyverse due to dependencies systemfonts, curl, xml2, openssl, and textshaping.
RUN apt-get install -y libfontconfig1-dev
RUN apt-get install -y libcurl4-openssl-dev
RUN apt-get install -y libxml2-dev
RUN apt-get install -y libssl-dev
RUN apt-get install -y libharfbuzz-dev
RUN apt-get install -y libfribidi-dev

# None
RUN R --slave -e 'install.packages("remotes")'
# jsonlite -> 1.8.7, findpython 1.0.8, R6 -> 2.5.1
RUN R --slave -e 'remotes::install_version("argparse", version = "2.2.2")'
# parallelly -> 1.36.0, listenv -> 0.9.0, globals -> 0.16.2, digest -> 0.6.33, future -> 1.33.0
RUN R --slave -e 'remotes::install_version("future.apply", version = "1.11.0")'
# None
RUN R --slave -e 'remotes::install_version("RColorBrewer", version = "1.1-3")'
# SnowballC -> 0.7.1
RUN R --slave -e 'remotes::install_version("lsa", version = "0.73.3")'
# None
RUN R --slave -e 'remotes::install_version("data.table", version = "1.14.8")'
# None
RUN R --slave -e 'remotes::install_version("xgboost", version = "1.7.5.1")' ### scDblFinder, depends on data.table
# vctrs -> 0.6.3, utf8 -> 1.2.3, rlang -> 1.1.1, lifecycle 1.0.3, glue -> 1.6.2, fansi -> 1.0.4
# cli -> 3.6.1, withr -> 2.5.1, pkgconfig -> 2.0.3, pillar -> 1.9.0, magrittr -> 2.0.3
# tidyselect -> 1.2.0, tibble -> 3.2.1, generics -> 0.1.3
RUN R --slave -e 'remotes::install_version("dplyr", version = "1.1.3")'
# stringi -> 1.7.12, cpp11 -> 0.4.6, stringr -> 1.5.0, purr -> 1.0.2
RUN R --slave -e 'remotes::install_version("tidyr", version = "1.3.0")' # depends on dplyr
# permute -> 0.9-7, Rcpp -> 1.0.11, viridisLite -> 0.4.2, vegan -> 2.6-4, pinfsc50 -> 1.2.0
# memuse -> 4.2-3, ape -> 5.7-1
RUN R --slave -e 'remotes::install_version("vcfR", version = "1.14.0")' # depends on dplyr
# colorspace -> 2.1-0, munsell -> 0.5.0, labeling -> 0.4.3, farver -> 2.1.1, scales -> 1.2.1
# isoband -> 0.2.7, gtable -> 0.3.4
RUN R --slave -e 'remotes::install_version("ggplot2", version = "3.4.3")'
# None
RUN R --slave -e 'remotes::install_version("ComplexUpset", version = "1.3.3")' # depends on ggplot2
# prettyunits -> 1.2.0, rematch2 -> 2.1.2, diffobj -> 0.3.5, pkgbuild -> 1.4.2, fs -> 1.6.3, crayon -> 1.5.2
# rprojroot -> 2.0.3, waldo -> 0.5.1, ps -> 1.7.5, processx -> 3.8.2, praise -> 1.0.0, pkgload -> 1.3.3
# evaluate -> 0.22, desc -> 1.4.2, callr -> 3.7.3, brio -> 1.1.3, testthat -> 3.1.10 ellipsis -> 0.3.2
# backports -> 1.4.1, MatrixModels (NA -> 0.5-2, SparseM -> 1.81, numDeriv -> 2016.8-1.1, broom -> 1.0.5
# RcppEigen -> 0.3.3.9.3, nloptr -> 2.0.3, minqa -> 1.2.6, lme4 -> 1.1-34, quantreg -> 5.97, pbkrtest -> 0.5.2
# abind -> 1.4-5, carData-> 3.0-5, car -> 3.1-2, corrplot -> 0.92, rstatix -> 0.7.2, polynom -> 1.4-1
# gridExtra -> 2.3, ggsignif -> 0.6.4, cowplot -> 1.1.1, ggsci -> 3.0.0, ggrepel -> 0.9.3
RUN R --slave -e 'remotes::install_version("ggpubr", version = "0.6.0")' # depends on ggplot2
# None
RUN R --slave -e 'remotes::install_version("cowplot", version = "1.1.1")' # depends on ggplot2
# hms -> 1.1.3, bit-> 4.0.5, progress-> 1.2.2, tzdb -> 0.4.0, bit64 -> 4.0.5, vroom -> 1.6.3, clipr -> 0.8.0)
RUN R --slave -e 'remotes::install_version("readr", version = "2.1.4")'
# rappdirs -> 0.3.3, highr -> 0.10, fastmap -> 1.1.1, sass -> 0.4.7, mime -> 0.12, memoise -> 2.0.1
# cachem -> 1.0.8, base64enc -> 0.1-3, yaml -> 2.3.7, xfun -> 0.40, tinytex -> 0.46, knitr -> 1.44
# jquerylib -> 0.1.4, htmltools -> 0.5.6, fontawesome -> 0.5.2, bslib -> 0.5.1, systemfonts -> 1.0.4
# sys -> 3.4.2, askpass -> 1.2.0, uuid -> 1.1-1, openssl -> 2.1.1, curl -> 5.0.2, httr -> 1.4.7
# gargle -> 1.5.2, rematch -> 2.0.0, xml2 -> 1.3.5, selectr -> 0.4-2, rstudioapi -> 0.15.0, rmarkdown -> 2.25
# cellranger -> 1.1.0, textshaping -> 0.3.6, timechange -> 0.2.0, forcats -> 1.0.0, ids -> 1.0.1
# googledrive -> 2.1.1, DBI -> 1.1.3, blob -> 1.2.4, rvest -> 1.0.3, reprex -> 2.0.2, readxl -> 1.4.3
# ragg -> 1.2.5, modelr -> 0.1.11, lubridate -> 1.9.3, haven -> 2.5.3, googlesheets4 -> 1.1.1, dtplyr -> 1.3.1
# dbplyr -> 2.3.4, conflicted -> 1.2.0
RUN R --slave -e 'remotes::install_version("tidyverse", version = "2.0.0")' # depends on dplyr, tidyr, ggplot2, readr
# bitops -> 1.0-7, caTools -> 1.18.2, gtools -> 3.9.4, dotCall64 -> 1.0-2, gplots -> 3.1.3
# maps -> 3.4.1, spam -> 2.9-1, ROCR -> 1.0-11, fields -> 15.2
RUN R --slave -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder@1b1d4e2d7f893a3552d9f8f791ab868ee4c782e6")' # Commit of Aug 18, 2023
# progressr -> 0.14.0, sp -> 2.0-0
RUN R --slave -e 'remotes::install_version("SeuratObject", version = "4.1.4")' # depends on future.apply
# sitmo -> 2.0.2, BH -> 1.81.0-1, spatstat.utils -> 3.0-3, tensor -> 1.5, polyclip -> 1.10-6, deldir -> 1.0-9
# spatstat.geom -> 3.2-5, spatstat.data -> 3.0-1, later -> 1.3.1, promises -> 1.2.1, plyr -> 1.8.8
# lazyeval -> 0.2.2, commonmark -> 1.9.0, sourcetools -> 0.1.7-1, xtable -> 1.8-4, httpuv -> 1.6.11
# png -> 0.1-8, here -> 1.0.1, RcppTOML -> 0.2.2, dqrng -> 0.3.1, RcppProgress -> 0.4.2, irlba -> 2.3.5.1
# RcppAnnoy -> 0.0.21, FNN -> 1.1.3.2, goftest -> 1.2-3, spatstat.sparse -> 3.0-2, spatstat.random -> 3.1-6
# RcppArmadillo -> 0.12.6.4.0, matrixStats -> 1.0.0, reshape2 -> 1.4.4, crosstalk -> 1.2.0, htmlwidgets -> 1.6.2
# shiny -> 1.7.5, zoo -> 1.8-12, igraph -> 1.5.1, reticulate -> 1.32.0, uwot -> 0.1.16, spatstat.explore -> 3.2-3
# sctransform -> 0.4.0, scattermore -> 1.2, Rtsne -> 0.16, RANN -> 2.6.1, plotly -> 4.10.2, pbapply -> 1.7-2
# miniUI -> 0.1.1.1, lmtest -> 0.9-40, leiden -> 0.4.3, ica -> 1.0-3, ggridges -> 0.5.4, fitdistrplus -> 1.1-11
RUN R --slave -e 'remotes::install_version("Seurat", version = "4.4.0")' # depends on future.apply, RColorBrewer, ggplot2, cowplot, SeuratObject
# textshaping -> 0.3.6, vipor -> 0.4.5, beeswarm -> 0.4.0, prismatic -> 1.1.1, snakecase -> 0.11.1
# ragg -> 1.2.5, ggbeeswarm -> 0.7.2, Cairo -> 1.6-1, shape -> 1.4.6, GlobalOptions -> 0.1.2, paletteer -> 1.5.0
# janitor -> 2.2.0, ggrastr -> 1.0.2, ggprism -> 1.0.4, circlize -> 0.4.15
RUN R --slave -e 'remotes::install_version ("scCustomize", version = "1.1.3")' # depends on Seurat

RUN R --slave -e 'install.packages("BiocManager")'
RUN R --slave -e 'BiocManager::install(version = "3.14")' # for R version 4.1
# rjson, BiocGenerics, S4Vectors, iterators, GetoptLong, clue, IRanges, foreach, doParallel
RUN R --slave -e 'BiocManager::install("ComplexHeatmap", version = "3.14")' # 2.16.0
# ????
RUN R --slave -e 'BiocManager::install("SingleCellExperiment", version = "3.14")' # 1.22.0
# locfit, sparseMatrixStats, RcppHNSW, ScaledMatrix, rsvd, edgeR, limma, statmod, DelayedMatrixStats
# metapod, viridis, BiocNeighbors, BiocSingular, scran, scater, scuttle, bluster
RUN R --slave -e 'BiocManager::install("scDblFinder", version = "3.14")' # Depends on SingleCellExperiment, 1.14.0
# pROC
RUN R --slave -e 'BiocManager::install("scds", version = "3.14")' # 1.16.0

################################
############# JAVA #############
################################

RUN apt-get install -y openjdk-17-jdk

#################################
############# OTHER #############
#################################

## This is only needed if you do not run the Python install section.
# RUN apt-get install -y python3 # bedtools2

## This is only needed if you do not run the Java install section.
# RUN apt-get install -y openjdk-17-jdk # picard

# zlib1g-dev
RUN apt-get install -y libz-dev # bedtools2, samtools, bcftools, minimap2, vcftools, htslib
# bzip2-doc libbz2-dev
RUN apt-get install -y libbz2-dev # bedtools2, samtools, bcftools, htslib
# liblzma-doc
RUN apt-get install -y liblzma-dev # bedtools2, samtools, bcftools, htslib
# libncurses-dev
RUN apt-get install -y libncurses-dev # samtools
# curl libcurl4
RUN apt-get install -y curl # cargo for souporcell, also Google Cloud Storage support for htslib
# git git-man less libcbor0.8 libcurl3-gnutls libedit2 liberror-perl libfido2-1 libnghttp2-14 librtmp1 libssh-4 libxext6
# libxmuu1 openssh-client xauth
RUN apt-get install -y git # souporcell, popscle
# libglib2.0-0 libglib2.0-data libicu70 libxml2 pkg-config shared-mime-info xdg-user-dirs
RUN apt-get install -y pkg-config # vcftools

## This is only needed if you do not run the R install section.
## libcurl4 libcurl4-openssl-dev
#RUN apt-get install -y libcurl4-openssl-dev # htslib
## cmake cmake-data dh-elpa-helper emacsen-common libarchive13 libcurl4 libicu70 libjsoncpp25 librhash0 libuv1 libxml2
#RUN apt-get install -y cmake # popscle
## libssl-dev
#RUN apt-get install -y libssl-dev # popscle
## unzip
#RUN apt-get install -y unzip # plink, plink2

# Requires Python to install.
RUN cd /opt \
    && wget https://github.com/arq5x/bedtools2/archive/refs/tags/v2.31.0.tar.gz \
    && tar -xzf v2.31.0.tar.gz \
    && rm v2.31.0.tar.gz \
    && cd bedtools2-2.31.0 \
      && make

RUN cd /opt \
    && wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 \
    && tar xjf samtools-1.18.tar.bz2 \
    && rm samtools-1.18.tar.bz2 \
    && cd samtools-1.18 \
      && ./configure \
      && make \
      && make install

RUN cd /opt \
    && wget https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2 \
    && tar xjf bcftools-1.18.tar.bz2 \
    && rm bcftools-1.18.tar.bz2 \
    && cd bcftools-1.18 \
      && ./configure \
      && make \
      && make install

RUN cd /opt \
    && wget https://github.com/freebayes/freebayes/archive/refs/tags/v1.3.7.tar.gz \
    && tar -xzf v1.3.7.tar.gz \
    && rm v1.3.7.tar.gz

RUN cd /opt \
    && wget https://github.com/10XGenomics/vartrix/archive/refs/tags/v1.1.22.tar.gz \
    && tar -xzf v1.1.22.tar.gz \
    && rm v1.1.22.tar.gz

RUN cd /opt \
    && wget https://github.com/lh3/minimap2/archive/refs/tags/v2.26.tar.gz \
    && tar -xzf v2.26.tar.gz \
    && rm v2.26.tar.gz \
    && cd minimap2-2.26 \
      && make

# Required for souporcell.
# TODO: specify version
# info: default host triple is x86_64-unknown-linux-gnu
# info: syncing channel updates for 'stable-x86_64-unknown-linux-gnu'
# info: latest update on 2023-09-19, rust version 1.72.1 (d5c2e9c34 2023-09-13)
RUN cd /opt \
    && CARGO_HOME=/opt/.cargo RUSTUP_HOME=/opt/.cargo bash -c 'curl https://sh.rustup.rs -sSf | sh -s -- -y' \
    && . /opt/.cargo/env \
    && rustup default stable

# Fork from commit 9fb5271ae9f2257ea9a8552dfda3d4b7080be194
# requires samtools, bcftools, htslib, python3, freebayes, vartrix, minimap2, Rust (i.e. cargo)
RUN cd /opt \
    && git clone --single-branch --branch RoyOelen/add-fastq-gzipping https://github.com/sc-eQTLgen-consortium/souporcell.git \
    && cd /opt/souporcell/souporcell \
      && cargo build --release \
    && cd /opt/souporcell/troublet \
      && cargo build --release

# RUN apt-get install -y minimac4=4.1.4
RUN cd /opt \
      && wget https://github.com/statgen/Minimac4/releases/download/v4.1.4/minimac4-4.1.4-Linux-x86_64.sh \
      && bash minimac4-4.1.4-Linux-x86_64.sh --skip-license TRUE

RUN cd /opt \
    && wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz \
    && tar -xzf vcftools-0.1.16.tar.gz \
    && rm vcftools-0.1.16.tar.gz \
    && cd vcftools-0.1.16 \
      && ./configure \
      && make \
      && make install

# Required for popscle.
RUN cd /opt \
    && wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 \
    && tar xjf htslib-1.18.tar.bz2 \
    && rm htslib-1.18.tar.bz2 \
    && cd htslib-1.18 \
      && ./configure \
      && make \
      && make install

## Commit of May 5, 2021
## TODO fix error. Perhaps downgrade htslib to 1.10.2?
## [CMakeFiles/popscle.dir/build.make:230: CMakeFiles/popscle.dir/cmd_plp_make_dge_matrix.cpp.o] Error 1
## [CMakeFiles/Makefile2:83: CMakeFiles/popscle.dir/all] Error 2
#RUN cd /opt \
#    && git clone https://github.com/statgen/popscle.git \
#    && cd popscle \
#      && git checkout da70fc78da385ef049e0e890342acfd62842cae0 \
#      && mkdir build \
#      && cd build \
#        && cmake ..  \
#        && make

# Used by minimap2 and picard? No clue but seems important.
RUN cd /opt \
    && wget https://ftp.gnu.org/gnu/parallel/parallel-20230922.tar.bz2 \
    && tar xjf parallel-20230922.tar.bz2 \
    && rm parallel-20230922.tar.bz2 \
    && cd parallel-20230922 \
      && ./configure \
      && make \
      && make install

RUN cd /opt \
    && wget https://alkesgroup.broadinstitute.org/Eagle/downloads/Eagle_v2.4.1.tar.gz \
    && tar -xzf Eagle_v2.4.1.tar.gz \
    && rm Eagle_v2.4.1.tar.gz

RUN cd /opt \
    && wget https://github.com/molgenis/systemsgenetics/releases/download/GH_1.4.27/GenotypeHarmonizer-1.4.27-SNAPSHOT-dist.tar.gz \
    && tar -xzf GenotypeHarmonizer-1.4.27-SNAPSHOT-dist.tar.gz  \
    && rm GenotypeHarmonizer-1.4.27-SNAPSHOT-dist.tar.gz

RUN cd /opt \
    && mkdir plink \
    && cd plink \
      && wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip \
      && unzip plink_linux_x86_64_20230116.zip \
      && rm plink_linux_x86_64_20230116.zip

RUN cd /opt \
    && mkdir plink2 \
    && cd plink2 \
      && wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20230927.zip \
      && unzip plink2_linux_x86_64_20230927.zip \
      && rm plink2_linux_x86_64_20230927.zip

# Requires openjdk-17-jdk.
RUN cd /opt \
    && git clone https://github.com/broadinstitute/picard \
    && cd picard \
      && git checkout tags/3.1.0 \
      && ./gradlew shadowJar

# Always get our own newest software.
# TODO: make sure you use the correct branch here.
RUN cd /opt \
    && git clone --single-branch --branch scMetaBrain https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC.git \
    && chmod 777 -R /opt/WG1-pipeline-QC/Demultiplexing/scripts \
    && chmod 777 -R /opt/WG1-pipeline-QC/Imputation/scripts


## TODO: this stuff should probably not be included in the image but downloaded separately.
## md5sum 77e1141441b2b443519249d331dab5ae 1000G.tar.gz
#RUN cd /opt \
#    && wget https://www.dropbox.com/s/xso2vt3p9h2rh8m/1000G.tar.gz \
#    && tar -xzf 1000G.tar.gz \
#    && rm 1000G.tar.gz
#
## md5sum 2a862e70fa1afc08f1f3e7a8cc5cda14 GRCh37_to_GRCh38.chain.gz
## md5sum 79f905be8ef01dce1196975cab239cca GRCh37_to_GRCh38.chain
#RUN cd /opt \
#    && wget https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz \
#    && gunzip GRCh37_to_GRCh38.chain.gz
#
## md5sum cb4ba08996a10b0ded40fae00987d1c8 hg18ToHg38.over.chain.gz
## md5sum da54dd9e82393c38512de0e63335a05d hg18ToHg38.over.chain
#RUN cd /opt \
#    && wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz \
#    && gunzip hg18ToHg38.over.chain.gz
#
## md5sum c323d4fb5da690ca6352a7ee8a14e022 hg38exonsUCSC.bed
#RUN cd /opt \
#    && wget https://www.dropbox.com/s/fvd4pl8no3ngg0l/hg38exonsUCSC.bed
