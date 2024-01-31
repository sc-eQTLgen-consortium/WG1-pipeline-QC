
################## BASE IMAGE ######################

FROM ubuntu:22.04

################## METADATA ######################

LABEL base_image="ubuntu:22.04"
LABEL version="1.0.0"
LABEL software="WG1 Pipeline"
LABEL about.summary="WG1 sceQTLGen Consortium Imputation, Demulitplexing, and Doublet Detection Pipeline"
LABEL about.documentation="https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC"
LABEL about.tags="Genomics"

# Build syntax: docker build ./ -t wg1-pipeline-qc:2023.11.23.0 --progress=plain > build.log 2>&1
# Total build takes 1 hour, and 48 minutes and has a size of 7.41 GB.
# Use dive wg1-pipeline-qc:1.0.0 to investigate memory usage.

################## MAINTAINER ######################

MAINTAINER Drew Neavin <d.neavin@garvan.org.au>, Martijn Vochteloo <m.vochteloo@umcg.nl>

################## INSTALLATION ######################

# Section build takes 23 seconds and has a size of 0.469 GB.

# Uses 74 MB.
ADD . /tmp/repo
WORKDIR /tmp/repo

ENV PATH=/opt:/usr/games:/opt/conda/envs/py311/bin:/opt/conda/bin:/opt/minimap2-2.26:/opt/bedtools2-2.31.0/bin:/opt/.cargo/bin:/opt/souporcell:/opt/souporcell/souporcell/target/release:/opt/souporcell/troublet/target/release:/opt/vartrix-1.1.22:/opt/freebayes-1.3.7:/opt/freebayes-1.3.7/scripts:/opt/popscle/bin:/opt/DoubletDetection:/opt/Eagle_v2.4.1:/opt/bin/:/opt/GenotypeHarmonizer-1.4.27:/opt/plink:/opt/plink2:$PATH
ENV BCFTOOLS_PLUGINS=/opt/bcftools-1.18/plugins
ENV SHELL=/bin/bash
ENV LC_ALL=C
ENV LANG=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

RUN echo 'alias python=python3' >> ~/.bashrc

# Needed to prevent asking for geographic location when installing things.
# Uses 0.000007 MB.
RUN export TZ=Europe/Amsterdam \
    && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone

# Uses 300 MB, mainly in /usr/lib/ (221 MB), /usr/bin/ (59 MB), and /var/ (50 MB)
RUN apt-get update -y \
    # libc-bin libc6 libsystemd0 libudev1
    && apt-get upgrade -y \
    # binutils binutils-common binutils-x86-64-linux-gnu build-essential bzip2 cpp
    # cpp-11 dpkg-dev g++ g++-11 gcc gcc-11 gcc-11-base libasan6 libatomic1
    # libbinutils libc-dev-bin libc6-dev libcc1-0 libcrypt-dev libctf-nobfd0
    # libctf0 libdpkg-perl libgcc-11-dev libgdbm-compat4 libgdbm6 libgomp1
    # libisl23 libitm1 liblsan0 libmpc3 libmpfr6 libnsl-dev libperl5.34
    # libquadmath0 libstdc++-11-dev libtirpc-dev libtsan0 libubsan1 linux-libc-dev
    # lto-disabled-list make patch perl perl-modules-5.34 rpcsvc-proto xz-utils
    && apt-get install -y --no-install-recommends build-essential \
    # ca-certificates openssl
    && apt-get install -y --no-install-recommends ca-certificates \
    # libpsl5 wget
    && apt-get install -y --no-install-recommends wget

##################################
############# PYTHON #############
##################################

# Section build takes 3 minutes and 27 seconds and has a size of 3.38 GB.

# Reduce conda size by preventing Python from recreating a corresponding bytecode cache file (*.pyc) at runtime.
ENV PYTHONDONTWRITEBYTECODE=true

# Install Python. Also required for bedtools2.
# libexpat1 libmpdec3 libpython3-stdlib libpython3.10-minimal
# libpython3.10-stdlib libreadline8 libsqlite3-0 media-types python3
# python3-minimal python3.10 python3.10-minimal readline-common
# Uses 30 MB, mainly in /usr/lib/ (23 MB) and /usr/bin/ (6 MB)
RUN apt-get install -y --no-install-recommends python3

# Install miniconda for the virtual environment.
# https://github.com/ContinuumIO/docker-images/blob/main/miniconda3/debian/Dockerfile
# Uses 294 MB, mainly in /opt/conda/lib/ (231 MB), /opt/conda/bin/ (37 MB), and /opt/conda/share/ (8.8 MB)
RUN cd /opt \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-py311_23.5.2-0-Linux-x86_64.sh -O miniconda.sh -q \
    && mkdir -p /opt \
    && bash miniconda.sh -b -p /opt/conda \
    && rm miniconda.sh \
    && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
    && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
    && /opt/conda/bin/conda clean -afy

# Create and activate virtual environment
# _libgcc_mutex-0.1, _openmp_mutex-5.1, ca-certificates-2023.08.22, ld_impl_linux-64-2.38, libffi-3.4.4
# libgcc-ng-11.2.0, libgomp-11.2.0, libstdcxx-ng-11.2.0 libuuid-1.41.5, ncurses-6.4, openssl-3.0.11,
# pip-23.2.1, python-3.11.5, readline-8.2, setuptools-68.0.0, sqlite-3.41.2 tk-8.6.12, tzdata-2023c
# wheel-0.41.2, xz-5.4.2, zlib-1.2.13
# Uses 3.1 GB, mainly in /opt/conda/envs/py311/lib/ (2.9 GB)
RUN eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)" \
    && conda create -n py311 python=3.11.5 \
    # MarkupSafe-2.1.3 appdirs-1.4.4 attrs-23.1.0 certifi-2023.7.22 charset-normalizer-3.2.0 configargparse-1.7
    # connection-pool-0.0.3 datrie-0.8.2 docutils-0.20.1 dpath-2.1.6 fastjsonschema-2.18.0 gitdb-4.0.10
    # gitpython-3.1.37 humanfriendly-10.0 idna-3.4 jinja2-3.1.2 jsonschema-4.19.1 jsonschema-specifications-2023.7.1
    # jupyter-core-5.3.2 nbformat-5.9.2 packaging-23.1 plac-1.4.0 platformdirs-3.10.0 psutil-5.9.5 pulp-2.7.0
    # pyyaml-6.0.1 referencing-0.30.2 requests-2.31.0 reretry-0.11.8 rpds-py-0.10.3 smart-open-6.4.0 smmap-5.0.1
    # snakemake-7.32.4 stopit-1.1.2 tabulate-0.9.0 throttler-1.2.2 toposort-1.10 traitlets-5.10.1 urllib3-2.0.5
    # wrapt-1.15.0 yte-1.5.1
    && /opt/conda/envs/py311/bin/pip install snakemake==7.32.4 \
    # Requirements for souporcell.
    # None
    && /opt/conda/envs/py311/bin/pip install numpy==1.26.0 \
    # None
    && /opt/conda/envs/py311/bin/pip install scipy==1.11.3 \
    # Requirements for Scrublet.
    # pandas-2.1.1 python-dateutil-2.8.2 pytz-2023.3.post1 tzdata-2023.3
    && /opt/conda/envs/py311/bin/pip install pandas==2.1.1 \
    # contourpy-1.1.1 cycler-0.12.0 fonttools-4.43.0 kiwisolver-1.4.5 matplotlib-3.8.0 pillow-10.0.1 pyparsing-3.1.1
    && /opt/conda/envs/py311/bin/pip install matplotlib==3.8.0 \
    # PyWavelets-1.4.1 annoy-1.17.3 cython-3.0.2 imageio-2.31.4 joblib-1.3.2 lazy_loader-0.3 llvmlite-0.41.0
    # networkx-3.1 numba-0.58.0 numpy-1.25.2 pynndescent-0.5.10 scikit-image-0.21.0 scikit-learn-1.3.1
    # scrublet-0.2.3 tbb-2021.10.0 threadpoolctl-3.2.0 tifffile-2023.9.26 tqdm-4.66.1 umap-learn-0.5.4
    && /opt/conda/envs/py311/bin/pip install scrublet==0.2.3 \
    # Requirements for CrossMap.
    # CrossMap-0.6.6 bx-python-0.10.0 pyBigWig-0.3.22 pysam-0.21.0
    && /opt/conda/envs/py311/bin/pip install CrossMap==0.6.6 \
    # Requirements for DoubletDetection.
    # anndata-0.9.2 asttokens-2.4.0 backcall-0.2.0 comm-0.1.4 decorator-5.1.1 doubletdetection-4.2
    # exceptiongroup-1.1.3 executing-1.2.0 igraph-0.10.8 ipython-8.16.0 ipywidgets-8.1.1 jedi-0.19.0
    # jupyterlab-widgets-3.0.9 leidenalg-0.10.1 louvain-0.8.1 matplotlib-inline-0.1.6 natsort-8.4.0
    # parso-0.8.3 patsy-0.5.3 pexpect-4.8.0 phenograph-1.5.7 pickleshare-0.7.5 prompt-toolkit-3.0.39
    # ptyprocess-0.7.0 pure-eval-0.2.2 pygments-2.16.1 scanpy-1.9.5 seaborn-0.12.2 session-info-1.0.0
    # stack-data-0.6.2 statsmodels-0.14.0 stdlib_list-0.9.0 texttable-1.6.7 wcwidth-0.2.7 widgetsnbextension-4.0.9
    && /opt/conda/envs/py311/bin/pip install doubletdetection==4.2 \
    # No clue why we do this.
    && sed -i 's/louvain.set_rng_seed(random_state)/partition_kwargs["seed"] = random_state/g'  /opt/conda/envs/py311/lib/python3.11/site-packages/scanpy/tools/_louvain.py \
    && /opt/conda/envs/py311/bin/pip cache purge


# Creating a conda environment for souporcell using the specific versions as specified by the author.
# numpy==1.16.4 (previous image used; 1.19.2)
# scipy==1.3.0 (previous image used; 1.5.2)
# pystan==2.19.1.1 (previous image used; 2.17.1.0)
# pyvcf==0.6.8 (previous image used; ?? .__version__ not implemented)
# pysam==0.15.2 (previous image used; 0.16.0.1)
# pyfaidx=0.5.8=py_1 (previous image used; 0.5.9.1)
RUN eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)" \
    && wget https://raw.githubusercontent.com/wheaton5/souporcell/master/souporcell_env.yaml \
    && conda env create -f souporcell_env.yaml \
    && rm souporcell_env.yaml

RUN conda clean -y --all

#############################
############# R #############
#############################

# Section build takes 41 minutes, and 14 seconds and has a size of 1.381 GB.

# ggpubr - nloptr.
# cmake cmake-data dh-elpa-helper emacsen-common libarchive13 libbrotli1
# libcurl4 libexpat1 libicu70 libjsoncpp25 libldap-2.5-0 libnghttp2-14
# librhash0 librtmp1 libsasl2-2 libsasl2-modules-db libssh-4 libuv1 libxml2
# Uses 386 MB, mainly in /usr/lib/ (328 MB), /usr/bin/ (37 MB), /usr/share/ (30 MB)
RUN apt-get install -y --no-install-recommends cmake \
    # install two helper packages we need: dirmngr and software-properties-common
    && apt-get install -y --no-install-recommends dirmngr \
    # dbus distro-info-data gir1.2-glib-2.0 gir1.2-packagekitglib-1.0 gpg gpgconf
    # iso-codes libapparmor1 libappstream4 libargon2-1 libassuan0 libcap2-bin
    # libcryptsetup12 libcurl3-gnutls libdbus-1-3 libdevmapper1.02.1 libdw1
    # libelf1 libgirepository-1.0-1 libglib2.0-0 libglib2.0-bin libglib2.0-data
    # libgstreamer1.0-0 libip4tc2 libjson-c5 libkmod2 libmpdec3
    # libpackagekit-glib2-18 libpam-systemd libpolkit-agent-1-0
    # libpolkit-gobject-1-0 libpython3-stdlib libpython3.10-minimal
    # libpython3.10-stdlib libreadline8 libsqlite3-0 libstemmer0d libunwind8
    # libxmlb2 libyaml-0-2 lsb-release media-types packagekit pkexec policykit-1
    # polkitd python-apt-common python3 python3-apt python3-blinker
    # python3-cffi-backend python3-cryptography python3-dbus python3-distro
    # python3-gi python3-httplib2 python3-importlib-metadata python3-jeepney
    # python3-jwt python3-keyring python3-launchpadlib python3-lazr.restfulclient
    # python3-lazr.uri python3-minimal python3-more-itertools python3-oauthlib
    # python3-pkg-resources python3-pyparsing python3-secretstorage python3-six
    # python3-software-properties python3-wadllib python3-zipp python3.10
    # python3.10-minimal readline-common software-properties-common systemd
    # systemd-sysv
    && apt-get install -y --no-install-recommends software-properties-common \
    # add the signing key (by Michael Rutter) for these repos
    # To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
    # Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
    && wget -qO- "https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc" \
    && tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    # add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
    # Install R and its dependencies. \
    # fontconfig fontconfig-config fonts-dejavu-core libblas3 libbsd0 libcairo2
    # libdatrie1 libdeflate0 libfontconfig1 libfreetype6 libfribidi0 libgfortran5
    # libgraphite2-3 libharfbuzz0b libice6 libjbig0 libjpeg-turbo8 libjpeg8
    # liblapack3 libmd0 libpango-1.0-0 libpangocairo-1.0-0 libpangoft2-1.0-0
    # libpaper-utils libpaper1 libpixman-1-0 libpng16-16 libsm6 libtcl8.6
    # libthai-data libthai0 libtiff5 libtk8.6 libwebp7 libx11-6 libx11-data
    # libxau6 libxcb-render0 libxcb-shm0 libxcb1 libxdmcp6 libxext6 libxft2
    # libxrender1 libxss1 libxt6 r-base r-base-core r-cran-boot r-cran-class
    # r-cran-cluster r-cran-codetools r-cran-foreign r-cran-kernsmooth
    # r-cran-lattice r-cran-mass r-cran-matrix r-cran-mgcv r-cran-nlme r-cran-nnet
    # r-cran-rpart r-cran-spatial r-cran-survival r-recommended tzdata ucf unzip
    # x11-common xdg-utils zip
    && apt-get install -y --no-install-recommends r-base \
    # gfortran gfortran-11 icu-devtools libblas-dev libbz2-dev libgfortran-11-dev
    # libicu-dev libjpeg-dev libjpeg-turbo8-dev libjpeg8-dev liblapack-dev
    # liblzma-dev libncurses-dev libncurses5-dev libpcre16-3 libpcre2-16-0
    # libpcre2-32-0 libpcre2-dev libpcre2-posix3 libpcre3-dev libpcre32-3
    # libpcrecpp0v5 libpng-dev libreadline-dev libxmuu1 pkg-config r-base-dev
    # xauth zlib1g-dev
    && apt-get install -y --no-install-recommends r-base-dev \
    # This needs to be in front of installing the R package tidyverse due to dependencies
    # systemfonts, curl, xml2, openssl, and textshaping.
    # libbrotli-dev libexpat1-dev libfontconfig-dev libfontconfig1-dev
    # libfreetype-dev libfreetype6-dev uuid-dev
    && apt-get install -y --no-install-recommends libfontconfig1-dev \
    # libcurl4-openssl-dev
    && apt-get install -y --no-install-recommends libcurl4-openssl-dev \
    # libxml2-dev
    && apt-get install -y --no-install-recommends libxml2-dev \
    # libssl-dev
    && apt-get install -y --no-install-recommends libssl-dev \
    # gir1.2-harfbuzz-0.0 libblkid-dev libffi-dev libglib2.0-dev
    # libglib2.0-dev-bin libgraphite2-dev libharfbuzz-dev libharfbuzz-gobject0
    # libharfbuzz-icu0 libmount-dev libselinux1-dev libsepol-dev python3-distutils
    # python3-lib2to3
    && apt-get install -y --no-install-recommends libharfbuzz-dev \
    # libfribidi-dev
    && apt-get install -y --no-install-recommends libfribidi-dev \
    # libdeflate-dev libjbig-dev libtiff-dev libtiffxx5
    && apt-get install -y --no-install-recommends libtiff-dev \
    # Required for hdf5r
    && apt-get install -y --no-install-recommends libhdf5-dev \
    # Required for scCustomize
    # libcairo-gobject2 libcairo-script-interpreter2 libcairo2-dev libice-dev liblzo2-2
    # libpixman-1-dev libpthread-stubs0-dev libsm-dev libx11-dev libxau-dev libxcb-render0-dev
    # libxcb-shm0-dev libxcb1-dev libxdmcp-dev libxext-dev libxrender-dev x11proto-dev
    # xorg-sgml-doctools xtrans-dev
    && apt-get install -y --no-install-recommends libcairo2-dev

# remotes_2.4.2.1
# Uses 791 MB, mainly in /usr/local/lib/R/BH/ (130 MB), xgboost (100 MB), RcppEigen (37 MB), lme4 (29MB), igraph (20 MB), vroom (20 MB), uwot (20 MB)
RUN R --slave -e 'install.packages("remotes")' \
    # jsonlite -> 1.8.7, findpython 1.0.8, R6 -> 2.5.1
    && R --slave -e 'remotes::install_version("argparse", version = "2.2.2", upgrade=FALSE)' \
    # parallelly -> 1.36.0, listenv -> 0.9.0, globals -> 0.16.2, digest -> 0.6.33, future -> 1.33.0
    && R --slave -e 'remotes::install_version("future.apply", version = "1.11.0", upgrade=FALSE)' \
    # None
    && R --slave -e 'remotes::install_version("R.utils", version = "2.12.2", upgrade=FALSE)' \
    # None
    && R --slave -e 'remotes::install_version("RColorBrewer", version = "1.1-3", upgrade=FALSE)' \
    # SnowballC -> 0.7.1
    && R --slave -e 'remotes::install_version("lsa", version = "0.73.3", upgrade=FALSE)' \
    # None
    && R --slave -e 'remotes::install_version("data.table", version = "1.14.8", upgrade=FALSE)' \
    ### scDblFinder, depends on data.table
    # None
    && R --slave -e 'remotes::install_version("xgboost", version = "1.7.5.1", upgrade=FALSE)' \
    # vctrs -> 0.6.3, utf8 -> 1.2.3, rlang -> 1.1.1, lifecycle 1.0.3, glue -> 1.6.2, fansi -> 1.0.4
    # cli -> 3.6.1, withr -> 2.5.1, pkgconfig -> 2.0.3, pillar -> 1.9.0, magrittr -> 2.0.3
    # tidyselect -> 1.2.0, tibble -> 3.2.1, generics -> 0.1.3
    && R --slave -e 'remotes::install_version("dplyr", version = "1.1.3", upgrade=FALSE)' \
    # depends on dplyr
    # stringi -> 1.7.12, cpp11 -> 0.4.6, stringr -> 1.5.0, purr -> 1.0.2
    && R --slave -e 'remotes::install_version("tidyr", version = "1.3.0", upgrade=FALSE)' \
    # depends on dplyr
    # permute -> 0.9-7, Rcpp -> 1.0.11, viridisLite -> 0.4.2, vegan -> 2.6-4, pinfsc50 -> 1.2.0
    # memuse -> 4.2-3, ape -> 5.7-1
    && R --slave -e 'remotes::install_version("vcfR", version = "1.14.0", upgrade=FALSE)' \
    # colorspace -> 2.1-0, munsell -> 0.5.0, labeling -> 0.4.3, farver -> 2.1.1, scales -> 1.2.1
    # isoband -> 0.2.7, gtable -> 0.3.4
    && R --slave -e 'remotes::install_version("ggplot2", version = "3.4.3", upgrade=FALSE)' \
    # depends on ggplot2
    # None
    && R --slave -e 'remotes::install_version("ComplexUpset", version = "1.3.3", upgrade=FALSE)' \
    # depends on ggplot2
    # Matrix -> 1.6-4, prettyunits -> 1.2.0, rematch2 -> 2.1.2, diffobj -> 0.3.5, pkgbuild -> 1.4.2, fs -> 1.6.3, crayon -> 1.5.2
    # rprojroot -> 2.0.3, waldo -> 0.5.1, ps -> 1.7.5, processx -> 3.8.2, praise -> 1.0.0, pkgload -> 1.3.3
    # evaluate -> 0.22, desc -> 1.4.2, callr -> 3.7.3, brio -> 1.1.3, testthat -> 3.1.10 ellipsis -> 0.3.2
    # backports -> 1.4.1, MatrixModels (NA -> 0.5-2, SparseM -> 1.81, numDeriv -> 2016.8-1.1, broom -> 1.0.5
    # RcppEigen -> 0.3.3.9.3, nloptr -> 2.0.3, minqa -> 1.2.6, lme4 -> 1.1-34, quantreg -> 5.97, pbkrtest -> 0.5.2
    # abind -> 1.4-5, carData-> 3.0-5, car -> 3.1-2, corrplot -> 0.92, rstatix -> 0.7.2, polynom -> 1.4-1
    # gridExtra -> 2.3, ggsignif -> 0.6.4, cowplot -> 1.1.1, ggsci -> 3.0.0, ggrepel -> 0.9.3
    && R --slave -e 'remotes::install_version("ggpubr", version = "0.6.0", upgrade=FALSE)' \
    # depends on ggplot2
    # None
    && R --slave -e 'remotes::install_version("cowplot", version = "1.1.1", upgrade=FALSE)' \
    # hms -> 1.1.3, bit-> 4.0.5, progress-> 1.2.2, tzdb -> 0.4.0, bit64 -> 4.0.5, vroom -> 1.6.3, clipr -> 0.8.0)
    && R --slave -e 'remotes::install_version("readr", version = "2.1.4", upgrade=FALSE)' \
    # depends on dplyr, tidyr, ggplot2, readr
    # rappdirs -> 0.3.3, highr -> 0.10, fastmap -> 1.1.1, sass -> 0.4.7, mime -> 0.12, memoise -> 2.0.1
    # cachem -> 1.0.8, base64enc -> 0.1-3, yaml -> 2.3.7, xfun -> 0.40, tinytex -> 0.46, knitr -> 1.44
    # jquerylib -> 0.1.4, htmltools -> 0.5.6, fontawesome -> 0.5.2, bslib -> 0.5.1, systemfonts -> 1.0.4
    # sys -> 3.4.2, askpass -> 1.2.0, uuid -> 1.1-1, openssl -> 2.1.1, curl -> 5.0.2, httr -> 1.4.7
    # gargle -> 1.5.2, rematch -> 2.0.0, xml2 -> 1.3.5, selectr -> 0.4-2, rstudioapi -> 0.15.0, rmarkdown -> 2.25
    # cellranger -> 1.1.0, textshaping -> 0.3.6, timechange -> 0.2.0, forcats -> 1.0.0, ids -> 1.0.1
    # googledrive -> 2.1.1, DBI -> 1.1.3, blob -> 1.2.4, rvest -> 1.0.3, reprex -> 2.0.2, readxl -> 1.4.3
    # ragg -> 1.2.5, modelr -> 0.1.11, lubridate -> 1.9.3, haven -> 2.5.3, googlesheets4 -> 1.1.1, dtplyr -> 1.3.1
    # dbplyr -> 2.3.4, conflicted -> 1.2.0
    && R --slave -e 'remotes::install_version("tidyverse", version = "2.0.0", upgrade=FALSE)' \
    # Using commit of Aug 18, 2023
    # bitops -> 1.0-7, caTools -> 1.18.2, gtools -> 3.9.4, dotCall64 -> 1.0-2, gplots -> 3.1.3
    # maps -> 3.4.1, spam -> 2.9-1, ROCR -> 1.0-11, fields -> 15.2
    && R --slave -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder@1b1d4e2d7f893a3552d9f8f791ab868ee4c782e6", upgrade=FALSE)' \
    # depends on future.apply
    && R --slave -e 'remotes::install_version("hdf5r", version = "1.3.8", upgrade=FALSE)' \
    # progressr -> 0.14.0, sp -> 2.0-0
    && R --slave -e 'remotes::install_version("SeuratObject", version = "4.1.4", upgrade=FALSE)' \
    # depends on future.apply, RColorBrewer, ggplot2, cowplot, SeuratObject
    # sitmo -> 2.0.2, BH -> 1.81.0-1, spatstat.utils -> 3.0-3, tensor -> 1.5, polyclip -> 1.10-6, deldir -> 1.0-9
    # spatstat.geom -> 3.2-5, spatstat.data -> 3.0-1, later -> 1.3.1, promises -> 1.2.1, plyr -> 1.8.8
    # lazyeval -> 0.2.2, commonmark -> 1.9.0, sourcetools -> 0.1.7-1, xtable -> 1.8-4, httpuv -> 1.6.11
    # png -> 0.1-8, here -> 1.0.1, RcppTOML -> 0.2.2, dqrng -> 0.3.1, RcppProgress -> 0.4.2, irlba -> 2.3.5.1
    # RcppAnnoy -> 0.0.21, FNN -> 1.1.3.2, goftest -> 1.2-3, spatstat.sparse -> 3.0-2, spatstat.random -> 3.1-6
    # RcppArmadillo -> 0.12.6.4.0, matrixStats -> 1.0.0, reshape2 -> 1.4.4, crosstalk -> 1.2.0, htmlwidgets -> 1.6.2
    # shiny -> 1.7.5, zoo -> 1.8-12, igraph -> 1.5.1, reticulate -> 1.32.0, uwot -> 0.1.16, spatstat.explore -> 3.2-3
    # sctransform -> 0.4.0, scattermore -> 1.2, Rtsne -> 0.16, RANN -> 2.6.1, plotly -> 4.10.2, pbapply -> 1.7-2
    # miniUI -> 0.1.1.1, lmtest -> 0.9-40, leiden -> 0.4.3, ica -> 1.0-3, ggridges -> 0.5.4, fitdistrplus -> 1.1-11
    && R --slave -e 'remotes::install_version("Seurat", version = "4.4.0", upgrade=FALSE)' \
    # depends on Seurat, prevent from upgrading Seurat to v5.0.0 since this screws with DoubletFinder and more things most likely
    # vipor -> 0.4.5, beeswarm -> 0.4.0, xfun 0.40 -> 0.41, evaluate 0.22 -> 0.23, prismatic -> 1.1.1, skecase -> 0.11.1
    # ggbeeswarm -> 0.7.2, Cairo -> 1.6-1, shape -> 1.4.6, GlobalOptions -> 0.1.2, paletteer -> 1.5.0
    # janitor -> 2.2.0, ggrastr -> 1.0.2, ggprism -> 1.0.4, circlize -> 0.4.15
    && R --slave -e 'remotes::install_version("scCustomize", version = "1.1.3", upgrade=FALSE)'

# BiocManager 1.30.22
# Uses 179 MB, mainly in /usr/local/lib/R/site-library/BiocNeighbors (18 MB), GenomeInfoDbData (11 MB), edgeR (11 MB), scran (12 MB), scuttle (12 MB)
RUN R --slave -e 'install.packages("BiocManager")' \
    # 3.14 is for R version 4.1
    # BiocVersion_3.14.0
    && R --slave -e 'BiocManager::install(version = "3.14")' \
    # rjson_0.2.21, BiocGenerics_0.40.0, S4Vectors_0.32.4, iterators_1.0.14, GetoptLong_1.0.5, clue_0.3-65
    # IRanges_2.28.0, foreach_1.5.2, doParallel_1.0.17, ComplexHeatmap_2.10.0
    && R --slave -e 'BiocManager::install("ComplexHeatmap", version = "3.14")' \
    # RCurl_1.98-1.12, GenomeInfoDbData_1.2.7, zlibbioc_1.40.0, MatrixGenerics_1.6.0, Biobase_2.54.0, GenomeInfoDb_1.30.1
    # XVector_0.34.0, SummarizedExperiment_1.24.0, GenomicRanges_1.46.1, DelayedArray_0.20.0, SingleCellExperiment_1.16.0
    && R --slave -e 'BiocManager::install("SingleCellExperiment", version = "3.14")' \
    # Depends on SingleCellExperiment, 1.14.0
    # formatR_1.14, lambda.r_1.2.4, futile.options_1.0.1, locfit_1.5-9.8, sparseMatrixStats_1.6.0futile.logger_1.4.3, snow_0.4-4
    # RcppHNSW_0.5.0, ScaledMatrix_1.2.0, rsvd_1.0.5, beachmat_2.10.0, edgeR_3.36.0, limma_3.50.3, statmod_1.5.0, DelayedMatrixStats_1.16.0
    # metapod_1.2.0, viridis_0.6.4, BiocParallel_1.28.3, BiocNeighbors_1.12.0, BiocSingular_1.10.0, scran_1.22.1, scater_1.22.0, scuttle_1.4.0
    # bluster_1.4.0, scDblFinder_1.8.0
    && R --slave -e 'BiocManager::install("scDblFinder", version = "3.14")' \
    # pROC_1.18.4, scds_1.10.0
    && R --slave -e 'BiocManager::install("scds", version = "3.14")'

#################################
############## JAVA #############
#################################

# Section build takes 19 seconds and has a size of 0.526 GB.

# adwaita-icon-theme ca-certificates-java fontconfig fontconfig-config
# fonts-dejavu-core gtk-update-icon-cache hicolor-icon-theme
# humanity-icon-theme java-common libasound2 libasound2-data libatk1.0-0
# libatk1.0-data libavahi-client3 libavahi-common-data libavahi-common3
# libbrotli1 libbsd0 libcairo2 libcups2 libdatrie1 libdbus-1-3 libdeflate0
# libdrm-amdgpu1 libdrm-common libdrm-intel1 libdrm-nouveau2 libdrm-radeon1
# libdrm2 libedit2 libelf1 libexpat1 libfontconfig1 libfreetype6 libfribidi0
# libgdk-pixbuf-2.0-0 libgdk-pixbuf2.0-common libgif7 libgl1 libgl1-mesa-dri
# libglapi-mesa libglib2.0-0 libglvnd0 libglx-mesa0 libglx0 libgraphite2-3
# libgtk2.0-0 libgtk2.0-common libharfbuzz0b libicu70 libjbig0 libjpeg-turbo8
# libjpeg8 liblcms2-2 libllvm15 libmd0 libnspr4 libnss3 libpango-1.0-0
# libpangocairo-1.0-0 libpangoft2-1.0-0 libpciaccess0 libpcsclite1
# libpixman-1-0 libpng16-16 libsensors-config libsensors5 libsqlite3-0
# libthai-data libthai0 libtiff5 libwebp7 libx11-6 libx11-data libx11-xcb1
# libxau6 libxcb-dri2-0 libxcb-dri3-0 libxcb-glx0 libxcb-present0
# libxcb-randr0 libxcb-render0 libxcb-shm0 libxcb-sync1 libxcb-xfixes0 libxcb1
# libxcomposite1 libxcursor1 libxdamage1 libxdmcp6 libxext6 libxfixes3 libxi6
# libxinerama1 libxml2 libxrandr2 libxrender1 libxshmfence1 libxtst6
# libxxf86vm1 openjdk-17-jdk openjdk-17-jdk-headless openjdk-17-jre
# openjdk-17-jre-headless shared-mime-info ubuntu-mono ucf x11-common
# Uses 468 MB, mainly in /usr/lib/ (221 MB), /var/ (50 MB), /usr/share/ (22 MB), and /usr/include/ (20 MB)
RUN apt-get install -y --no-install-recommends openjdk-17-jdk

##################################
############## OTHER #############
##################################

# Section build takes 12 minutes and 31 seconds and has a size of 1.931 GB.

# cargo for souporcell
# cargo libssh2-1 libstd-rust-1.66 libstd-rust-dev rustc
# Uses 316 MB
RUN apt-get install -y --no-install-recommends cargo \
    # souporcell, popscle
    # git git-man less libcbor0.8 libcurl3-gnutls libedit2 liberror-perl libfido2-1 libnghttp2-14 librtmp1 libssh-4 libxext6
    # libxmuu1 openssh-client xauth
    && apt-get install -y --no-install-recommends git

## This is only needed if you do not run the Python install section.
## bedtools2
#RUN apt-get install -y --no-install-recommends python3

## This is only needed if you do not run the R install section.
#RUN apt-get install -y --no-install-recommends libz-dev \
#    # bedtools2, samtools, bcftools, htslib
#    # bzip2-doc libbz2-dev
#    && apt-get install -y --no-install-recommends libbz2-dev \
#    # bedtools2, samtools, bcftools, htslib
#    # liblzma-doc
#    && apt-get install -y --no-install-recommends liblzma-dev \
#    # samtools
#    # libncurses-dev
#    && apt-get install -y --no-install-recommends libncurses-dev \
#    # vcftools
#    # libglib2.0-0 libglib2.0-data libicu70 libxml2 pkg-config shared-mime-info xdg-user-dirs
#    && apt-get install -y --no-install-recommends pkg-config \
#    # htslib
#    && apt-get install -y --no-install-recommends libcurl4-openssl-dev \
#    # popscle
#    && apt-get install -y --no-install-recommends cmake \
#    # popscle
#    && apt-get install -y --no-install-recommends libssl-dev \
#    # plink, plink2
#    && apt-get install -y --no-install-recommends unzip

# Requires Python to install.
# Uses 178 MB, mainly in /opt/bedtools2-2.31.0/obj/ (117 MB), /opt/bedtools2-2.31.0/bin/ (41 MB), and /opt/bedtools2-2.31.0/src/ (13 MB)
RUN cd /opt \
    && wget https://github.com/arq5x/bedtools2/archive/refs/tags/v2.31.0.tar.gz \
    && tar -xzf v2.31.0.tar.gz \
    && rm v2.31.0.tar.gz \
    && cd bedtools2-2.31.0 \
      && make \
    && rm -rf /opt/bedtools2-2.31.0/test

# Uses 55 MB, mainly in /opt/samtools-1.18/htslib-1.18/ (32 MB)
RUN cd /opt \
    && wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 \
    && tar xjf samtools-1.18.tar.bz2 \
    && rm samtools-1.18.tar.bz2 \
    && cd samtools-1.18 \
      && ./configure \
      && make \
      && make install \
    && rm -rf /opt/samtools-1.18/test

# Uses 69 MB, mainly in /opt/bcftools-1.18/htslib-1.18/ (32 MB), /opt/bcftools-1.18/bcftools/ (9 MB), /opt/bcftools-1.18/plugins/ (4.2 MB)
RUN cd /opt \
    && wget https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2 \
    && tar xjf bcftools-1.18.tar.bz2 \
    && rm bcftools-1.18.tar.bz2 \
    && cd bcftools-1.18 \
      && ./configure \
      && make \
      && make install \
    && rm -rf /opt/bcftools-1.18/test

# Uses 5.7 MB, mainly in /opt/freebayes-1.3.7/src (3 MB)
RUN cd /opt \
    && wget https://github.com/freebayes/freebayes/archive/refs/tags/v1.3.7.tar.gz \
    && tar -xzf v1.3.7.tar.gz \
    && rm v1.3.7.tar.gz

#
RUN cd /opt \
    && mkdir vartrix-1.1.22 \
    && cd vartrix-1.1.22 \
    && wget https://github.com/10XGenomics/vartrix/releases/download/v1.1.22/vartrix_linux \
    && mv vartrix_linux vartrix \
    && chmod 775 vartrix

# Uses 5.4 MB, mainly in /opt/minimap2-2.26/libminimap2.a/ (1.7 MB)
RUN cd /opt \
    && wget https://github.com/lh3/minimap2/archive/refs/tags/v2.26.tar.gz \
    && tar -xzf v2.26.tar.gz \
    && rm v2.26.tar.gz \
    && cd minimap2-2.26 \
      && make

## Required for souporcell.
## Reducing size using https://github.com/rust-lang/rustup/issues/837.
## Uses 1.2 GB in /opt/.cargo/, mainly in /opt/.cargo/toolchain/1.72.1-x86_64-unknown-linux-gnu/share/ (582 MB), /opt/.cargo/toolchain/1.72.1-x86_64-unknown-linux-gnu/lib/ (501 MB), /opt/.cargo/toolchain/1.72.1-x86_64-unknown-linux-gnu/bin/ (83 MB)
## Uses 0.57 GB /root/.rustup/, mainly in /root/.rustup/toolchains/stable-x86_64-unknown-linux-gnu/lib/ (501 MB)
## Not using this code but instead added apt-get install -y --no-install-recommends cargo above since it uses less space.
#RUN cd /opt \
#    && CARGO_HOME=/opt/.cargo RUSTUP_HOME=/opt/.cargo bash -c 'curl https://sh.rustup.rs -sSf | sh -s -- --default-toolchain=1.72.1 -y' \
#    && . /opt/.cargo/env \
#    && rustup set profile minimal \
#    && rustup default stable-x86_64-unknown-linux-gnu \
#    && rm -rf /root/.rustup/toolchains/stable-x86_64-unknown-linux-gnu/share/doc \
#    && rm -rf /root/.rustup/toolchains/stable-x86_64-unknown-linux-gnu/share/man \
#    && rm -rf /root/.rustup/toolchains/stable-x86_64-unknown-linux-gnu/share/zsh

# Fork from commit 9fb5271ae9f2257ea9a8552dfda3d4b7080be194
# requires samtools, bcftools, htslib, python3, freebayes, vartrix, minimap2, Rust (i.e. cargo)
# Using tips from https://github.com/johnthagen/min-sized-rust to reduce the binary size.
# Uses 529 MB, mainly in /root/.cargo/registry/index/ (384 MB)
RUN cd /opt \
    && git clone --single-branch --branch MVochteloo/gzipped-barcodes-input https://github.com/sc-eQTLgen-consortium/souporcell.git \
    && cd /opt/souporcell/souporcell \
      && echo "\n[profile.release]\nopt-level = 'z'\nlto = true\ncodegen-units = 1\npanic = 'abort'\nstrip = true" >> Cargo.toml \
      && cargo build --release \
    && cd /opt/souporcell/troublet \
      && echo "\n[profile.release]\nopt-level = 'z'\nlto = true\ncodegen-units = 1\npanic = 'abort'\nstrip = true" >> Cargo.toml \
      && cargo build --release \
    && rm -rf /opt/souporcell/.git

# RUN apt-get install -y minimac4=4.1.4
# Uses 14 MB, mainly in /opt/bin/ (11 MB)
RUN cd /opt \
    && wget https://github.com/statgen/Minimac4/releases/download/v4.1.4/minimac4-4.1.4-Linux-x86_64.sh \
    && bash minimac4-4.1.4-Linux-x86_64.sh --prefix=/opt/ --skip-license \
    && rm minimac4-4.1.4-Linux-x86_64.sh

# Uses 84 MB, mainly in /opt/vcftools-0.1.16/src (64 MB)
RUN cd /opt \
    && wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz \
    && tar -xzf vcftools-0.1.16.tar.gz \
    && rm vcftools-0.1.16.tar.gz \
    && cd vcftools-0.1.16 \
      && ./configure \
      && make \
      && make install

# Used by minimap2 and picard? No clue but seems important.
# Uses 11 MB, mainly in /opt/parallel-20230922/src/ (5.9 MB)
RUN cd /opt \
    && wget https://ftp.gnu.org/gnu/parallel/parallel-20230922.tar.bz2 \
    && tar xjf parallel-20230922.tar.bz2 \
    && rm parallel-20230922.tar.bz2 \
    && cd parallel-20230922 \
      && ./configure \
      && make \
      && make install

# Uses 6 MB, mainly in /opt/Eagle_v2.4.1/eagle/ (4.6 MB)
# We can delete the tables since we will supply our own reference map.
RUN cd /opt \
    && wget https://alkesgroup.broadinstitute.org/Eagle/downloads/Eagle_v2.4.1.tar.gz \
    && tar -xzf Eagle_v2.4.1.tar.gz \
    && rm Eagle_v2.4.1.tar.gz \
    && rm -rf /opt/Eagle_v2.4.1/tables/

# Uses 38 MB, mainly in /opt/GenotypeHarmonizer-1.4.27-SNAPSHOT/lib/ (32 MB) and /opt/GenotypeHarmonizer-1.4.27-SNAPSHOT/exampleData/ (5.9 MB)
RUN cd /opt \
    && wget https://github.com/molgenis/systemsgenetics/releases/download/GH_1.4.27/GenotypeHarmonizer-1.4.27-SNAPSHOT-dist.tar.gz \
    && tar -xzf GenotypeHarmonizer-1.4.27-SNAPSHOT-dist.tar.gz  \
    && rm GenotypeHarmonizer-1.4.27-SNAPSHOT-dist.tar.gz \
    && chmod 777 /opt/GenotypeHarmonizer-1.4.27-SNAPSHOT/GenotypeHarmonizer.jar

# Uses 42 MB, mainly in /opt/plink/plink/ (42 MB)
RUN cd /opt \
    && mkdir plink \
    && cd plink \
      && wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip \
      && unzip -q plink_linux_x86_64_20230116.zip \
      && rm plink_linux_x86_64_20230116.zip

# Uses 42 MB, mainly in /opt/plink2/plink2/ (42 MB)
RUN cd /opt \
    && mkdir plink2 \
    && cd plink2 \
      && wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20231005.zip \
      && unzip -q plink2_linux_x86_64_20231005.zip \
      && rm plink2_linux_x86_64_20231005.zip

# Requires openjdk-17-jdk.
# Uses 62 MB, mainly in /opt/picard-3.1.0/build/libs/picard.jar (62 MB)
# I am downloading the jar file directly since it uses less memory that building a new one.
RUN cd /opt \
    && mkdir -p picard-3.1.0/build/libs/ \
    && cd picard-3.1.0/build/libs/ \
    && wget https://github.com/broadinstitute/picard/releases/download/3.1.0/picard.jar \
    && chmod 777 /opt/picard-3.1.0/build/libs/picard.jar

# Switch GCC to version 10.5 since popscle doesn't build with version 11.4:
# stl_algo.h:3455:5: note: 'std::min' declared here
# https://forum.qt.io/topic/59045/build-qt-static-make-error-solved/18
# Note that we do this last so that the rest is build on version 11.4.
# cpp-10 g++-10 gcc-10 gcc-10-base libgcc-10-dev libstdc++-10-dev
# Uses 140 MB, mainly in /usr/lib/ (100 MB), /usr/bin/ (28 MB), and /usr/include/ (10 MB)
RUN apt-get install -y g++-10 gcc-10 \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 80 --slave /usr/bin/g++ g++ /usr/bin/g++-11 --slave /usr/bin/gcov gcov /usr/bin/gcov-11 \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 100 --slave /usr/bin/g++ g++ /usr/bin/g++-10 --slave /usr/bin/gcov gcov /usr/bin/gcov-10 \
    && update-alternatives --config gcc

# Required for popscle.
# Uses 109 MB, mainly in /opt/htslib-1.18/ (71 MB), /usr/local/lib/ (19 MB), /usr/local/bin/ (18 MB)
RUN cd /opt \
    && wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 \
    && tar xjf htslib-1.18.tar.bz2 \
    && rm htslib-1.18.tar.bz2 \
    && cd htslib-1.18 \
      && ./configure \
      && make \
      && make install \
    && rm -rf /opt/htslib-1.18/test

# Commit of May 5, 2021
# Source: https://github.com/statgen/popscle/issues/21
# Uses 7 MB, mainly in /opt/popscle/tutorials/ (2.9 MB), /opt/popscle/build/ (2.3 MB)
RUN cd /opt \
    && git clone https://github.com/statgen/popscle.git \
    && cd popscle \
      && git checkout da70fc78da385ef049e0e890342acfd62842cae0 \
      && mkdir build \
      && cd build \
        && cmake ..  \
        && make \
    && rm -rf /opt/popscle/.git

####################################
################ CLEAN #############
####################################

RUN apt-get clean \
    && apt-get autoremove -y
