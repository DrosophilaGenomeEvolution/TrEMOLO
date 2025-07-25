Bootstrap: docker
From: ubuntu:20.04
%help
    Container for TrEMOLO v2.4
    https://github.com/DrosophilaGenomeEvolution/TrEMOLO
    Includes
        Blast 2.2+
        Bedtools2
        RaGOO v1.1
        Assemblytics
        Snakemake 7.32.4
        Minimap2 2.28+
        Samtools 1.15.1
        Samtools 1.20
        Sniffles 1.0.12b
        SVIM 1.4.2+
        libfontconfig1-dev
        Liftoff
        Python libs
            Biopython 1.79
            Pandas 1.5.3
            Numpy 1.21.2
            matplotlib 3.5.1
            pysam 0.20.0
            intervaltree 2.1.0
            scipy 1.10.1
        R libs
            knitr 1.38
            rmarkdown 2.26
            bookdown 0.38
            viridis 0.6.2
            viridisLite 0.4.0
            rjson 0.2.20
            ggthemes 4.2.4
            forcats 0.5.1
            reshape2 1.4.4
            dplyr 1.0.8
            kableExtra 1.3.4
            extrafont 0.17
            ggplot2 3.4.2
            RColorBrewer 1.1-2
        Perl v5.26.2

%labels
    VERSION "TrEMOLO v2.4"
    Maintainer Francois Sabot <francois.sabot@ird.fr>
    Oct, 2023

%post
    # faster apt downloads
    export DEBIAN_FRONTEND=noninteractive
    export LC_ALL=C
    (
        . /etc/os-release
        cat << _EOF_ > mirror.txt
deb mirror://mirrors.ubuntu.com/mirrors.txt ${UBUNTU_CODENAME} main restricted universe multiverse
deb mirror://mirrors.ubuntu.com/mirrors.txt ${UBUNTU_CODENAME}-updates main restricted universe multiverse
deb mirror://mirrors.ubuntu.com/mirrors.txt ${UBUNTU_CODENAME}-backports main restricted universe multiverse
deb mirror://mirrors.ubuntu.com/mirrors.txt ${UBUNTU_CODENAME}-security main restricted universe multiverse

_EOF_
        mv /etc/apt/sources.list /etc/apt/sources.list.bak
        cat mirror.txt /etc/apt/sources.list.bak > /etc/apt/sources.list
    )

    set -e
    # apt dependencies
    apt update
    apt install -y ca-certificates

    apt-get install -y \
        r-base r-base-dev \
        libcurl4-openssl-dev libssl-dev libxml2-dev \
        libfontconfig1-dev libcairo2-dev libpango1.0-dev \
        libpng-dev libjpeg-dev libtiff5-dev libicu-dev \
        libfreetype6-dev libharfbuzz-dev libfribidi-dev 
    
    apt install -y \
    	apt-utils \
        autoconf \
        automake \
        cmake \
        gcc \
        build-essential \
        software-properties-common \
        tar \
        unzip \
        wget \
        zlib1g-dev \
        sudo \
        git-core \
        locales \
        python3-pip \
        ncbi-blast+ \
        bedtools \
        assemblytics \
        perl \
        pandoc-citeproc \
        curl \
        bc \
        texlive-full \
        snakemake
    
    mkdir -p /usr/local/lib/R/site-library && chmod 775 /usr/local/lib/R/site-library

    R -e "install.packages('remotes', repos='https://cloud.r-project.org')"

    R -e "remotes::install_version('knitr', version='1.38', repos='https://cloud.r-project.org')"
    R -e "remotes::install_version('rmarkdown', version='2.26', repos='https://cloud.r-project.org')"
    R -e "remotes::install_version('bookdown', version='0.38', repos='https://cloud.r-project.org')"
    
    R -e "remotes::install_version('gtable',    version='0.3.0', repos='https://cloud.r-project.org')"
    R -e "remotes::install_version('gridExtra', version='2.3',   repos='https://cloud.r-project.org')"
    R -e "remotes::install_version('scales',    version='1.1.1', repos='https://cloud.r-project.org')"
    R -e "remotes::install_version('purrr',     version='0.3.4', repos='https://cloud.r-project.org')"

    R -e "remotes::install_version('ggplot2', version='3.3.6', repos='https://cloud.r-project.org')"
    
    R -e "remotes::install_version('viridis',    version='0.6.2', repos='https://cloud.r-project.org')"
    R -e "remotes::install_version('viridisLite',version='0.4.0', repos='https://cloud.r-project.org')"
    R -e "remotes::install_version('rjson',      version='0.2.20',repos='https://cloud.r-project.org')"
    R -e "remotes::install_version('ggthemes',   version='4.2.4', repos='https://cloud.r-project.org')"
    R -e "remotes::install_version('forcats',    version='0.5.1', repos='https://cloud.r-project.org')"
    R -e "remotes::install_version('reshape2',   version='1.4.4', repos='https://cloud.r-project.org')"
    R -e "remotes::install_version('dplyr',      version='1.0.8', repos='https://cloud.r-project.org')"
    R -e "remotes::install_version('extrafont',  version='0.17',  repos='https://cloud.r-project.org')"
    # R -e "remotes::install_version('ggplot2',    version='3.4.2', repos='https://cloud.r-project.org')"
    R -e "remotes::install_version('RColorBrewer', version='1.1-2', repos='https://cloud.r-project.org')"

    R -e "remotes::install_version('scales', version='1.1.1', repos='https://cloud.r-project.org')"
    R -e "remotes::install_version('systemfonts', version='1.0.4', repos='https://cloud.r-project.org')"
    R -e "remotes::install_version('textshaping', version='0.3.6', repos='https://cloud.r-project.org')"
    # R -e "remotes::install_version('svglite', version='2.1.0', repos='https://cloud.r-project.org')"
    R -e "remotes::install_version('cpp11', version='0.4.6', repos='https://cloud.r-project.org')"

    # Installer les dépendances nécessaires en précisant les versions compatibles
    R -e "remotes::install_version('systemfonts', version='1.0.4', repos='https://cloud.r-project.org')"
    R -e "remotes::install_version('textshaping', version='0.3.6', repos='https://cloud.r-project.org')"
    R -e "remotes::install_version('svglite', version='2.1.0', repos='https://cloud.r-project.org')"

    # Installer finalement kableExtra
    R -e "remotes::install_version('kableExtra', version='1.3.4', repos='https://cloud.r-project.org', dependencies=NA)"

    R -e "remotes::install_version('kableExtra', version='1.3.4', repos='https://cloud.r-project.org', dependencies=NA)" 2>&1 | tee -a /tmp/kableExtra_install.log

    # Test de vérification : doit retourner TRUE, sinon échec
    R -e "if (!('kableExtra' %in% rownames(installed.packages()))) { cat('kableExtra non installé\n') }"
    


    # R dependencies
    # R --slave -e 'install.packages("knitr")'
    # R --slave -e 'install.packages("rmarkdown")'
    # R --slave -e 'install.packages("bookdown")'
    # R --slave -e 'install.packages("viridis")'
    # R --slave -e 'install.packages("viridisLite")'
    # R --slave -e 'install.packages("rjson")'
    # R --slave -e 'install.packages("forcats")'
    # R --slave -e 'install.packages("ggthemes")'
    # R --slave -e 'install.packages("reshape2")'
    # R --slave -e 'install.packages("dplyr")'
    # R --slave -e 'install.packages("kableExtra")'
    # R --slave -e 'install.packages("extrafont")'
    # R --slave -e 'install.packages("ggplot2")'
    # R --slave -e 'install.packages("RColorBrewer")'

    #Good Version
    #rmarkdown_2.22     knitr_1.36         bookdown_0.25      viridis_0.6.2     
    #viridisLite_0.4.0  rjson_0.2.20       ggthemes_4.2.4     forcats_0.5.1     
    #reshape2_1.4.4     dplyr_1.0.8        kableExtra_1.3.4   extrafont_0.17    
    #ggplot2_3.4.2      RColorBrewer_1.1-2


    pip install pulp==2.6.0
    # # pip install snakemake==7.32.4

    ##minimap2-2.28
    cd /usr/bin
    #git clone https://github.com/lh3/minimap2
    wget https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28.tar.bz2
    tar -vxjf minimap2-2.28.tar.bz2
    cd minimap2-2.28 && make


    ##samtools1.15.1
    cd /usr/bin
    wget https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2
    tar -vxjf htslib-1.15.1.tar.bz2
    cd htslib-1.15.1
    make

    cd ..
    wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2
    tar -vxjf samtools-1.15.1.tar.bz2
    cd samtools-1.15.1
    make

    cd ..
    wget https://github.com/samtools/bcftools/releases/download/1.15.1/bcftools-1.15.1.tar.bz2
    tar -vxjf bcftools-1.15.1.tar.bz2
    cd bcftools-1.15.1
    make

    ##samtools 1.20
    cd /usr/bin
    wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2
    tar -vxjf htslib-1.20.tar.bz2
    cd htslib-1.20
    make

    cd ..
    wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2
    tar -vxjf samtools-1.20.tar.bz2
    cd samtools-1.20
    make

    cd ..
    wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2
    tar -vxjf bcftools-1.20.tar.bz2
    cd bcftools-1.20
    make

    #Get vcfutils.pl
    cd ../
    git clone https://github.com/lh3/samtools.git
    cd samtools/bcftools/
    make

    cd

    #Python libs
    python3 -m pip install biopython==1.79 pandas==1.5.3 matplotlib==3.5.1 svim==1.4.2 intervaltree==2.1.0 scipy==1.10.1 numpy==1.24.1 pysam==0.23.3

    # python3 -m pip install Liftoff
    # build variables
    export TOOLDIR=/opt/tools

    #Preparing Directories
    mkdir -p $TOOLDIR

    #liftoff
    # mkdir -p env/python
    # python3 -m venv env/python/liftoff_env
    # source env/python/liftoff_env/bin/activate

    # python3 -m pip install Liftoff
    # python3 -m pip install numpy==1.21 pysam==0.16.0.1 biopython==1.76 requests==2.20.1

    # source env/python/liftoff_env/bin/deactivate
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda

    # Ajouter conda au PATH
    /opt/conda/bin/conda init bash

    # Créer un environnement conda pour Liftoff
    /opt/conda/bin/conda create --name liftoff_env python=3.8.10 numpy=1.20 -y

    # Activer l'environnement et installer Liftoff
    /bin/bash -c "source /opt/conda/bin/activate liftoff_env"
    /opt/conda/envs/liftoff_env/bin/pip install liftoff
    /opt/conda/envs/liftoff_env/bin/pip install numpy==1.21 pysam==0.16.0.1 biopython==1.76 requests==2.20.1

    /bin/bash -c "source /opt/conda/bin/deactivate"


    #installing Sniffles 1.0.12b
    cd $TOOLDIR
    wget https://github.com/fritzsedlazeck/Sniffles/archive/refs/tags/v1.0.12b.tar.gz -O Sniffles.tar.gz
    tar xzvf Sniffles.tar.gz
    cd Sniffles-1.0.12b/
    mkdir -p build/
    cd build/
    cmake ..
    make
    cd ../bin/sniffles*
    ln -s $PWD/sniffles /usr/bin/sniffles


    #install RaGOO
    cd $TOOLDIR
    git clone https://github.com/malonge/RaGOO.git
    cd RaGOO
    python3 setup.py install

    #For ragoo running
    ln -s /usr/bin/python3 /usr/bin/python

    #install TrEMOLO
    cd $TOOLDIR
    git clone https://github.com/DrosophilaGenomeEvolution/TrEMOLO.git

    # Installation de nvm (Node Version Manager)
    mkdir -p /opt/nvm/
    export NVM_DIR=/opt/nvm
    curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.39.1/install.sh | bash
    
    . "$NVM_DIR/nvm.sh"
    
    # Installation Node.js v18.8.0
    nvm install 18.8.0
    nvm alias default 18.8.0
    nvm use default

    # add nvm and node all usr
    echo 'export NVM_DIR="/opt/nvm"' >> /etc/profile.d/nvm.sh
    echo '[ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh"' >> /etc/profile.d/nvm.sh
    echo '[ -s "$NVM_DIR/bash_completion" ] && \. "$NVM_DIR/bash_completion"' >> /etc/profile.d/nvm.sh

    #Must be compiling at the end
    R --save -e 'install.packages("stringi")'
    R --save -e 'install.packages("stringr")'

    #Force R path container
    mkdir -p /opt/init-file/
    echo '.libPaths(c("/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library"))' > /opt/init-file/.Rprofile

%environment
    export LC_ALL=C
    export TOOLDIR=/opt/tools
    export PATH=$TOOLDIR/RaGOO/:$TOOLDIR/TrEMOLO/:$PATH
    export PATH="$PATH:/usr/bin/bcftools-1.20"
    export PATH="$PATH:/usr/bin/samtools-1.20"
    export PATH="$PATH:/usr/bin/htslib-1.20"
    export PATH="$PATH:/usr/bin/minimap2-2.28"

    export SAMTOOLS_1_9="/usr/bin/samtools-1.20"
    export SAMTOOLS_1_15_1="/usr/bin/samtools-1.15.1"
    # export PATH="$PATH:/usr/bin/bcftools-1.15.1"
    # export PATH="$PATH:/usr/bin/samtools-1.15.1"
    # export PATH="$PATH:/usr/bin/htslib-1.15.1"
    export PATH="$PATH:/usr/bin/samtools/bcftools/"

    export PATH=/opt/conda/bin:$PATH
    export PATH=/usr/bin:$PATH

    export PATH="/opt/nvm/versions/node/v18.8.0/bin:$PATH"
    # Fix force R path container
    export R_PROFILE=/opt/init-file/.Rprofile
    export R_LIBS="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"
    export R_LIBS_SITE="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"

%runscript
    exec "$@"
