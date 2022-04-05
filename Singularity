Bootstrap: docker
From: ubuntu:20.04
%help
    Container for TrEmolo v2.0
    https://github.com/DrosophilaGenomeEvolution/TrEMOLO
    Includes
        Blast 2.2+
        Bedtools2
        RaGOO v1.1
        Assemblytics
        Snakemake 5.5.2+
        Minimap2 2.16+
        Samtools 1.9
        Sniffles 1.0.10+
        SVIM 1.4.2+
        libfontconfig1-dev
        Python libs
            Biopython
            Pandas
            Numpy
            pylab
            intervaltree
        R libs
            knitr 1.38
            rmarkdown 2.13
            bookdown 0.25
            viridis 0.6.2
            viridisLite 0.4.0
            rjson 0.2.20
            ggthemes 4.2.4
            forcats 0.5.1
            reshape2 1.4.4
            dplyr 1.0.8
            kableExtra 1.3.4
            extrafont 0.17
            ggplot2 3.3.5
            RColorBrewer 1.1-2
        Perl v5.26.2

%labels
    VERSION "TrEMOLO v2.0"
    Maintainer Francois Sabot <francois.sabot@ird.fr>
    March, 2021

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

    # apt dependencies
    apt update
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
        minimap2 \
        snakemake \
        assemblytics \
        r-base \
        perl \
        pandoc-citeproc \
        libfontconfig1-dev \
        libxml2-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        curl


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
    #rmarkdown_2.13     knitr_1.36         bookdown_0.25      viridis_0.6.2     
    #viridisLite_0.4.0  rjson_0.2.20       ggthemes_4.2.4     forcats_0.5.1     
    #reshape2_1.4.4     dplyr_1.0.8        kableExtra_1.3.4   extrafont_0.17    
    #ggplot2_3.3.5      RColorBrewer_1.1-2

    #R --slave -e 'devtools::install_github("yihui/knitr@v1.36")'
    R --slave -e 'require(devtools); install_version("knitr", version = "1.38")'
    R --slave -e 'require(devtools); install_version("rmarkdown", version = "2.13")'
    R --slave -e 'require(devtools); install_version("bookdown", version = "0.25")'
    R --slave -e 'require(devtools); install_version("viridis", version = "0.6.2")'
    R --slave -e 'require(devtools); install_version("viridisLite", version = "0.4.0")'
    R --slave -e 'require(devtools); install_version("rjson", version = "0.2.20")'
    R --slave -e 'require(devtools); install_version("ggthemes", version = "4.2.4")'
    R --slave -e 'require(devtools); install_version("forcats", version = "0.5.1")'
    R --slave -e 'require(devtools); install_version("reshape2", version = "1.4.4")'
    R --slave -e 'require(devtools); install_version("dplyr", version = "1.0.8")'
    R --slave -e 'require(devtools); install_version("kableExtra", version = "1.3.4")'
    R --slave -e 'require(devtools); install_version("extrafont", version = "0.17")'
    R --slave -e 'require(devtools); install_version("ggplot2", version = "3.3.5")'
    R --slave -e 'require(devtools); install_version("RColorBrewer", version = "1.1-2")'    
    

    #samtools
    cd /usr/bin
    wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
    tar -vxjf htslib-1.9.tar.bz2
    cd htslib-1.9
    make

    cd ..
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
    tar -vxjf samtools-1.9.tar.bz2
    cd samtools-1.9
    make

    cd ..
    wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
    tar -vxjf bcftools-1.9.tar.bz2
    cd bcftools-1.9
    make

    #Get vcfutils.pl
    cd ../
    git clone https://github.com/lh3/samtools.git
    cd samtools/bcftools/
    make

    cd
    
    #Python libs
    python3 -m pip install biopython pandas numpy==1.21.2 matplotlib svim==1.4.2 intervaltree scipy pysam

    # build variables
    export TOOLDIR=/opt/tools

    #Preparing Directories
    mkdir -p $TOOLDIR

    #installing SamTools 1.9
    cd $TOOLDIR
    mkdir samtools1.9
    cd samtools1.9
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
    tar -xf samtools-1.9.tar.bz2
    cd samtools-1.9
    ./configure
    make all all-htslib
    make install install-htslib

    #installing Sniffles 1.0.12+
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


    #installing Flye 2.8+
    # cd $TOOLDIR
    # git clone https://github.com/fenderglass/Flye
    # cd Flye
    # python3 setup.py install

    #install RaGOO
    cd $TOOLDIR
    git clone https://github.com/malonge/RaGOO.git
    cd RaGOO
    python3 setup.py install

    #install WTDBG2
    # cd $TOOLDIR
    # git clone https://github.com/ruanjue/wtdbg2
    # cd wtdbg2
    # make

    #install TrEMOLO
    cd $TOOLDIR
    git clone https://github.com/DrosophilaGenomeEvolution/TrEMOLO.git

    #Must be compiling at the end
    R --save -e 'install.packages("stringi")'
    R --save -e 'install.packages("stringr")'

%environment
    export LC_ALL=C
    #export PATH=$TOOLDIR/wtdbg2/:$TOOLDIR/RaGOO/:$TOOLDIR/Flye/:$TOOLDIR/TrEMOLO/:$PATH
    export PATH=$TOOLDIR/RaGOO/:$TOOLDIR/TrEMOLO/:$PATH
    export PATH="$PATH:/usr/bin/bcftools-1.9"
    export PATH="$PATH:/usr/bin/samtools-1.9"
    export PATH="$PATH:/usr/bin/htslib-1.9"
    export PATH="$PATH:/usr/bin/samtools/bcftools/"

%runscript
    exec "$@"
