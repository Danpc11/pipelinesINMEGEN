FROM ubuntu:20.04

LABEL "container.pipelines_inmegen"="VCgermiline1"
LABEL version="0.1"

# Building takes several hours! Be sure to symlink /var/lib/docker to a larger disk.
# Build the image: docker build -t fungi .
# Run testcases:   docker run --rm -ti fungi tests/run-tests.sh


ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    automake \
    build-essential \
    bzip2 \
    ca-certificates \
    cmake \
    curl \
    default-jdk \
    gawk \
    gcc \
    git \
    gmap \
    gzip \
    g++ \
    hmmer \
    libbz2-dev \
    libc6-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libdb-dev \
    libipc-run-perl \
    liblzma-dev \
    libncurses5-dev \
    libtbb2 \
    libtbb-dev \
    libxml2-dev \
    liblapack-dev \
    libopenblas-dev \
    make \
    parallel \
    pigz \
    pkg-config \
    python3-dev \
    python-dev \
    python-numpy \
    python-setuptools \
    rna-star \
    r-base \
    r-base-dev \
    tar \
    unzip \
    wget \
    zip \
    zlib1g \
    zlib1g-dev \
  && apt-get clean

RUN curl https://bootstrap.pypa.io/pip/3.6/get-pip.py | python3 && \
    pip3 install \
      Cython \
      numpy \
      biopython \
      xlrd \
      openpyxl \
      pandas \
      igv-reports==0.2.1 \
    && \
    rm -rdf /root/.cache/pip

RUN curl -jL https://github.com/pypa/pip/archive/20.2.2.tar.gz | tar xvz && \
    cd pip-20.2.2 && python2 setup.py install && \
        pip2 install \
        numpy \
        biopython==1.76 \
        xlrd \
        openpyxl \
        pandas \
        scikit-learn \
    && \
    rm -rdf /root/.cache/pi

ENV PERL_MM_USE_DEFAULT=1
RUN cpan App::cpanminus
RUN cpanm DB_File \
    URI::Escape \
    Set::IntervalTree \
    Carp::Assert \
    JSON::XS \
    PerlIO::gzip \
    Parallel::ForkManager \
    Bio::Seq

RUN Rscript -e 'install.packages("circlize", repos="http://cran.r-project.org"); \
        install.packages("BiocManager"); \
        install.packages("ranger"); \
        install.packages("tidyverse"); \
        install.packages("vroom"); \
        install.packages("gsalib"); \
        install.packages("ggplot2"); \
        install.packages("reshape"); \
        install.packages("gplots"); \
        BiocManager::install(c("GenomicRanges", "GenomicAlignments", "rlang"))'

WORKDIR /usr/local/lib/

####### STAR ######
RUN git clone https://github.com/alexdobin/STAR.git 
RUN ln -s /usr/local/lib/STAR/bin/Linux_x86_64_static/STAR /usr/local/bin/STAR

###### Samtools & friends #####
RUN git clone --recurse-submodules https://github.com/samtools/htslib.git
RUN git clone https://github.com/samtools/bcftools.git && \
    cd bcftools && \
    make
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
    tar xjvf samtools-1.16.1.tar.bz2 && \
    cd samtools-1.16.1 && \
    ./configure --prefix=/usr/local/ && \
    make && \
    make install

###### BWA #####
RUN git clone https://github.com/lh3/bwa.git
    cd bwa; make

###### otros ####

RUN apt-get install -y --no-install-recommends gffread tzdata libvcflib-tools libvcflib-dev
ENV TZ="America/Mexico_City"

CMD ["bash", "-l"]
