# install python dependencies
FROM ubuntu:20.04 AS cat-python

RUN apt update && apt install -y --no-install-recommends \
  gcc \
  python3-dev \
  python3-pip \
  wget

COPY ./setup.py /

RUN wget https://ont-research.s3-eu-west-1.amazonaws.com/parasail-1.1.17-py2.py3-none-manylinux1_x86_64.whl
RUN pip3 install parasail-1.1.17-py2.py3-none-manylinux1_x86_64.whl

RUN mkdir cat tools \
  && python3 setup.py egg_info \
  && pip3 install -r cat.egg-info/requires.txt

COPY ./ /cat

RUN cd /cat \
  && sed -i'' "s#'augustus_cfgs/#'/opt/augustus/config/extrinsic/#g" cat/__init__.py \
  && python3 setup.py install

########################################

FROM curlimages/curl:7.70.0 AS cat-binaries

USER root

WORKDIR /binaries

# Need >= v395 for clusterGenes -minOverlappingBases option
RUN curl -LO http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/{axtChain,bamToPsl,bedSort,bedToBigBed,chainMergeSort,clusterGenes,faToTwoBit,genePredToBed,genePredToFakePsl,genePredToGtf,gff3ToGenePred,gtfToGenePred,pslCDnaFilter,pslCheck,pslMap,pslMapPostChain,pslPosTarget,pslRecalcMatch,pslToBigPsl,transMapPslToGenePred,wigToBigWig,hgLoadChain} \
  && chmod a+x /binaries/*

RUN set -o pipefail && curl -L https://github.com/biod/sambamba/releases/download/v0.7.1/sambamba-0.7.1-linux-static.gz \
    | gzip -d > /binaries/sambamba && chmod a+x /binaries/sambamba

# CAT v2.1.0 needs more recent hal2fasta supporting the --onlySequenceNames option
#RUN set -o pipefail && curl -L https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v1.0.0/cactus-bin-v1.0.0.tar.gz \
#    | tar -C /tmp -xzf - \
#        cactus-bin-v1.0.0/bin/hal2maf \
#        cactus-bin-v1.0.0/bin/hal2fasta \
#        cactus-bin-v1.0.0/bin/halLiftover \
#        cactus-bin-v1.0.0/bin/halStats \
#  && mv /tmp/cactus-bin-v1.0.0/bin/* /binaries && chmod a+x /binaries/hal*

########################################

FROM ubuntu:20.04 AS cat-augustus

# Install required packages
RUN apt update && DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends \
  autoconf \
  build-essential \
  ca-certificates \
  curl \
  libbamtools-dev \
  libboost-iostreams-dev \
  libgsl-dev \
  libhts-dev \
  liblpsolve55-dev \
  libsqlite3-dev \
  libsuitesparse-dev \
  zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*

# 2020-07-03 snapshot + bam2wig build simplification
# https://github.com/Gaius-Augustus/Augustus/pull/153
RUN mkdir /src && cd /src \
  && curl -L https://github.com/harvardinformatics/Augustus/archive/08b7d320cbee586ebfbee410aeae90d81ce03f1e.tar.gz \
    | tar --strip-components=1 -xzf - \
  && make INSTALLDIR=/opt/augustus MYSQL=false HTSLIBS='-lhts' \
  && make -j install \
  && mv /opt/augustus-* /opt/augustus \
  && rm -rf /src

########################################
# https://github.com/Ensembl/WiggleTools/blob/597d84/Dockerfile

FROM ubuntu:20.04 AS cat-wiggletools

RUN apt update && apt install -y --no-install-recommends \
  ca-certificates \
  libgsl-dev \
  libhts-dev \
  libbigwig-dev \
  libcurl4-openssl-dev \
  gcc \
  python \
  make

WORKDIR /build

# 2020-06-02 snapshot
ADD https://github.com/Ensembl/WiggleTools/archive/c1daac89e3775bc8f96376fc1ed7f7e645ce168c.tar.gz wiggletools.tar.gz

RUN tar --strip-components=1 -xzf wiggletools.tar.gz \
  && make LIBS='-lwiggletools -lBigWig -lcurl -lhts -lgsl  -lgslcblas -lz -lpthread -lm -llzma'

########################################

FROM ubuntu:20.04 AS cat-hal

RUN apt update && apt install -y \
  libhdf5-dev \
  g++ \
  make \
  zlib1g-dev

WORKDIR /sonLib
# 2020-06-16 snapshot
ADD https://github.com/ComparativeGenomicsToolkit/sonLib/archive/ea0b939828ba24d998a7c1aa407ff5a016912f56.tar.gz sonLib.tar.gz
RUN tar --strip-components=1 -xzf sonLib.tar.gz
RUN make -j

WORKDIR /hal
# 2020-07-08 snapshot
ADD https://github.com/ComparativeGenomicsToolkit/hal/archive/cb7c044731271ec41640db71f5694af53a0ead57.tar.gz hal.tar.gz
RUN tar --strip-components=1 -xzf hal.tar.gz \
  && make -j \
  && mkdir /binaries \
  && mv bin/hal2fasta bin/hal2maf bin/halStats bin/halLiftover /binaries \
  && strip /binaries/* \
  && rm -rf /hal

########################################

FROM ubuntu:20.04 AS final

RUN apt update && DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends \
  bamtools \
  bedtools \
  exonerate \
  libbamtools2.5.1 \
  libbigwig0 \
  libboost-iostreams1.71.0 \
  libcolamd2 \
  libcurl4 \
  libgsl23 \
  libhts3 \
  libsqlite3-0 \
  libsz2 \
  libsuitesparseconfig5 \
  python3-pip \
  samtools \
  && rm -rf /var/lib/apt/lists/*

COPY --from=cat-python /usr/local /usr/local
COPY --from=cat-binaries /binaries /usr/local/bin
COPY --from=cat-hal /binaries /usr/local/bin
COPY --from=cat-wiggletools /build/bin/wiggletools /usr/local/bin
COPY --from=cat-augustus /opt/augustus /opt/augustus

# (2020-06-04) augustus_cfgs/log_reg_parameters_default.cfg identical to Augustus config/cgp/log_reg_parameters_default.cfg
COPY ./augustus_cfgs/*extrinsic*.cfg /opt/augustus/config/extrinsic/

# luigi looks for luigi.cfg in /etc/luigi/luigi.cfg by default
COPY ./logging.cfg ./luigi.cfg /etc/luigi/

# but need to tell luigi to look for logging.cfg at /etc/luigi/logging.cfg
RUN sed -i'' '/logging_conf_file/s#.*#logging_conf_file=/etc/luigi/logging.cfg#' /etc/luigi/luigi.cfg

ENV PATH=/opt/augustus/bin:/opt/augustus/scripts:${PATH}
