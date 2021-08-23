FROM ubuntu:18.04 AS builder
ARG AUGUSTUS_COMMIT=36ae43d
RUN apt-get update
RUN apt-get install -y build-essential libssl-dev libncurses5-dev libcurl4-openssl-dev liblzma-dev libbz2-dev \
libboost-all-dev sqlite3 libsqlite3-0 libsqlite3-dev libgsl0-dev lp-solve liblpsolve55-dev libbamtools-dev wget git

# htslib
RUN git clone --recursive git://github.com/samtools/htslib.git
RUN cd htslib && make install

# bcftools
RUN git clone git://github.com/samtools/bcftools.git
RUN cd bcftools && make

# samtools
RUN git clone git://github.com/samtools/samtools
RUN cd samtools && make && make install

# MOVE Directories INTO $HOME/tool
RUN mkdir /root/tools
RUN mv samtools /root/tools
RUN mv htslib /root/tools
RUN mv bcftools /root/tools

# Augustus
RUN git clone https://github.com/Gaius-Augustus/Augustus augustus
RUN cd augustus && git reset --hard ${AUGUSTUS_COMMIT}
RUN echo 'COMPGENEPRED = true' >> augustus/common.mk
RUN echo 'SQLITE = true' >> augustus/common.mk
RUN cd augustus/auxprogs/homGeneMapping/src && sed 's/# BOOST = true/BOOST = true/g' -i Makefile && sed 's/# SQLITE = true/SQLITE = true/g' -i Makefile
RUN cd augustus && make

# HDF5
RUN wget -q http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz
RUN tar xzf hdf5-1.10.1.tar.gz
RUN cd hdf5-1.10.1 && ./configure --enable-cxx --prefix=/usr
RUN cd hdf5-1.10.1 && make && make install

# sonLib
RUN git clone git://github.com/ComparativeGenomicsToolkit/sonLib.git

# HAL
RUN git clone git://github.com/ComparativeGenomicsToolkit/hal.git
RUN cd sonLib && make
RUN cd hal && make

# LibBigWig
RUN git clone https://github.com/dpryan79/libBigWig.git
RUN cd libBigWig && make install

# WiggleTools
RUN git clone https://github.com/dahlo/WiggleTools
# Their makefile now hardcodes /bin/cc as compiler :(
RUN ln -s /usr/bin/cc /bin/cc
RUN cd WiggleTools && make

# sambamba
RUN wget -q https://github.com/biod/sambamba/releases/download/v0.6.7/sambamba_v0.6.7_linux.tar.bz2
RUN tar xvjf sambamba_v0.6.7_linux.tar.bz2

# Slimmer final Docker image

FROM ubuntu:18.04
RUN apt-get update
RUN apt-get install -y wget bedtools bamtools samtools sqlite3 libgsl0-dev libcolamd2 software-properties-common libcurl4-openssl-dev exonerate
RUN add-apt-repository -y ppa:deadsnakes/ppa
RUN apt-get install -y python3.7 python3-pip
# Kent
RUN for i in wigToBigWig faToTwoBit gff3ToGenePred genePredToBed genePredToFakePsl bamToPsl transMapPslToGenePred \
pslPosTarget axtChain chainMergeSort pslMap pslRecalcMatch pslMapPostChain gtfToGenePred genePredToGtf bedtools \
pslCheck pslCDnaFilter clusterGenes pslToBigPsl bedSort bedToBigBed hgLoadChain; do \
wget -q http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/$i -O /bin/$i ; chmod +x /bin/$i ; done

COPY --from=builder /hal/bin/* /bin/
COPY --from=builder /sambamba /bin/
COPY --from=builder /augustus/bin/* /bin/
COPY --from=builder /augustus/scripts/* /bin/
COPY --from=builder /WiggleTools/bin/* /bin/

RUN mkdir -p /augustus
COPY --from=builder /augustus/config /augustus/config

# Python deps
RUN pip3 install bd2k-python-lib toil[all]==5.0 pyfasta numpy matplotlib pandas==1.0

# make Python 3 primary python
RUN rm /usr/bin/python
RUN ln -s /usr/bin/python3.7 /usr/bin/python

ENV AUGUSTUS_CONFIG_PATH=/augustus/config/
