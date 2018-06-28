FROM ubuntu:18.04
RUN apt-get update
RUN apt-get install wget -y
RUN apt-get install git -y

# Kent
RUN for i in faToTwoBit gff3ToGenePred genePredToBed genePredToFakePsl bamToPsl transMapPslToGenePred pslPosTarget axtChain chainMergeSort pslMap pslRecalcMatch pslMapPostChain gtfToGenePred genePredToGtf bedtools blat pslCheck pslCDnaFilter clusterGenes pslToBigPsl bedSort bedToBigBed ; do wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/$i -O /bin/$i ; chmod +x /bin/$i ; done

# bedtools
RUN apt-get install bedtools -y

# Augustus dependencies
RUN apt-get install libboost-all-dev sqlite3 libsqlite3-0 libsqlite3-dev libgsl0-dev lp-solve liblpsolve55-dev bamtools libbamtools-dev -y
# zlib
RUN apt-get install liblzma-dev -y
RUN apt-get install libbz2-dev -y
# libcurl
RUN apt-get install libcurl4-openssl-dev -y
# htslib
RUN git clone git://github.com/samtools/htslib.git
RUN cd htslib && make install
# bcftools
RUN git clone git://github.com/samtools/bcftools.git 
RUN cd bcftools && make
# cureses
RUN apt-get install libncurses5-dev -y
# samtools
RUN git clone git://github.com/samtools/samtools
RUN cd samtools && make
# MOVE Directories INTO $HOME/tool
RUN mkdir /root/tools
RUN mv samtools /root/tools
RUN mv htslib /root/tools
RUN mv bcftools /root/tools
# Augustus
RUN wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.3.1.tar.gz
RUN tar -xzf augustus-3.3.1.tar.gz
RUN echo 'COMGENEPRED = true' >> augustus-3.3.1/common.mk
RUN cd augustus-3.3.1 && make
#### Unsure pathing ####
# RUN PATH=$PATH:~/augustus/bin:~/augustus/scripts
# RUN export AUGUSTUS_CONFIG_PATH=/my_path_to_AUGUSTUS/augustus/config/
RUN cd augustus-3.3.1/src && make clean all

# HDF5
RUN wget http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz
RUN tar xzf hdf5-1.10.1.tar.gz
RUN cd hdf5-1.10.1 && ./configure --prefix=/usr
RUN cd hdf5-1.10.1 && make && make install
# sonLib
RUN git clone git://github.com/benedictpaten/sonLib.git
# HAL
RUN git clone git://github.com/glennhickey/hal.git
RUN /bin/bash -c "pushd sonLib && make && popd"
RUN cd hal && make
#### Unsure pathing ####
# RUN export PATH=hal/bin:${PATH} 

# LibBigWig
RUN git clone https://github.com/dpryan79/libBigWig.git
RUN cd libBigWig && make install
# GSL
RUN wget ftp://www.mirrorservice.org/sites/ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz 
RUN tar -xvzpf gsl-latest.tar.gz
RUN cd gsl* && ./configure
RUN cd gsl* && make && make install
# WiggleTools
RUN git clone https://github.com/Ensembl/WiggleTools.git
RUN cd WiggleTools && make

# sambamba
RUN wget https://github.com/biod/sambamba/releases/download/v0.6.7/sambamba_v0.6.7_linux.tar.bz2
RUN tar xvjf sambamba_v0.6.7_linux.tar.bz2
