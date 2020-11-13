#Build with:
#sudo docker build -t davidebolo1993/tricolor .

FROM ubuntu:20.04

# File author/maintainer info
MAINTAINER Davide Bolognini <davidebolognini7@gmail.com>

# Install dependencies
ENV DEBIAN_FRONTEND=noninteractive 
RUN apt-get update && apt-get install -y nano curl git build-essential g++ cmake libz-dev libcurl4-openssl-dev libssl-dev libbz2-dev  liblzma-dev libncurses5-dev && apt-get clean
RUN curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda
RUN conda create -y -n tricolorenv python=3.8
RUN echo "source activate tricolorenv" > ~/.bashrc
ENV PATH /miniconda/envs/tricolorenv/bin:$PATH

#get htslib
RUN curl -LO https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 && tar -vxjf htslib-1.11.tar.bz2 && cd htslib-1.11 && make && make install
#get samtools
RUN curl -LO https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && tar -vxjf samtools-1.11.tar.bz2 && cd samtools-1.11 && make && make install
#get bcftools
RUN curl -LO https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2 && tar -vxjf bcftools-1.11.tar.bz2 && cd bcftools-1.11 && make && make install

#install bedtools through miniconda
RUN conda install -y -n tricolorenv -c bioconda bedtools

#get TRiCoLOR
RUN git clone --recursive https://github.com/davidebolo1993/TRiCoLOR.git && cd TRiCoLOR && ./configure && python setup.py install

#Pull with:
#sudo docker pull davidebolo1993/tricolor

#Then run:
#sudo docker run davidebolo1993/tricolor TRiCoLOR --help

#Or load the environment
#sudo docker run -ti davidebolo1993/tricolor
#$(tricolorenv) TRiCoLOR --help
