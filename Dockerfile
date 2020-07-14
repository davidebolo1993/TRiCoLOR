#Build with:
#sudo docker build -t davidebolo1993/tricolor .

FROM ubuntu:18.04

# File author/maintainer info
MAINTAINER Davide Bolognini <davidebolognini7@gmail.com>

# Install dependencies
RUN apt-get update && apt-get install -y nano curl git build-essential g++ cmake && apt-get clean
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
RUN bash Miniconda-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda
RUN conda create -y -n tricolorenv python=3.7
RUN echo "source activate tricolorenv" > ~/.bashrc
ENV PATH /miniconda/envs/tricolorenv/bin:$PATH
RUN conda install -y -n tricolorenv -c bioconda samtools bedtools bedops minimap2 bcftools pysam pyfaidx cyvcf2
RUN git clone --recursive https://github.com/davidebolo1993/TRiCoLOR.git && cd TRiCoLOR && ./configure && python setup.py install

#Pull with:
#sudo docker pull davidebolo1993/tricolor

#Then run:
#sudo docker run davidebolo1993/tricolor TRiCoLOR --help

#Or load the environment
#sudo docker run -ti davidebolo1993/tricolor
#$(tricolorenv) TRiCoLOR --help
