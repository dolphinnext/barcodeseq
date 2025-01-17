# Don't upgrade nfcore/base, it creates "Kernel too old" error for singularity (because of the debian image)
FROM nfcore/base:1.7 

LABEL author="alper.kucukural@umassmed.edu" description="Docker image containing all requirements for the dolphinnext/barcodseq pipeline"

RUN apt-get update && apt-get install -y gcc libltdl7 libtbb-dev wget
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN mkdir -p /project /nl /mnt /share
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-x64-linux.tar.gz -O /tmp/blast.tar.gz && tar zxvf /tmp/blast.tar.gz -C /opt/ && rm /tmp/blast.tar.gz

ENV PATH /opt/conda/envs/dolphinnext-barcodeseq-1.0/bin:$PATH:/opt/ncbi-blast-2.11.0+/bin
COPY install_packages.R /
RUN Rscript /install_packages.R