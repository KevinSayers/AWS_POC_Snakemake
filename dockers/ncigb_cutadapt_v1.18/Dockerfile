FROM nciccbr/ccbr_ubuntu_base_20.04:v1.4

# install cutadapt
RUN pip3 install cutadapt==1.18
COPY TruSeq_and_nextera_adapters.consolidated.fa /opt2

# cleanup etc
COPY Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2
WORKDIR /data2
RUN apt-get clean
