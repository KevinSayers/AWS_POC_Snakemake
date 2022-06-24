FROM nciccbr/ccbr_ubuntu_base_20.04:v1.4

# install multiqc
RUN pip3 install multiqc==1.12

# cleanup etc
COPY Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2
WORKDIR /data2
RUN apt-get clean
