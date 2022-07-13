FROM nciccbr/ccbr_ubuntu_base_20.04:v1.4

# install star
WORKDIR /opt2
RUN wget https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz
RUN tar -xzf 2.7.10a.tar.gz
WORKDIR /opt2/STAR-2.7.10a/source
RUN make STAR
ENV PATH="/opt2/STAR-2.7.10a/source":$PATH

# cleanup etc
COPY Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2
WORKDIR /data2
RUN apt-get clean
