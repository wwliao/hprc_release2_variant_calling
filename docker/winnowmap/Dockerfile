FROM wwliao/hpp_base:latest
LABEL maintainer="wen-wei.liao@yale.edu"

USER root

# Winnowmap v2.03
WORKDIR /opt/winnowmap
RUN wget https://github.com/marbl/Winnowmap/archive/refs/tags/v2.03.tar.gz && \
    tar jxf v2.03.tar.gz && \
    rm v2.03.tar.gz && \
    cd Winnowmap-2.03 && \
    make -j8
ENV PATH="/opt/winnowmap/Winnowmap-2.03/bin:$PATH"

WORKDIR /data
