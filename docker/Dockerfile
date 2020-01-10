FROM ubuntu:18.04
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ARG git_commit

# update and install dependencies
RUN apt-get update && \
    apt-get -y install time git make wget autoconf gcc g++ sudo && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# install cmake
WORKDIR /tmp
RUN mkdir /opt/cmake && \
    wget https://cmake.org/files/v3.11/cmake-3.11.4-Linux-x86_64.sh && \
    sh /tmp/cmake-3.11.4-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && \
    ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake

# get shasta code
WORKDIR /opt
RUN git clone https://github.com/chanzuckerberg/shasta.git && \
    cd /opt/shasta && \
    git fetch && \
    git checkout $git_commit && \
    git submodule update --init && \
    rm -rf /opt/shasta/.git

# install prereqs and shasta
WORKDIR /opt/shasta/build
RUN ../scripts/InstallPrerequisites-Ubuntu.sh
RUN cmake .. && make all && make install/strip && /opt/shasta/build/shasta-install/bin/shasta -h

# setup entrypoint
COPY wrapper.sh /opt/wrapper.sh
WORKDIR /data
ENTRYPOINT ["sh", "/opt/wrapper.sh"]