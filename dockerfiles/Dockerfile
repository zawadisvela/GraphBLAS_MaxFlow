#Build with
#docker build /path/to/this/folder -t hzsvela_groute
#Run with
#nvidia-docker run --rm -v /home/hzsvela/inputs/:/inputs/ -it hzsvela_groute /bin/bash
ARG DEBIAN_FRONTEND=noninteractive

FROM nvidia/cuda:10.2-devel

RUN apt-get update && apt-get install -y \
        wget \
        vim \
        build-essential \
        libboost-all-dev \
        python3 \
	valgrind \
        #devscripts \
        #debhelper \
        #fakeroot \
        #libssl-dev \
        git

RUN version=3.11 && \
        build=1 && \
        wget https://cmake.org/files/v$version/cmake-$version.$build.tar.gz && \
        tar -xzvf cmake-$version.$build.tar.gz && \
        cd cmake-$version.$build/ && \
        ./bootstrap && \
        make -j4 && \
        make install

RUN 	apt update && \
	apt install -y software-properties-common --fix-missing && \
	add-apt-repository -y ppa:ubuntu-toolchain-r/test && \
	apt update && \
	apt install -y gcc-10 && \
	apt install -y g++-10

RUN 	git clone --recursive https://github.com/DrTimothyAldenDavis/GraphBLAS.git && \
	cd GraphBLAS && \
    	make jobs=16 && \
	make install

COPY install_graphblas.sh /install_grapblas.sh

ENV CC=/usr/bin/gcc-10
ENV CXX=/usr/bin/g++-10
