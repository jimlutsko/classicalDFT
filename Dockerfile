# Base: Get the base Ubuntu image from Docker Hub
FROM debian:stretch-slim

# To avoid interactive
ENV DEBIAN_FRONTEND=noninteractive

ENV APP_DIR=/usr/app
ENV LIB_DIR=dft_lib
ENV BUILD_SH=buildlib.sh

# Set up:
RUN apt-get -qq -y update \
    && apt-get -qq -y install curl \
    && apt-get install -y wget \
    && apt-get -qq install -y build-essential \
    && wget https://github.com/microsoft/CMake/releases/download/v3.17.3587832/cmake-3.17.3587832-MSVC_2-Linux-x64.sh \
    && chmod a+x cmake-3.17.3587832-MSVC_2-Linux-x64.sh \
    && ./cmake-3.17.3587832-MSVC_2-Linux-x64.sh --skip-license --prefix=/usr

# Dependencies:
RUN apt-get -qq -y install libgsl-dev \
    && apt-get -qq -y install libboost-all-dev \
    && apt-get -qq -y install libfftw3-dev libfftw3-doc \
    && apt-get -qq -y install grace \
    && apt-get -qq -y install gnuplot \
    && apt-get -qq -y install libarmadillo-dev \
    && apt-get -qq -y install libgtest-dev

# Config gtest:
RUN cd /usr/src/googletest/googletest \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && cp libgtest* /usr/lib/ \
    && cd .. \
    && rm -rf build

RUN mkdir /usr/local/lib/googletest \
    && ln -s /usr/lib/libgtest.a /usr/local/lib/googletest/libgtest.a \
    && ln -s /usr/lib/libgtest_main.a /usr/local/lib/googletest/libgtest_main.a

WORKDIR ${APP_DIR}
COPY ${LIB_DIR}/cmake ./${LIB_DIR}/cmake
COPY ${LIB_DIR}/core ./${LIB_DIR}/core
COPY ${LIB_DIR}/tests ./${LIB_DIR}/tests
COPY ${LIB_DIR}/examples ./${LIB_DIR}/examples
COPY ${LIB_DIR}/CMakeLists.txt ./${LIB_DIR}/CMakeLists.txt

COPY ${BUILD_SH} .

# Testing:
WORKDIR ${APP_DIR}
ENTRYPOINT ["./buildlib.sh", "test"]