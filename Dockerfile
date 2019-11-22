FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

MAINTAINER mm6@sanger.ac.uk

RUN apt-get -qq update
RUN apt-get install --no-install-recommends -y
RUN apt-get install --quiet --assume-yes autoconf automake make clang git curl zlib1g-dev liblzma-dev libbz2-dev pkg-config openssl libssl-dev
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
RUN . $HOME/.cargo/env
RUN export PATH="$HOME/.cargo/bin:$PATH"
RUN mkdir -p /opt
ENV  BUILD_DIR /opt/genedbot_rs
COPY . ${BUILD_DIR}
WORKDIR ${BUILD_DIR}
RUN ~/.cargo/bin/cargo build --release

ENV PATH=$PATH:${BUILD_DIR}/target/release
CMD genedbot
