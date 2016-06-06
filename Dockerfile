FROM feelpp/develop:latest
MAINTAINER Feel++ Support <support@feelpp.org>

USER feelpp
ENV HOME /home/feelpp

# VMTK
ENV VMTK_VERSION 1.3

WORKDIR /tmp
RUN git clone https://github.com/vmtk/vmtk.git

WORKDIR /tmp/vmtk
RUN git checkout refs/tags/v1.3

RUN mkdir -p /tmp/vmtk/build
WORKDIR /tmp/vmtk/build
RUN sudo apt-get update && sudo apt-get -y install libosmesa6-dev
RUN sudo cmake /tmp/vmtk -DSUPERBUILD_INSTALL_PREFIX=/usr/local
RUN sudo make -j 24

# Install AngioTK
COPY . /tmp/angiotk
RUN sudo chown -R feelpp /tmp/angiotk
RUN mkdir -p /tmp/angiotk/build

WORKDIR /tmp/angiotk/build
#RUN cmake /tmp/angiotk/

# COPY WELCOME $HOME/WELCOME
USER root
ENTRYPOINT ["/sbin/my_init","--quiet","--","sudo","-u","feelpp","/bin/sh","-c"]
CMD ["/bin/bash"]
