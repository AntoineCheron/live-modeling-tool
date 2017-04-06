FROM python
LABEL name="Python OSUR"
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update
RUN apt-get install -y apt-utils
RUN pip3 install matplotlib numpy scipy cython
RUN pip3 install flopy
RUN apt-get install -y cython gfortran cmake unzip
RUN mkdir sundials
WORKDIR sundials
RUN mkdir srcdir builddir instdir
WORKDIR srcdir
RUN wget http://computation.llnl.gov/projects/sundials/download/sundials-2.7.0.tar.gz
RUN tar -xzvf sundials-2.7.0.tar.gz
RUN rm sundials-2.7.0.tar.gz
WORKDIR ../instdir
RUN cmake -DCMAKE_INSTALL_PREFIX=/sundials/instdir -DEXAMPLES_INSTALL_PATH=/sundials/instdir/examples ../srcdir/sundials-2.7.0
RUN make
RUN make install
WORKDIR ../../
RUN mkdir assimulo
WORKDIR assimulo
RUN wget https://pypi.python.org/packages/4c/c0/19a54949817204313efff9f83f1e4a247edebed0a1cc5a317a95d3f374ae/Assimulo-2.9.zip#md5=3f28fd98011d2ec7a01703a1ef1dff45
RUN unzip Assimulo-2.9.zip
RUN rm Assimulo-2.9.zip
WORKDIR Assimulo-2.9
RUN python setup.py install --sundials-home=/sundials/instdir
ENV LD_LIBRARY_PATH=/sundials/instdir/lib
COPY entrypoint.sh /
COPY hs1d/src_python/ /app/src_python/
COPY hs1d/test_case/ /app/test_case/
WORKDIR /
