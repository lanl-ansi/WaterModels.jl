#!/bin/bash

# Obtain and compile Bonmin.
wget https://www.coin-or.org/Tarballs/Bonmin/Bonmin-1.8.6.tgz
tar -xzf Bonmin-1.8.6.tgz && rm Bonmin-1.8.6.tgz && cd Bonmin-1.8.6/ThirdParty
cd Blas && ./get.Blas && cd ..
cd Lapack && ./get.Lapack && cd ..
cd Metis && ./get.Metis && cd ..
cd Mumps && ./get.Mumps && cd ..
cd .. && mkdir build && cd build && ../configure
make && make install && export PATH="${PATH}:${PWD}/bin"
