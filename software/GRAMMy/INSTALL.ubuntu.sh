#!/bin/bash

#Compiler Dependencies:
sudo apt-get install build-essential libstdc++6 gfortran
#Python Dependencies:
sudo apt-get install python2.7 python-numpy python-scipy python-setuptools python-dev python-biopython
#Get Latest Release:
wget https://bitbucket.org/charade/grammy/get/release.tar.gz
#Install:
tar -zxvf charade-grammy-release.tar.gz
cd charade-grammy-release
python setup.py install

