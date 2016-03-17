#!/bin/sh

##############
# Shell script for installing nexus class library (NCL)
##############

# Downloading and unpacking ncl-2.1.18.tar.gz from http://sourceforge.net/projects/ncl/
cd $HOME
curl -L http://sourceforge.net/projects/ncl/files/latest/download > ncl-2.1.18.tar.gz
tar zxvf ncl-2.1.18.tar.gz

# Configuring and installing
cd $HOME/ncl-2.1.18
./configure
sudo make install

# Non-sudo installation
#cd $HOME/ncl-2.1.18
#./configure --prefix=$HOME/nclib
#make install

# Cleaning
cd ..
rm -R ncl-2.1.18/
rm ncl-2.1.18.tar.gz
