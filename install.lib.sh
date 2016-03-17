#!/bin/sh

##############
# Installing morphy dependencies libraries
##############

#~~~~~~~~~
# NCL
#~~~~~~~~~

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

#~~~~~~~~~
# GSL
#~~~~~~~~~

# Downloading and unpacking
cd $HOME
curl -L http://ftp.heanet.ie/mirrors/gnu/gsl/gsl-2.1.tar.gz  > gsl-2.1.tar.gz
tar zxvf gsl-2.1.tar.gz

# Configuring and installing
cd $HOME/gsl-2.1
./configure
make
make check
sudo make install
make installcheck

# Cleaning
rm -R $HOME/ncl-2.1.18
rm $HOME/ncl-2.1.18.tar.gz
rm -R $HOME/gsl-2.1
rm $HOME/gsl-2.1.tar.gz
