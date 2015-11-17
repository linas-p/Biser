[![Build Status](https://travis-ci.org/linas-p/Biser.svg?branch=master)](https://travis-ci.org/linas-p/Biser)

Requirements
====

 * c++11 i.e g++4.8/clang3.3 or better
 * R package

Install c++4.9(optional)
```sh
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test # gcc-4.9
sudo apt-get update
sudo apt-get install g++-4.9

```
To test cpp sample:

```sh
$ mkdir build
$ cd build
$ cmake ..
$ make
$ ./BiserLikeModel

```

To change params change file
```sh
$ config/params.json
```

R compatibility
====

```sh
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo add-apt-repository ppa:marutter/rdev -y # R 3.2
sudo apt-get update

apt-get install r-base
```

Install R packages(first time only)
```sh
R
install.packages("Rcpp", dependencies=TRUE)
install.packages("shiny", dependencies=TRUE)
quit()
```

To compile type
```sh
make
```

Usage 
```sh
R
shiny::runApp()
```



