[![Build Status](https://travis-ci.org/linas-p/Biser.svg?branch=master)](https://travis-ci.org/linas-p/Biser)

Demo
===
Implementation can be access from <a href="https://elainas.shinyapps.io/biser">https://elainas.shinyapps.io/biser</a>.




Requirements
====


To test cpp sample:

```sh
$ mkdir build
$ cd build
$ cmake ..
$ make
$ ./BiserLikeModel

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


 * c++11 i.e g++4.9/clang3.3 or better(optional params from JSON)

```

To change params change file params.json
```sh
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test # gcc-4.9
sudo apt-get update
sudo apt-get install g++-4.9
$ config/params.json
```





