#!/bin/tcsh -f

gcc -I/usr/local/gsl/include -c $1.c
gcc -L/usr/local/gsl/lib -o $2 $1.o -lgsl -lgslcblas -lm
