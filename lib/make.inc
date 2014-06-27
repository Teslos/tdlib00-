CC = gcc
CXX = g++
F77 = gfortran
ARCH = ar

CPPFLAGS = -O3 -s -I. -I../../include

F77FLAGS = -O3

ADDLIB = 

#change the path to LAPACK 
LAPACK = -llapack -lblas

G2C = -lf2c -lgfortran

