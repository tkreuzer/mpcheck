CC=gcc

# GMP is the directory where GMP is installed
# gmp.h should be in $(GMP)/include and libgmp.a in $(GMP)/lib
GMP=/usr

# MPFR is the directory where MPFR is installed
# mpfr.h should be in $(MPFR)/include and libmpfr.a in $(MPFR)/lib
MPFR=/usr

# the following is for Solaris
#FENV=/local/logiciels/WorkShop-6/SUNWspro/WS6/include/cc
#M9X=/local/logiciels/WorkShop-6/SUNWspro/lib

# flags for Alpha: -mfp-rounding-mode=d -mieee-with-inexact
# need -ffloat-store on ia64-hpux for gcc
CFLAGS= -O2 -g -ffloat-store

mpcheck: mpcheck.c makefile
	$(CC) $(CFLAGS) -I$(GMP)/include -L$(GMP)/lib -I$(MPFR)/include -L$(MPFR)/lib mpcheck.c -o mpcheck -lmpfr -lgmp -lm

clean:
	rm mpcheck


