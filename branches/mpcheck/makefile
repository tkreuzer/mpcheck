ARCH=ia64
CC=gcc
GMP=/users/spaces/logiciels/gmp-snap/$(ARCH)
MPFR=/users/spaces/logiciels/mpfr-2.0.1/$(ARCH)

# the following is for Solaris
#FENV=/local/logiciels/WorkShop-6/SUNWspro/WS6/include/cc
#M9X=/local/logiciels/WorkShop-6/SUNWspro/lib

# flags for Alpha: -mfp-rounding-mode=d -mieee-with-inexact
# need -ffloat-store on ia64-hpux for gcc
CFLAGS= -O2 -g -ffloat-store

mpcheck: mpcheck.c makefile
	$(CC) $(CFLAGS) -I$(GMP)/include -L$(GMP)/lib -I$(MPFR) -L$(MPFR) mpcheck.c -o mpcheck -lmpfr -lgmp -lm


