#FC = gfortran
FC = ifort

FFLAGS =  -O3 -heap-arrays -fp-model precise -xSSE4.2
#FFLAGS = -O2 -fbounds-check

# The location of the libraries such as LAPACK and optimized BLAS
#
# 
MKLROOT =/opt/intel/composer_xe_2013.0.079/mkl
IFORTROOT= /opt/intel/composer_xe_2013.0.079/mkl/lib/intel64

#LAPACK = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm 

LAPACK=-Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm
LIn= -openmp -I$(MKLROOT)/include -L$(IFORTROOT)

#LAPACK=-L/home/lsg/lapack-3.5.0 -llapack -lrefblas

RANLIB = /usr/bin/ranlib
