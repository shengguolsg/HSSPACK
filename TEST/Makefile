
include ../make.inc

HSSPACK = ../LIB/libhsspack.a

# The test routines for different problems
RRLU_TEST = testCauchyHss_serial.o
BIDC_TEST = testvaldqds.o 
EIG_TEST = Testedc.o 
EIG_SB   = test_sbed.o
CHOL_TEST = testhsschol.o
BDC_TEST = test_bdc1.o
MMM_TEST = testfmml.o

all: rrlu bidc eig chol bdc mmm

rrlu: $(RRLU_TEST) 
	$(FC) $(FFLAGS) -o a.out $(RRLU_TEST) $(HSSPACK) $(LAPACK) $(LIn)
	./a.out

bidc: $(BIDC_TEST) 
	$(FC) $(FFLAGS) -o a.out $(BIDC_TEST) $(HSSPACK) $(LAPACK) $(LIn)
	./a.out

eig: $(EIG_TEST)
	$(FC) $(FFLAGS) -o a.out $(EIG_TEST) $(HSSPACK) $(LAPACK) $(LIn)
	./a.out

chol: $(CHOL_TEST)
	$(FC) $(FFLAGS) -o a.out $(CHOL_TEST) $(HSSPACK) $(LAPACK) $(LIn)
	./a.out

bdc: $(BDC_TEST)
	$(FC) $(FFLAGS) -o a.out $(BDC_TEST) $(HSSPACK) $(LAPACK) $(LIn)
	./a.out

sbdc: $(EIG_SB)
	$(FC) $(FFLAGS) -o a.out $(EIG_SB) $(HSSPACK) $(LAPACK) $(LIn)
	./a.out

mmm: $(MMM_TEST)
	$(FC) $(FFLAGS) $(LIn) -o a.out $(MMM_TEST) $(HSSPACK) $(LAPACK) $(LIn)
	./a.out


%.mod: %.f90
	$(FC) -c $(FFLAGS) $(LIn) *.f90

%.o: %.f90
	$(FC) -c $(FFLAGS) $(LIn) $<

%.o: %.f
	$(FC) -c $(FFLAGS) $(LIn) $<


#Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD a.out
veryclean:
	rm *~ 
