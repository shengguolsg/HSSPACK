
include ../make.inc

MAKEFILE = Makefile

HSSMODULES = comp_Wtime.o aux_hss.o BasicMM.o ConstructHssd.o transp2d.o torthsvd.o \
             CauchyHssvd_VPS.o DHSSMML.o DHSSMM.o SvdHssUpdat.o testOrthAll.o

CHOLOBJS = mat2hsssvd.o HssCholSolver.o 

SVDOBJS  = CompSvdup.o mdlasd31.o mdlasd20.o

RRLUOBJS = Cauchylike.o CauchyHss_CP.o CauchyHss_VP2.o cauchyhssmm.o

BIDCOBJS = mdbdsdc.o mdlasd0.o mdlasd1.o mdlasd2.o mdlasd3.o mdlasd8.o 

EIGOBJS = CauchyHssEig.o mdlaed0.o mdlaed1.o mdlaed2.o mdlaed3.o mdlaed33.o \
          mdstedc.o

BDCOBJS = TestRoutine.o mdbdsvd.o mdbdsvd1.o mdbdsvd2.o mdbdsvd8.o mdbdsvd3.o \
          mdbdsvd22.o mdbdsvd4.o mdbdsvdt.o transp2di.o SFtransf.o

SBED = mdsbed1.o  mdsbed2.o  mdsbed22.o  mdsbed3.o  mdsbed4.o  mdsbed5.o  mdsbedc.o  mdsbevd.o


ALLROUTINES = $(HSSMODULES) $(SVDOBJS) $(RRLUOBJS) $(BIDCOBJS) $(EIGOBJS) $(CHOLOBJS) $(BDCOBJS) $(SBED)

libhsspack.a: $(ALLROUTINES)
	ar rc $@ $(ALLROUTINES)
	$(RANLIB) $@
	mv libhsspack.a ../LIB

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
