#
#
# Top Level Makefile for HSSPACK
# Version 0.0.1
# March 2014
#

include make.inc

all: hsspack_lib 
# hsspack_testing

clean: cleanlib cleantesting 

hsspack_lib: 
	(cd SRC; $(MAKE)  libhsspack.a )

hsspack_testing:
	(cd TEST; $(MAKE) )

#Utility targets
.PHONY: clean veryclean

cleanlib:
	( cd SRC; $(MAKE) clean ) 
	( cd LIB; rm libhsspack.a )

cleantesting: 
	( cd TEST; $(MAKE) clean ) 

veryclean:
	rm *~
	( cd SRC;  $(MAKE) veryclean )
	( cd TEST; $(MAKE) veryclean )
