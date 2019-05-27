


SHELL = /bin/bash


COMPILERC   = icc
#COMPILERF   = ifort 
COMPILERF   = mpif90 
FLAG  		= -openmp
BUG		= -check all -traceback


CODE 		= McCrit_beta.out 

include FILE.list

.SUFFIXES: .f90 .o .c
.c.o:
	$(COMPILERC) $(FLAG) -c  $*.c

.f90.o:
	$(COMPILERF) $(FLAG) $(BUG) -c  $*.f90

	
%.o: %.mod

#a.f90.mod:
#	$(COMPILER) $(FLAG) -c $<

# Make CODE
$(CODE) : $(FOBJS) $(COBJS)
	$(COMPILERF) $(FLAG) -o  $(CODE) $(FOBJS) $(COBJS)
	
## Make CODE - openmp
#omp : $(CODE)
#
#$(CODE) : $(FOBJS) $(COBJS)
#	$(COMPILERF) $(FLAG) -o  $(CODE) $(FOBJS) $(COBJS)

# Make clean
clean :
	rm  $(FOBJS) $(COBJS)  *.mod  -f
