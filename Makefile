SHELL = /bin/bash

COMPILERC   = icc
COMPILERF   = mpif90 
FLAG  		= -openmp
#BUG		= -check all -traceback

### TAU UTIL ### 
#  #export PATH=$PATH:/home/guest/HyeonTae/Traiing/tau-2.28.1/bin
#  export TAU_PROFILE=/home/guest/HyeonTae/Training/tau-2.28.1/x86_64/lib/Makefile.tau-icpc-mpi-pdt-openmp
#  #TAU_PROFILE ?=/home/guest/HyeonTae/Training/tau-2.28.1/include/Makefile
#  include $(TAU_PROFILE) 
#  COMPILERF = TAU_PROFILE=$(TAU_PROFILE) /home/guest/HyeonTae/Training/tau-2.28.1/x86_64/bin/tau_f90.sh
#  COMPILERC = TAU_PROFILE=$(TAU_PROFILE) /home/guest/HyeonTae/Training/tau-2.28.1/x86_64/bin/tau_cc.sh

CODE 		= McBOX_new.out 

include FILE.list

.SUFFIXES: .f90 .o .c
.c.o:
	$(COMPILERC) $(FLAG) $(BUG) $(PROFILE) -c  $*.c

.f90.o:
	$(COMPILERF) $(FLAG) $(BUG) $(PROFILE) -c  $*.f90

	
%.o: %.mod

#a.f90.mod:
#	$(COMPILERF) $(FLAG) -c $<

# Make CODE
# all : $(CODE)

$(CODE) : $(FOBJS) $(COBJS)
	$(COMPILERF) $(FLAG) $(PROFILE) -o  $(CODE) $(FOBJS) $(COBJS)
	
# Make CODE - openmp
#omp : $(CODE)
#
#$(CODE) : $(FOBJS) $(COBJS)
#	$(COMPILERF) $(FLAG)  -o $(CODE) $(FOBJS) $(COBJS)
#
#noomp : $(CODE) 
#	$(COMPILERF) -o $(CODE) $(FOBJS) $(COBJS)

# Make clean
clean :
	rm  $(FOBJS) $(COBJS) $(CODE) *.mod  -f
