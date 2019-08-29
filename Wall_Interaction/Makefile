all: prog

INCLUDE=${HOME}/classicalDFT/include   -I../
DEBUG= -O3

# Name of text file containing build number.
BUILD_NUMBER_FILE=build-number.txt

# Include build number rules.
include local.mak

###########################################################################
# gcc version egcs-2.95.2 19991024 (release) on a Linux system
# Intel Pentium PC, use -m486 for 486 CPU or -mpentium for Pentium Systems
###########################################################################

## What processor?

I686=i686

MACHINE := $(shell uname -m)	

ifeq ($(strip $(MACHINE)),$(I686))
	PROCESSOR=pentium
	COPTS= $(DEBUG)  -march=${PROCESSOR} 
	CCC= g++
	CCOPTS= $(DEBUG) -std=gnu++11 -march=${PROCESSOR} 
	LIBINTEL= iode_ia32
else
	PROCESSOR=opteron
	COPTS= $(DEBUG)  -march=${PROCESSOR} -fopenmp
	CCC= g++
	CCOPTS= $(DEBUG) -std=gnu++11 -march=${PROCESSOR}  -fopenmp  -D USE_OMP 
	LIBINTEL =iode_intel64
endif


###########################################################################
### general rules							  #
###########################################################################

.SUFFIXES: .cpp .o .d

.c.o: 
	${CC} -c ${COPTS}  -I ${INCLUDE} $<

.cpp.o: 
	@echo 'building' $<
#	sed -i -r version_script &<
	${CCC} -c ${CCOPTS}  -I ${INCLUDE} $<

.cpp.d:
	g++ -MMD -E -march=${PROCESSOR} -std=gnu++11 -I ${INCLUDE} $< > /dev/null

.c.d:
	gcc -MMD -E -march=${PROCESSOR} -std=gnu++11 -I ${INCLUDE} $< > /dev/null

sourcescpp := $(wildcard ./*.cpp)
sourcesc := $(wildcard ./*.c)
sources := ${sourcescpp} ${sourcesc}

depscpp = $(sourcescpp:.cpp=.d)
depsc= $(sourcesc:.c=.d)

-include $(depscpp) $(depsc)

COMMON = ../

#LIBS= ${HOME}/dev/physics.a ${HOME}/dev/applications/Interface/interface.a
LIBS= 

prog : $(OBS) $(LIBS) $(BUILD_NUMBER_FILE) common
	${CCC} $(DEBUG)  -march=${PROCESSOR} $(BUILD_NUMBER_LDFLAGS) -o $(PROG)   ../DFT.o ../DFT_Surfactant.o ../FMT.o ../Species.o ../Density.o  ../Minimizer.o  ../options.o ../VDW1.o  ../Grace.o  $(OBS) $(LIBS)  -lgomp -lgsl -lgslcblas  -lgrace_np -larmadillo -lmgl  -lfftw3_omp -lfftw3 -lnlopt -lpthread -lm


.PHONY : clean
clean :
	-rm -f *.o
	-rm -f *.d
	-rm -f *~

.PHONY: common $(COMMON)
common: $(COMMON)
$(COMMON):
	$(MAKE) -C $@


# Include build number rules.
include buildnumber.mak



