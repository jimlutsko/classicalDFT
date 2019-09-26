INCLUDE= ./include
SRC_DIR= ./src
OBJ_DIR= .
DEP_DIR= .

SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

all: $(OBJ_FILES)

DEBUG= -O3

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
	CCC= mpic++
	CCOPTS= $(DEBUG) -std=gnu++11 -march=${PROCESSOR} -MMD
#	LIBINTEL= iode_ia32
else
	PROCESSOR=opteron
	COPTS= $(DEBUG)  -march=${PROCESSOR} -fopenmp -MMD
	CCC= mpic++
	CCOPTS= $(DEBUG) -std=gnu++11 -march=${PROCESSOR}  -MMD -fopenmp  -D USE_OMP 
#	LIBINTEL =iode_intel64
endif

###########################################################################
### general rules							  #
###########################################################################

.SUFFIXES: .cpp .o .d

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	${CCC} $(CCOPTS) -I ${INCLUDE} -c -o $@ $<

$(DEP_DIR)/%.d: $(SRC_DIR)/%.cpp
	${CCC} $(CCOPTS) -I ${INCLUDE} -c -o $@ $<

.cpp.d:
	mpic++ -MMD -E -march=${PROCESSOR} -std=gnu++11 -I ${INCLUDE} $< > /dev/null

sourcescpp := $(wildcard ./src/*.cpp)
sourcesc := $(wildcard ./*.c)
sources := ${sourcescpp} ${sourcesc}

depscpp = $(sourcescpp:.cpp=.d)
depsc= $(sourcesc:.c=.d)

-include $(depscpp) $(depsc)

LIBS= 

.PHONY : clean
clean :
	-rm -f ${OBJ_DIR}/*.o
	-rm -f *.d
	-rm -f *~
