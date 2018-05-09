#
# setup control
#

#CXX           = g++
#CXXFLAGS      = -O -Wall -fPIC
#LD            = g++
#LDFLAGS       = -O
#SOFLAGS       = -shared


DEBUG = TRUE

TOP     := $(shell pwd)/
OBJ     := $(TOP)obj/
BIN     := $(TOP)bin/
LIB     := $(TOP)lib/
SRC     := $(TOP)src/
MAIN    := $(TOP)main/
INC     := $(TOP)starincludes/
STAF    := $(TOP)staf/
STAFINC := $(STAF)inc/

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTINC       = $(shell root-config --incdir)

#
# set up dependencies
#
vpath %.h   $(ROOTINC):$(INC)
vpath %.inc $(INC):$(STAFINC)
vpath %.cc  $(SRC):$(MAIN)
vpath %.cxx  $(SRC):$(MAIN)
vpath %.c   $(SRC):$(MAIN)
vpath %.C   $(SRC):$(MAIN)
vpath %.F   $(SRC)
vpath %.o   $(BIN)
#	
# set up compilers
#	

CPP     = g++
CPPFLAGS = -O -Wall -fPIC -I $(ROOTINC) -I$(INC)
#CPPFLAGS = -m32 -O -Wall -fPIC -I $(ROOTINC) -I$(INC)
CC      = gcc
ifeq ($(OPT),TRUE)	
#CPPFLAGS += -g 
endif	
ifeq ($(DEBUG), TRUE)
#CPPFLAGS += -g
endif
#
CFLAGS   = -w -I$(INC) -I$(STAFINC) 
ifeq ($(OPT), TRUE)
#CFLAGS += -fast
endif
ifeq ($(DEBUG), TRUE)
#CFLAGS += -g	
endif	
#
FF      = f77 
FFLAGS   = -c -e -xl -PIC -DSunOS -Dsun4os5 -I$(INC) -I$(STAFINC) 
ifeq ($(OPT), TRUE)
FFLAGS += -fast
endif
ifeq ($(DEBUG), TRUE)
FFLAGS += -g
endif
#
#CERNLIBS := $(shell $(CERN_ROOT)/bin/cernlib  pawlib graflib packlib mathlib kernlib)
#

ifeq ($(PURIFY), TRUE)
LD = purify $(CPP)
else
LD     = $(CPP)
endif
#LDFLAGS = -t -z muldefs -DSUN -DSOLARIS -Dsun 
#LDFLAGS = -m32
ifeq ($(OPT), TRUE)
#LDFLAGS += -fast
endif
#
CSRCS := $(notdir $(wildcard $(SRC)*.c)) 
FSRCS := $(notdir $(wildcard $(SRC)*.F))
INCS  := $(notdir $(wildcard $(INC)*.h)) $(notdir $(wildcard $(INC)*.inc))
MAINS := $(notdir $(wildcard $(MAIN)*.c))
OBJS  := $(FSRCS:.F=.o) $(CSRCS:.c=.o)
#
NOOPTOBJS :=  $(BIN)Raw.o $(BIN)TrkSort.o
#	
SOFLAGS = -G

############## Make Executables ####################################
all: genMassHypo
#all: ds kaka pika pipi lambdac charged d0rho thetapp l1520 thetapp2 d0


genMassHypo : $(OBJ)genMassHypo.o $(OBJ)event.o $(OBJ)hit.o $(OBJ)material.o
	$(LD) $(LDFLAGS) $(LDINC) $(XM_LIBS) $^ $(LDLIB) $(ROOTLIBS) $(ROOTGLIBS) -o $(BIN)$(notdir $@)
	@echo


######################################################
$(OBJ)%.o : 	$(MAIN)%.c $(INCS)
	$(CC)     $(CFLAGS) -c $(MAIN)$(notdir $<) -o $(OBJ)$(notdir $@)
	@echo	
$(OBJ)%.o : 	$(MAIN)%.cc $(INCS)
	$(CPP)     $(CPPFLAGS) -c $(MAIN)$(notdir $<) -o $(OBJ)$(notdir $@)
	@echo	
$(OBJ)%.o : 	$(SRC)%.cxx $(INCS)
	$(CPP)     $(CPPFLAGS) -c $(SRC)$(notdir $<) -o $(OBJ)$(notdir $@)
	@echo	
$(OBJ)%.o : 	$(MAIN)%.C $(INCS)
	$(CPP)     $(CPPFLAGS) -c $(MAIN)$(notdir $<) -o $(OBJ)$(notdir $@)
	@echo	
$(OBJ)%.o : 	$(SRC)%.c $(INCS)
	$(CC)     $(CFLAGS) -c $(SRC)$(notdir $<) -o $(OBJ)$(notdir $@)
	@echo	
$(OBJ)%.o : 	$(SRC)%.C $(INCS)
	$(CPP)     $(CPPFLAGS) -c $(SRC)$(notdir $<) -o $(OBJ)$(notdir $@)
	@echo	
$(OBJ)%.o : 	$(SRC)%.F $(INCS)
	$(FF)     $(FFLAGS) -c $(SRC)$(notdir $<) -o $(OBJ)$(notdir $@)
	@echo	
foo:
	@echo $(CSRCS)
	@echo
	@echo $(FSRCS)
	@echo
	@echo $(MAINS)
	@echo
	@echo $(OBJS)
	@echo
	@echo $(INCS)
	@echo
clean:
	rm -f $(OBJ)*.o 
	rm -f $(BIN)*
