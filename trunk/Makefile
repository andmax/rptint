#
# Makefile for RPTINT
#

#SHELL=/bin/csh

CXX = g++
#CXX = icpc

LDIR = $(HOME)/lcgtk

KDIR = $(LDIR)/glslKernel

LIBDIR = -L$(LDIR)/lib

SRC =	src
OBJ =	obj

INCLUDES =	-I$(LDIR) -I$(KDIR) -Iinclude

OBJS =	$(OBJ)/ftrackball.o $(OBJ)/render_volume_gpu.o \
	$(OBJ)/transferFunction.o $(OBJ)/volume.o

SRCS =	$(OBJ)/ftrackball.cc $(OBJ)/render_volume_gpu.cc \
	$(OBJ)/transferFunction.cc $(OBJ)/volume.cc

APP = rptint

DEBUG_FLAGS = #-g
OPT_FLAGS = -O3 -ffast-math
ICPC_FLAGS = -D_GLIBCXX_GTHREAD_USE_WEAK=0 -pthread

CXX_FLAGS = -Wall -Wno-deprecated $(DEBUG_FLAGS) $(INCLUDES) $(OPT_FLAGS) $(ICPC_FLAGS)

LNK_FLAGS = -fPIC $(OPT_FLAGS) $(ICPC_FLAGS)

LIBS =	-lGLee -lglut -lGL -lglslKernel

#-----------------------------------------------------------------------------

$(APP): $(OBJS)
	@echo "Linking ..."
	$(CXX) $(LNK_FLAGS) -o $(APP) $(OBJS) $(LIBDIR) $(LIBS)

depend:
	rm -f .depend
	$(CXX) -M $(CXX_FLAGS) $(SRCS) > .depend

$(OBJ)/%.o: $(SRC)/%.cc
	@echo "Compiling ..."
	$(CXX) $(CXX_FLAGS) -c $< -o $@

clean:
	rm -f $(OBJ)/* $(SRC)/*~ $(APP) .depend

ifeq (.depend,$(wildcard .depend))
include .depend
endif
