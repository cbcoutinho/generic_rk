# Makefile for generic_rk project

current_dir = $(shell pwd)
SRC		= $(current_dir)/src
OBJ		= $(current_dir)/obj
BIN		= $(current_dir)/bin
TEST	= $(SRC)/test

# Compiler
FF = gfortran
FFLAGS = -Wall -std=f2008 -Wextra -fPIC -fmax-errors=1 -Wimplicit-interface
# Debug flags:
FFLAGS += -g -fcheck=all -fbacktrace #-ffpe-trap=zero,overflow,underflow
# Release flags:
# FFLAGS += -O3 -march=native -ffast-math -funroll-loops

# FLIBS = -L/home/chris/Software/intel/mkl/lib/intel64/
FLIBS += -lblas -llapack

.DEFAULT_GOAL := $(OBJ)/runge_kutta.o

# Dependencies of main program
objects = $(OBJ)/runge_kutta.o \
	$(OBJ)/rk_constants.o

# Modules
$(OBJ)/runge_kutta.o: $(SRC)/runge_kutta.f90 $(OBJ)/rk_constants.o
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<
$(OBJ)/rk_constants.o: $(SRC)/rk_constants.f90
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<

clean:
	rm -f $(OBJ)/*.o $(OBJ)/*.mod
	$(MAKE) -C $(TEST) clean

test: $(OBJ)/runge_kutta.o
	$(MAKE) -C $(TEST)

plot: $(OBJ)/runge_kutta.o
	$(MAKE) -C $(TEST) plot
