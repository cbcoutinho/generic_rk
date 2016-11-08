# Makefile for generic_rk project

current_dir = $(shell pwd)
SRC=$(current_dir)/src
OBJ=$(current_dir)/obj
BIN=$(current_dir)/bin
FORTRANLIB_SRC=$(current_dir)/src/fortranlib/src

# Compiler
FF = gfortran
FFLAGS = -Wall -std=f2008 -Wextra -fPIC -fmax-errors=1 -Wimplicit-interface
# Debug flags:
FFLAGS += -g -fcheck=all -fbacktrace #-ffpe-trap=zero,overflow,underflow
# Release flags:
# FFLAGS += -O3 -march=native -ffast-math -funroll-loops

FLIBS = -lblas -llapack

.DEFAULT_GOAL := $(BIN)/main

# Dependencies of main program
objects=$(OBJ)/misc.o \
	$(OBJ)/runge_kutta.o \
	$(OBJ)/rk_constants.o \
	$(OBJ)/lib_array.o \
	$(OBJ)/lib_constants.o

# Fortran library
$(OBJ)/lib_array.o: $(FORTRANLIB_SRC)/lib_array.f90
	$(FF) -J$(OBJ) -c -o $@ $<
$(OBJ)/lib_constants.o: $(FORTRANLIB_SRC)/lib_constants.f90
	$(FF) -J$(OBJ) -c -o $@ $<

# Modules
$(OBJ)/misc.o: $(SRC)/misc.f90 $(OBJ)/lib_constants.o
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<
$(OBJ)/runge_kutta.o: $(SRC)/runge_kutta.f90 $(OBJ)/rk_constants.o
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<
$(OBJ)/rk_constants.o: $(SRC)/rk_constants.f90
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<

# Main program
$(OBJ)/main.o: $(SRC)/main.f90 $(objects)
	$(FF) $(FFLAGS) -I$(OBJ) -c -o $@ $<
$(BIN)/main: $(OBJ)/main.o $(objects)
	$(FF) $(FFLAGS) -o $@ $+ $(FLIBS)

clean:
	rm -f $(OBJ)/*.o $(OBJ)/*.mod $(BIN)/main

run: $(BIN)/main
	$(BIN)/main

plot: clean $(BIN)/main run
	python plotter.py

debug: clean $(BIN)/main
	/usr/bin/valgrind --track-origins=yes --leak-check=full $(BIN)/main
