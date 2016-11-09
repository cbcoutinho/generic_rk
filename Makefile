# Makefile for generic_rk project

current_dir = $(shell pwd)
SRC=$(current_dir)/src
OBJ=$(current_dir)/obj
BIN=$(current_dir)/bin
TEST=$(SRC)/test
FORTRANLIB_SRC=$(SRC)/fortranlib/src

# Compiler
FF = gfortran
FFLAGS = -Wall -std=f2008 -Wextra -fPIC -fmax-errors=1 -Wimplicit-interface
# Debug flags:
FFLAGS += -g -fcheck=all -fbacktrace #-ffpe-trap=zero,overflow,underflow
# Release flags:
# FFLAGS += -O3 -march=native -ffast-math -funroll-loops

FLIBS = -lblas -llapack

.DEFAULT_GOAL := $(OBJ)/runge_kutta.o

# Dependencies of main program
objects=$(OBJ)/runge_kutta.o \
	$(OBJ)/rk_constants.o \
	$(OBJ)/lib_array.o \
	$(OBJ)/lib_constants.o

# Fortran library
$(OBJ)/lib_array.o: $(FORTRANLIB_SRC)/lib_array.f90
	$(FF) -J$(OBJ) -c -o $@ $<
$(OBJ)/lib_constants.o: $(FORTRANLIB_SRC)/lib_constants.f90
	$(FF) -J$(OBJ) -c -o $@ $<

# Modules
$(OBJ)/runge_kutta.o: $(SRC)/runge_kutta.f90 $(OBJ)/rk_constants.o
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<
$(OBJ)/rk_constants.o: $(SRC)/rk_constants.f90
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<

# Module containing test functions
$(OBJ)/misc.o: $(TEST)/misc.f90 $(OBJ)/lib_constants.o
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<

# Main test program
$(OBJ)/main.o: $(TEST)/main.f90 $(OBJ)/misc.o $(objects)
	$(FF) $(FFLAGS) -I$(OBJ) -c -o $@ $<
$(TEST)/main: $(OBJ)/main.o $(OBJ)/misc.o $(objects)
	$(FF) $(FFLAGS) -o $@ $+ $(FLIBS)

clean:
	rm -f $(OBJ)/*.o $(OBJ)/*.mod
	rm -f $(TEST)/main $(TEST)/*.o $(TEST)/*.mod
	rm -f data.out raw.out

test: $(TEST)/main
	$(TEST)/main

plot: clean $(TEST)/main test
	python plotter.py

debug: clean $(BIN)/main
	/usr/bin/valgrind --track-origins=yes --leak-check=full $(BIN)/main
