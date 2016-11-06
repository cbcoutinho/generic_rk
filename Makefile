# Makefile for generic_rk project

current_dir = $(shell pwd)
SRC=$(current_dir)/src
OBJ=$(current_dir)/obj
BIN=$(current_dir)/bin

# Compiler
FF = gfortran
FFLAGS = -Wall -std=f2008 -Wextra -fPIC -fmax-errors=1 -Wimplicit-interface
# Debug flags:
FFLAGS += -O0 -g -fcheck=all -fbacktrace #-ffpe-trap=zero,overflow,underflow
# Release flags:
# FFLAGS += -O3 -march=native -ffast-math -funroll-loops

FLIBS = -lblas -llapack

.DEFAULT_GOAL := $(BIN)/main

# Dependencies of main program
objects=$(OBJ)/misc.o \
	$(OBJ)/runge_kutta.o

# Modules
$(OBJ)/misc.o: $(SRC)/misc.f90
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<
$(OBJ)/runge_kutta.o: $(SRC)/runge_kutta.f90
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
