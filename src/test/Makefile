# Makefile for test of generic_rk project

current_dir = $(shell pwd)
SRC		= $(current_dir)/src
OBJ		= $(current_dir)/obj
BIN		= $(current_dir)/bin
RK_OBJ=$(current_dir)/../../obj
FORTRANLIB_SRC=$(current_dir)/../../src/fortranlib/src

FF = gfortran
FFLAGS = -Wall -std=f2008 -Wextra -fPIC -fmax-errors=1 -Wimplicit-interface
# Debug flags:
FFLAGS += -g -fcheck=all -fbacktrace #-ffpe-trap=zero,overflow,underflow
# Release flags:
# FFLAGS += -O3 -march=native -ffast-math -funroll-loops

FLIBS += -lblas -llapack

.DEFAULT_GOAL := test

rk_objects = $(RK_OBJ)/runge_kutta.o \
	$(RK_OBJ)/rk_constants.o
objects = $(OBJ)/lib_array.o \
	$(OBJ)/misc.o

$(OBJ)/lib_array.o: $(FORTRANLIB_SRC)/lib_array.f90
	$(FF) -J$(OBJ) -c -o $@ $<

$(OBJ)/misc.o: $(SRC)/misc.f90
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<

$(BIN)/main: $(SRC)/main.f90 $(objects) $(rk_objects)
	$(FF) $(FFLAGS) -I$(RK_OBJ) -I$(OBJ) -o $@ $+ $(FLIBS)

clean:
	rm -f $(OBJ)/*.o $(OBJ)/*.mod
	rm -f $(BIN)/main
	rm -f data.out raw.out plot.png

test: clean $(BIN)/main
	$(BIN)/main

plot: test
	python plotter.py
