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
# FFLAGS += -g -fcheck=all -fbacktrace #-ffpe-trap=zero,overflow,underflow
# Release flags:
FFLAGS += -O3 -march=native -ffast-math -funroll-loops

# export OMP_NUM_THREADS=4
FLIBS = -lblas -llapack -fopenmp

# Modules
runge_kutta.o: $(SRC)/runge_kutta.f90 rk_constants.o nonlinalg.o
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $(OBJ)/$@ $< $(FLIBS)
rk_constants.o: $(SRC)/rk_constants.f90
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $(OBJ)/$@ $<
linalg.o: $(SRC)/linalg.f90
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $(OBJ)/$@ $< $(FLIBS)
nonlinalg.o: $(SRC)/nonlinalg.f90 linalg.o
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $(OBJ)/$@ $< $(FLIBS)

clean:
	rm -f $(OBJ)/*.o $(OBJ)/*.mod
	$(MAKE) -C $(TEST) clean

test: clean runge_kutta.o nonlinalg.o linalg.o
	$(MAKE) -C $(TEST)

debug: clean runge_kutta.o nonlinalg.o linalg.o
	$(MAKE) -C $(TEST) debug

plot: clean runge_kutta.o
	$(MAKE) -C $(TEST) plot
