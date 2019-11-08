PROG = main
SRCS =  main.f90 mesh.f90 p4est_binding.f90 p4est_wrapper.c constantes.f90 vtk.f90
OBJS =  main.o mesh.o p4est_binding.o p4est_wrapper.o constantes.o vtk.o

INC  = -I/net/m/tverdier001/3A_CHP/TER/p4est-2.2/local/include/
LIBP = -L/net/m/tverdier001/3A_CHP/TER/p4est-2.2/local/lib/
LIBS = -lp4est -lsc

CC = gcc

F90 = gfortran
F90FLAGS = -O0 -pedantic -Wall -ffpe-trap=invalid,zero,overflow,underflow -g -fcheck=all -fbacktrace -fno-underscoring
CCFLAGS = -std=c99 -O0 -Wall -D_DEBUG -g

include /net/m/tverdier001/3A_CHP/TER/p4est-2.2/Makefile.p4est.mk

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) -g -o $@ $(OBJS) $(INC) $(LIBP) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .c

.c.o:
	$(CC) $(CCFLAGS) -c $<  $(INC)

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

main.o: constantes.o mesh.o vtk.o

p4est_binding.o: p4est_wrapper.o

mesh.o: p4est_binding.o constantes.o
