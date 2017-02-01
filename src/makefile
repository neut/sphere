#Makefile for parallel compiling
#=====================================================
# Note: use -C=all to catch array out of bounds errors
#============================================================================ 
# Compiler options
#============================================================================
# IFORT
#FC = ifort -r4 -O3 -ip -ipo -fpp -xAVX -axAVX

#FC += -openmp
#FC += -debug all -warn all -check all -g -traceback
#FC += -fpe0 -ftrapuv

# GNU
FC = gfortran -Ofast -cpp 
#-fdefault-real-8
FC += -fopenmp

# MINGW
#FC = /usr/local/i386-mingw32-4.3.0/bin/i386-mingw32-gfortran

#FC += -g -fbacktrace -fbounds-check -Wall
#FC += -DDEBUG

LINKS =
#============================================================================ 
# make PROGRAM   
#============================================================================

PROGRAM = main 

OBJECTS = find_rotation.o param.o rotated_pattern.o pattern.o \
          compare_frames.o projection.o load_img.o drawpattern.o nelmin.o \
          solver.o 

MODULES = param.o

#============================================================================ 
# Linking    
#============================================================================

$(PROGRAM) : $(MODULES) $(OBJECTS)
	$(FC) $(OBJECTS) $(LINKS) -o $@

#============================================================================
#  Dependencies
#============================================================================

param.o: param.F
	$(FC) -c  $(OP_COMP) param.F

%.o:   %.F $(MODULES)
	$(FC) -c $(OP_COMP) $< 
        

#============================================================================
#  Clean up
#============================================================================

clean:

	rm *.o *.mod



