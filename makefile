PROG = wells
# SPECIFY THE "C" and "Fortran" compilers
# Program needs GSLIB Library. Please specify the correct path where
# these libraries are installed.
ifeq ($(OSTYPE),linux)
	ifeq ($(LANL_ROOT),/usr/lanl) 
		DIRS = -I/home/monty/local/include/ -L/home/monty/local/lib
		GFORTRAN = /usr/lanl/bin/gfortran
		CC = /usr/lanl/bin/gcc
		LG = -lgfortran
	else
		DIRS = -I/opt/local/include/ -L/opt/local/lib
		CC = gcc
		LG = -lgfortran
	endif
else
        DIRS = -I/opt/local/include/ -L/opt/local/lib
	CC = gcc
	GFORTRAN = gfortran
endif
# Following flags are valid for gfortran 4.6 or higher
# CFLAGS = -O3 -msse2 -ftree-vectorize -ftree-vectorizer-verbose=4  -ffast-math -Wall $(DIRS) 
CFLAGS = -static -Wall $(DIRS) 
LDLIBS = -lgsl -lm -lgslcblas $(DIRS)

OBJ = wells.o invlap-lib.o complex_bessels.o
MOD = complex_bessel.mod
SOURCE = $(OBJ:%.o=%.c)
SOURCESTYLE = $(OBJ:%.o=%.c)

$(PROG): $(OBJ)
	$(GFORTRAN) $(OBJ) -o wells $(LDLIBS) 

wells.o: wells.c wells.h design.h
invlap-lib.o: invlap-lib.c wells.h design.h
complex_bessel.mod: complex_bessels.f90

complex_bessels.o: complex_bessels.f90 
	$(GFORTRAN) -c complex_bessels.f90 

clean:
	rm -f $(OBJ) $(MOD)

cleaner:
	rm -f $(PROG) $(OBJ) $(MOD)

astyle:
	astyle $(SOURCESTYLE)
	rm -f $(SOURCESTYLE:%c=%c.orig)

tar:
	tar -cvzf wells.tgz `hg st -c -m | awk '{print $$2}'` `find . \( -name "*.[ch]" -o -name "[Mm]akef*" -name "[Rr]eadme" \) -print | sed 's/\.\///'` .hg

tarf:
	tar -cvzf wells.tgz `hg st -c -m | awk '{print $$2}'` `find . \( -name "*.[ch]" -o -name "[Mm]akef*" -name "[Rr]eadme" \) -print | sed 's/\.\///'`
