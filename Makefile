
F90 = gfortran

FCFLAGS = -O3 -cpp

# required libraries
LIB =-L/mnt/data/home/xieyc/lib/lapack-3.12.0/build/lib -llapack -lblas

# required includes
INCLUDE=

#									     #
#----------------------------------------------------------------------------#
#									     #

# object to build
OBJ = Ellison.o Guillou.o Lawrence.o Liu.o Missner.o Ocean_Permittivity.o\
	foam_emis_model.o main.o


# RULES
#
default: foam_test


foam_test: $(OBJ)
		$(F90) -o $@ $(FCFLAGS) $(OBJ) $(LIB) $(LFLAGS)


clean:
	rm -f foam_test *.o *.mod *.MOD

#									     #
#----------------------------------------------------------------------------#
# GENERIC RULES                                                              #


%.o : %.f90
	$(F90) $(FCFLAGS) $(INCLUDE) -c $< -o $@

%.o : %.F90
	$(F90) $(FCFLAGS) $(INCLUDE) -c $< -o $@

%.o : %.F
	$(F90) $(FCFLAGS) $(INCLUDE) -c $< -o $@

%.o : %.f
	$(F90) $(FCFLAGS) $(INCLUDE) -c $< -o $@

