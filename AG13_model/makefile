SHELL = /bin/sh


FC = ifort
FL = ifort

FC_FLAG = -mcmodel=large -c -g 
FL_FLAG = -o




# include path
INCLUDE = -I /mnt/data/home/xieyc/bck/gsdart-mul-docker/app/netcdf-4/include

#-L/mnt/data/home/xieyc/lib/lapack-3.12.0/build/lib -llapack -lblas  --gfortan
# library path
LIBRARY = -L /mnt/data/home/xieyc/work/sciatran/lapack-3.12.0 -llapack -lblas \
		-L /mnt/data/home/xieyc/bck/gsdart-mul-docker/app/netcdf-4/lib -lnetcdff -lnetcdf


EXE = test
OBJ =Ellison.o Guillou.o Lawrence.o Liu.o Missner.o Ocean_Permittivity.o \
	mod_foam_emiss.o main.o


.SUFFIXES:
.SUFFIXES: .F90 .f90 .o .f


all: $(EXE)


$(EXE): $(OBJ)
	$(FL) $(FL_FLAG) $@ $(OBJ) ${LIBRARY} ${INCLUDE}

.f90.o:
	$(FC) ${LIBRARY} ${INCLUDE} $(FC_FLAG) $<

.PHONY:clean
clean:
	rm -f $(EXE) $(OBJ) *.mod *.o
