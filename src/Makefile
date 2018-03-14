# Makefile created by mkmf $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 



include template


.DEFAULT:
	-echo $@ does not exist.
all: ../bin/imtest.app
FV_Solve.o: ./FV_Solve.f90 decimal.o FV_types.o scheme_utils.o linearsystems.o morse.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./FV_Solve.f90
FV_types.o: ./FV_types.f90 decimal.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./FV_types.f90
IMEX_scheme.o: ./IMEX_scheme.f90 decimal.o FV_types.o scheme_utils.o eno_weno.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./IMEX_scheme.f90
decimal.o: ./decimal.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./decimal.f90
eno_weno.o: ./eno_weno.f90 decimal.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./eno_weno.f90
linearsystems.o: ./linearsystems.f90 decimal.o morse.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./linearsystems.f90
main.o: ./main.f90 example1.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./main.f90
morse.o: ./morse.f90 decimal.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./morse.f90
rwdata.o: ./rwdata.f90 decimal.o util.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./rwdata.f90
scheme_utils.o: ./scheme_utils.f90 decimal.o FV_types.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./scheme_utils.f90
util.o: ./util.f90 decimal.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./util.f90
example1.o: ./examples/example1.f90 decimal.o rwdata.o FV_types.o IMEX_scheme.o FV_Solve.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./examples/example1.f90
SRC = ./rwdata.f90 ./scheme_utils.f90 ./decimal.f90 ./FV_types.f90 ./FV_Solve.f90 ./morse.f90 ./main.f90 ./linearsystems.f90 ./IMEX_scheme.f90 ./eno_weno.f90 ./util.f90 ./examples/example1.f90
OBJ = rwdata.o scheme_utils.o decimal.o FV_types.o FV_Solve.o morse.o main.o linearsystems.o IMEX_scheme.o eno_weno.o util.o example1.o
clean: neat
	-rm -f .../bin/imtest.app.cppdefs $(OBJ) ../bin/imtest.app *.mod
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
../bin/imtest.app: $(OBJ) 
	$(LD) $(OBJ) -o ../bin/imtest.app  $(LDFLAGS)
