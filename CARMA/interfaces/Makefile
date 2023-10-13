
INTERFACE_OBJ = carma_box.o

carma_box.o : carma_box.F90 carma_mod.mod carma_constants_mod.mod carma_precision_mod.mod
	$(FORTRAN) $(FFLAGS) -c $<
