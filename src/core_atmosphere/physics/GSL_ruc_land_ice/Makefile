.SUFFIXES: .F90 .o

all: dummy ruc_land_ice

dummy:
	echo "****** compiling RUC land and sea ice models ******"

OBJS = \
	module_ruc_land.o    \
	module_ruc_ice.o

land_ice: $(OBJS)
	ar -ru ./../libphys.a $(OBJS)

# DEPENDENCIES:
module_mp_thompson.o: \
	module_ruc_land.o \
	module_ruc_ice.o

clean:
	$(RM) *.f90 *.o *.mod
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F90.o:
ifeq "$(GEN_F90)" "true"
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I.. -I../physics_wrf -I../physics_mmm -I../../../framework -I../../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../physics_wrf -I../physics_mmm -I../../../framework -I../../../external/esmf_time_f90
endif
