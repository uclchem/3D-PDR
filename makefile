#Makefile for 3DPDR
#Written by T. G. Bisbas
#University College London

#compiler options
F90               = gfortran
CPPFLAGS          = -cpp
INCLUDES          = -I/home/tb/Codes/sundials/include
LIBRARIES	  = -L/home/tb/Codes/sundials/lib 
LIBS	          = -lsundials_cvode -lsundials_nvecserial -lm
OPTIMISE          = 0
OPENMP            = 1
#main code
DIMENSIONS        = 1
NETWORK		  = REDUCED
DUST              = 2
GUESS_TEMP        = 1
THERMALBALANCE    = 1
TEMP_FIX          = 1
CO_FIX            = 1
H2FORM		  = 0

#Compilers options and flags

#F95 compiler
ifeq ($(F90),f95) 
  ifeq ($(OPENMP),1)
    OPT += -fopenmp -DOPENMP
  endif
#ifort compiler
else ifeq ($(F90),ifort)
  ifeq ($(OPENMP),1) 
    CFLAGS += -DOPENMP
    OPT += -openmp
  endif
#gfortran compiler
else ifeq ($(F90),gfortran)
  ifeq ($(OPENMP),1)
    CFLAGS += -DOPENMP
    OPT += -fopenmp
  endif
#gfortran44 compiler
else ifeq ($(F90),gfortran44)
  ifeq ($(OPENMP),1)
    CFLAGS += -DOPENMP
    OPT += -fopenmp
  endif
endif

#Other flags
ifeq ($(THERMALBALANCE),1)
  CFLAGS += -DTHERMALBALANCE
endif
ifeq ($(GUESS_TEMP),1)
  CFLAGS += -DGUESS_TEMP
endif
ifeq ($(DIMENSIONS),1)
  CFLAGS += -DPSEUDO_1D
endif
ifeq ($(DIMENSIONS),2)
  CFLAGS += -DPSEUDO_2D
endif
ifeq ($(OPTIMISE),0)
  OPT += -O0
endif
ifeq ($(OPTIMISE),1)
  OPT += -O1
endif
ifeq ($(OPTIMISE),2)
  OPT += -O2
endif
ifeq ($(OPTIMISE),3)
  OPT += -O3
endif
ifeq ($(OPTIMISE),4)
  OPT += -fast
endif

ifeq ($(DUST),1)
  CFLAGS += -DDUST
endif
ifeq ($(DUST),2)
  CFLAGS += -DDUST2
endif
ifeq ($(CO_FIX),1)
  CFLAGS += -DCO_FIX
endif
ifeq ($(H2FORM),1)
  CFLAGS += -DH2FORM
endif

ifeq ($(TEMP_FIX),1)
  CFLAGS += -DTEMP_FIX
endif

ifeq ($(NETWORK),REDUCED)
  CFLAGS += -DREDUCED
endif
ifeq ($(NETWORK),FULL)
  CFLAGS += -DFULL
endif
ifeq ($(NETWORK),MYNETWORK)
  CFLAGS += -DMYNETWORK
endif



MODULE_OBJ += definitions.o healpix_types.o modules.o 
CODE_OBJ += healpix.o input_parameters.o solvlevpop.o \
read_species.o read_rates.o heapsort.o calc_reac_rates.o shield.o \
spline.o escape_probability.o eval_points.o\
calculate_abundances.o sub_calculate_heating.o\
find_Ccoeff.o read_input.o 
ifeq ($(NETWORK),REDUCED)
CODE_OBJ += odes_reduced.o jacobian_reduced.o
endif
ifeq ($(NETWORK),FULL)
CODE_OBJ += odes_full.o jacobian_full.o
endif
ifeq ($(NETWORK),MYNETWORK)
CODE_OBJ += odes_mynetwork.o jacobian_mynetwork.o
endif
ifeq ($(DUST),2)
CODE_OBJ += dust_t.o
endif
ifeq ($(H2FORM),1)
CODE_OBJ += h2_form.o
endif

OBJ += $(MODULE_OBJ) $(CODE_OBJ)
#OPT += -check bounds

%.o: %.F90 definitions.o healpix_types.o modules.o
	$(F90) $(CPPFLAGS) $(OPT) $(CFLAGS) -c $<

%.o: %.F
	$(F90) $(OPT) $(CFLAGS) -c $<

ifeq ($(OPENMP),1) 
%.o: %.c
	gcc -fopenmp -O0 -Wall $(INCLUDES) $(CFLAGS) -c $<
else
%.o: %.c
	gcc -O0 $(INCLUDES) $(CFLAGS) -c $<
endif

3DPDR :: $(OBJ) 3DPDR.o
	$(F90) $(OPT) $(CFLAGS) $(LIBRARIES) -o 3DPDR $(OBJ) 3DPDR.o $(LIBS)
clean :: 
	rm *.mod *.o cvode.log main.log fort.*
compress ::
	tar cvzf 3DPDR.tgz 12c.dat 12c+.dat 16o.dat 3DPDR.F90 analyse_chem.F90 calc_reac_rates.F90 calculate_abundances.c 12co.dat definitions.F90 dust_t.F90 escape_probability.F90 eval_points.F90 find_Ccoeff.F90 h2_form.F90 healpix.F90 healpix_types.F90 heapsort.F90 input_parameters.F90 jacobian_full.c jacobian_mynetwork.c jacobian_reduced.c makefile modules.F90 odes_full.c odes_mynetwork.c odes_reduced.c params.dat rates_full.d rates_mynetwork.d rates_reduced.d read_input.F90 read_rates.F90 read_species.F90 shield.F90 solvlevpop.F90 species_full.d species_mynetwork.d species_reduced.d spline.F90 sub_calculate_heating.F90 sphere.dat 12c+_nometa.dat 1Dn*.dat

definitions.o : definitions.F90
	$(F90) $(CPPFLAGS) $(OPT) $(CFLAGS) -c $<

healpix_types.o : healpix_types.F90
	$(F90) $(CPPFLAGS) $(OPT) $(CFLAGS) -c $<

modules.o : modules.F90
	$(F90) $(CPPFLAGS) $(OPT) $(CFLAGS) -c $<

healpix_types.o : definitions.o
modules.o : definitions.o healpix_types.o
