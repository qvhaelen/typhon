#############################################################
###
###        PARAMETERS THAT HAVE TO BE DEFINED MANUALLY
###
#############################################################
default:
	@make typhon

##### PLEASE UNCOMMENT THE DESIRED PRECISION

#REALTYPE   = single
REALTYPE   = double

#############################################################
###
###        DIRECTORIES
###
#############################################################

objdir = ./OBJ/
srcdir = ./SRC/


#############################################################
###
###        PARAMETERS DEPENDING ON THE SIMPLE/DOUBLE
###
#############################################################

ifeq ($(REALTYPE),single) 
  MESSAGE = echo !!!!!!!!!!!!!!!!!!!!!!  CAUTION : This code is compiled in simple precision
endif
ifeq ($(REALTYPE),double)
  MESSAGE = echo !!!!!!!!!!!!!!!!!!!!!!  CAUTION : This code is compiled in double precision
endif

MESSAGE+=;echo !!!!!!!!!!!!!!!!!!!!!!  CAUTION : If you change the simple/double precision option, it is highly 
MESSAGE+=;echo !!!!!!!!!!!!!!!!!!!!!! '         ' recommended to recompile everything : use 'make clean_all'

###############################################
###  PARAMETERS FOR POWERBOOK-G4
###############################################
LIBDIR = $(HOME)/Unix/fftw/lib
CMPL   = gfortran  -O3 -ffree-line-length-none -fdefault-double-8 -fdefault-integer-8 -fdefault-real-8 -funroll-loops -g     
#-fbacktrace -Wall -v -da -Q 
#CMPL   = ifort -O1 -C -WB -traceback -g
DBG    = -w
 #ifeq ($(REALTYPE),double)
 # REALSIZE =-r8
 #  DEFINES= -DDOUBLE
 #endif
MESSAGE+=;echo !!!!!!!!!!!!!!!!!!!!!!  CAUTION : This code is for host = MacBook Pro
CMPLFLAGS = $(REALSIZE) $(DBG) -I$(objdir)
#CMPLFLAGS = $(REALSIZE) $(DBG) -I$(moddir)


#############################################################
###
###        COMPILE THE CODES 
###
#############################################################
	
vpath %.f90 ./SRC
vpath %.F90 ./SRC
vpath %.o ./OBJ
vpath %.mod ./OBJ

mod_src = common_mod.F90
mod_obj = $(patsubst %.F90,%.o,$(mod_src))

$(mod_obj): $(mod_src)
	$(CMPL) $(CMPLFLAGS) $(DEFINES) -c $^ -o $@
	@mv $@ $(objdir)$@
	@mv *.mod $(objdir)

typhon_sources = SSA_input_data.f90 SSA_scaling_set_up.f90 SSA_output_data.f90  SSA_central_core_partial_summation_DG.f90   initialisation_setup.f90 output_statistics.f90
typhon_objects = $(patsubst %.f90,%.o,$(typhon_sources))
$(typhon_objects) : %.o : %.f90 $(mod_obj) 
	$(CMPL) $(CMPLFLAGS) -c $(srcdir)$(patsubst %.o,%.f90,$@) -o $(objdir)$(notdir $@)

typhon	        :
	@$(MESSAGE)
	@make typhon
		
typhon		:	$(typhon_objects) typhon.f90
	$(CMPL) $(CMPLFLAGS) $(srcdir)typhon.f90 $(patsubst %.o,$(objdir)%.o,$(typhon_objects)) $(objdir)$(mod_obj) -o typhon

	@echo !!!!!!!!!!!!!!!!!!!!!!  SUCCESS.

clean		:
	@rm -f $(objdir)*.o *.o sph $(objdir)*.mod *.mod 
	@echo !!!!!!!!!!!!!!!!!!!!!!  CLEANING DONE.
	

