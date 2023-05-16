# MakeFile written by Luka Poniatowski 
# for line 3D based on Levin Heneckers original script

###########################HDF5 ENVIRONMENT VARIABLE#####################

ifndef LIB_HDF5
   $(info LIB_HDF5 Not found)
   $(error Undefined variable LIB_HDF5: export LIB_HDF5=<location of hdf5 lib>)
endif

#########################COMPILER ENVIRONMENT VARIABLE###################

ifndef COMPILER
   $(info COMPILER  Not defined)
   $(error Please define environment variable COMPILER: export COMPILER=<location of sutable compiler>)
endif

#########################################################################

#operating system
UNAME := $(shell uname)

.SUFFIXES: .f90

#####################GFORTRAN COMPILER OPTIONS ETC#######################

ifneq (,$(findstring gfortran,$(COMPILER)))
   F90 = COMPILER
   LD = COMPILER

   ifeq ($(UNAME), Linux)
      $(info OS detected: Linux)
      DIR_MOD_HDF5=$(LIB_HDF5)/include
      DIR_LIB_HDF5=$(LIB_HDF5)/lib                      \
                   $(LIB_HDF5)/lib/libhdf5_fortran.so   \
                   $(LIB_HDF5)/lib/libhdf5.so           \
                   $(LIB_HDF5)/lib/libhdf5_hl.so        \
                   $(LIB_HDF5)/lib/libhdf5hl_fortran.so
   else ifeq ($(UNAME), Darwin)
      $(info OS detected: MAC)
      DIR_MOD_HDF5=$(LIB_HDF5)/include
      DIR_LIB_HDF5=$(LIB_HDF5)/lib                         \
                   $(LIB_HDF5)/lib/libhdf5_fortran.dylib   \
                   $(LIB_HDF5)/lib/libhdf5.dylib           \
                   $(LIB_HDF5)/lib/libhdf5_hl.dylib        \
                   $(LIB_HDF5)/lib/libhdf5hl_fortran.dylib
   else
      $(error No OS detected, please have a look into the Makefile)
   endif

   OMP_FLAG=-fopenmp

   ifeq ($(DEBUG),true)
      DEBUG_FLAGS = -Og -Wall -Wextra -Warray-temporaries -Wconversion -pedantic-errors -fcheck=all -fbacktrace
   else
      DEBUG_FLAGS =
   endif

   CFLAGS = -ffree-line-length-none -c -O3 $(OMP_FLAG) $(DEBUG_FLAGS)
   MFLAGS = -J

   CDEFINED = true

endif

#####################IFORT COMPILER OPTIONS ETC#########################

ifneq (,$(findstring ifort,$(COMPILER)))

   F90 = COMPILER
   LD = COMPILER

   ifeq ($(UNAME), Linux)
      $(info OS detected: Linux)
      DIR_MOD_HDF5=$(LIB_HDF5)/include
      DIR_LIB_HDF5=$(LIB_HDF5)/lib                      \
                   $(LIB_HDF5)/lib/libhdf5_fortran.so   \
                   $(LIB_HDF5)/lib/libhdf5.so           \
                   $(LIB_HDF5)/lib/libhdf5_hl.so        \
                   $(LIB_HDF5)/lib/libhdf5hl_fortran.so
   else ifeq ($(UNAME), Darwin)
      $(info OS detected: MAC)
      DIR_MOD_HDF5=$(LIB_HDF5)/include
      DIR_LIB_HDF5=$(LIB_HDF5)/lib                         \
                   $(LIB_HDF5)/lib/libhdf5_fortran.dylib   \
                   $(LIB_HDF5)/lib/libhdf5.dylib           \
                   $(LIB_HDF5)/lib/libhdf5_hl.dylib        \
                   $(LIB_HDF5)/lib/libhdf5hl_fortran.dylib
   else
      $(error No OS detected, please have a look into the Makefile)
   endif

   OMP_FLAG=-fopenmp

   ifeq ($(DEBUG),true)
      DEBUG_FLAGS = -traceback -fpe0 -mcmodel=medium -shared-intel -check all -debug all
   else
      DEBUG_FLAGS =
   endif

   CFLAGS = -w -c -O3 -mp -assume byterecl -r8 $(OMP_FLAG) $(DEBUG_FLAG)
   MFLAGS = -module

   CDEFINED = true

endif

#####################MPIF909 COMPILER OPTIONS ETC########################

ifneq (,$(findstring mpif90,$(COMPILER)))

   F90 = COMPILER
   LD = COMPILER

   ifeq ($(UNAME), Linux)
      $(info OS detected: Linux)
      DIR_MOD_HDF5=$(LIB_HDF5)/include
      DIR_LIB_HDF5=$(LIB_HDF5)/lib                      \
                   $(LIB_HDF5)/lib/libhdf5_fortran.so   \
                   $(LIB_HDF5)/lib/libhdf5.so           \
                   $(LIB_HDF5)/lib/libhdf5_hl.so        \
                   $(LIB_HDF5)/lib/libhdf5hl_fortran.so
   else ifeq ($(UNAME), Darwin)
      $(info OS detected: MAC)
      DIR_MOD_HDF5=$(LIB_HDF5)/include
      DIR_LIB_HDF5=$(LIB_HDF5)/lib                         \
                   $(LIB_HDF5)/lib/libhdf5_fortran.dylib   \
                   $(LIB_HDF5)/lib/libhdf5.dylib           \
                   $(LIB_HDF5)/lib/libhdf5_hl.dylib        \
                   $(LIB_HDF5)/lib/libhdf5hl_fortran.dylib
   else
      $(error No OS detected, please have a look into the Makefile)
   endif

   OMP_FLAG=-fopenmp

   ifeq ($(DEBUG),true)
      DEBUG_FLAGS = -Og -Wall -Wextra -Warray-temporaries -Wconversion -pedantic-errors -fcheck=all -fbacktrace
   else
      DEBUG_FLAGS =
   endif

   CFLAGS = -ffree-line-length-none -c -O3 $(OMP_FLAG) $(DEBUG_FLAGS)
   MFLAGS = -J

   CDEFINED = true

endif

#########################################################################

ifneq ($(CDEFINED),true)
   ifneq ($(MAKECMDGOALS),clean)
      $(error compiler $(COMPILER) not defined)
   endif
endif


########################################################################

DIR_SRC = src
DIR_OBJ = objects/objects
DIR_MOD = modules/modules

DIR_SRC_MODEL = src_model
DIR_OBJ_MODEL = objects/objects_model
DIR_MOD_MODEL = modules/modules_model

DIR_SRC_SC3D = src_sc3d
DIR_OBJ_SC3D = objects/objects_sc3d
DIR_MOD_SC3D = modules/modules_sc3d

DIR_SRC_PHOTPROF = src_photprof
DIR_OBJ_PHOTPROF = objects/objects_photprof
DIR_MOD_PHOTPROF = modules/modules_photprof

DIR_SRC_SPEC = src_spec
DIR_OBJ_SPEC = objects/objects_spec
DIR_MOD_SPEC = modules/modules_spec

DIR_SRC_MODELSPEC = src_modelspec
DIR_OBJ_MODELSPEC = objects/objects_modelspec
DIR_MOD_MODELSPEC = modules/modules_modelspec

DIR_SRC_SPECVBIN = src_spec_vbin
DIR_OBJ_SPECVBIN = objects/objects_spec_vbin
DIR_MOD_SPECVBIN = modules/modules_spec_vbin

DIR_SRC_MODELSPECVBIN = src_modelspec_vbin
DIR_OBJ_MODELSPECVBIN = objects/objects_modelspec_vbin
DIR_MOD_MODELSPECVBIN = modules/modules_modelspec_vbin

DIR_SRC_OPAL = src_opal
DIR_OBJ_OPAL = objects/objects_opal
DIR_MOD_OPAL = modules/modules_opal

DIR_SRC_LTE = src_lte
DIR_OBJ_LTE = objects/objects_lte
DIR_MOD_LTE = modules/modules_lte


#module files
OBJSM_TYPE = $(DIR_OBJ)/mod_type.o

OBJSM_MODEL = $(DIR_OBJ)/mod_type.o                   \
              $(DIR_OBJ)/mod_sort.o                   \
              $(DIR_OBJ)/mod_grid.o                   \
              $(DIR_OBJ)/mod_integ1d.o                \
              $(DIR_OBJ)/mod_interp1d.o               \
              $(DIR_OBJ)/mod_interp2d.o               \
              $(DIR_OBJ)/mod_interp3d.o               \
              $(DIR_OBJ)/mod_opacities.o              \
              $(DIR_OBJ)/mod_iline.o                  \
              $(DIR_OBJ)/mod_amrvac_reader.o          \
              $(DIR_OBJ_LTE)/mod_lte.o                \
              $(DIR_OBJ_MODEL)/mod_model.o

OBJSM_PHOTPROF = $(DIR_OBJ)/mod_type.o                \
                 $(DIR_OBJ)/mod_grid.o                \
                 $(DIR_OBJ)/mod_interp1d.o            \
                 $(DIR_OBJ_PHOTPROF)/mod_photprof.o

OBJSM_SC3D = $(DIR_OBJ)/mod_type.o                    \
             $(DIR_OBJ)/mod_sort.o                    \
             $(DIR_OBJ)/mod_grid.o                    \
             $(DIR_OBJ)/mod_integ1d.o                 \
             $(DIR_OBJ)/mod_interp1d.o                \
             $(DIR_OBJ)/mod_interp2d.o                \
             $(DIR_OBJ)/mod_interp3d.o                \
             $(DIR_OBJ)/mod_opacities.o               \
             $(DIR_OBJ)/mod_iline.o                   \
             $(DIR_OBJ_LTE)/mod_lte.o                 \
             $(DIR_OBJ_SC3D)/mod_sc3d.o

OBJSM_SPEC = $(DIR_OBJ)/mod_type.o                    \
             $(DIR_OBJ)/mod_grid.o                    \
             $(DIR_OBJ)/mod_sort.o                    \
             $(DIR_OBJ)/mod_integ1d.o                 \
             $(DIR_OBJ)/mod_interp1d.o                \
             $(DIR_OBJ)/mod_interp2d.o                \
             $(DIR_OBJ)/mod_interp3d.o                \
             $(DIR_OBJ)/mod_opacities.o               \
             $(DIR_OBJ)/mod_iline.o                   \
             $(DIR_OBJ_LTE)/mod_lte.o                 \
             $(DIR_OBJ_PHOTPROF)/mod_photprof.o       \
             $(DIR_OBJ_SPEC)/mod_spec.o

OBJSM_MODELSPEC = $(DIR_OBJ)/mod_type.o               \
                  $(DIR_OBJ)/mod_sort.o               \
                  $(DIR_OBJ)/mod_grid.o               \
                  $(DIR_OBJ)/mod_integ1d.o            \
                  $(DIR_OBJ)/mod_interp1d.o           \
                  $(DIR_OBJ)/mod_interp2d.o           \
                  $(DIR_OBJ)/mod_interp3d.o           \
                  $(DIR_OBJ)/mod_opacities.o          \
                  $(DIR_OBJ)/mod_iline.o              \
                  $(DIR_OBJ_LTE)/mod_lte.o            \
                  $(DIR_OBJ_MODELSPEC)/mod_modelspec.o

OBJSM_SPECVBIN = $(DIR_OBJ)/mod_type.o                \
                 $(DIR_OBJ)/mod_sort.o                \
                 $(DIR_OBJ)/mod_grid.o                \
                 $(DIR_OBJ)/mod_integ1d.o             \
                 $(DIR_OBJ)/mod_interp1d.o            \
                 $(DIR_OBJ)/mod_interp2d.o            \
                 $(DIR_OBJ)/mod_interp3d.o            \
                 $(DIR_OBJ)/mod_opacities.o           \
                 $(DIR_OBJ)/mod_iline.o               \
                 $(DIR_OBJ_LTE)/mod_lte.o             \
                 $(DIR_OBJ_PHOTPROF)/mod_photprof.o   \
                 $(DIR_OBJ_SPECVBIN)/mod_spec.o

OBJSM_MODELSPECVBIN = $(DIR_OBJ)/mod_type.o           \
                      $(DIR_OBJ)/mod_sort.o           \
                      $(DIR_OBJ)/mod_grid.o           \
                      $(DIR_OBJ)/mod_integ1d.o        \
                      $(DIR_OBJ)/mod_interp1d.o       \
                      $(DIR_OBJ)/mod_interp2d.o       \
                      $(DIR_OBJ)/mod_interp3d.o       \
                      $(DIR_OBJ)/mod_opacities.o      \
                      $(DIR_OBJ)/mod_iline.o          \
                      $(DIR_OBJ_LTE)/mod_lte.o        \
                      $(DIR_OBJ_MODELSPECVBIN)/mod_modelspec.o

OBJSM_OPAL = $(DIR_OBJ)/mod_type.o \
             $(DIR_OBJ_OPAL)/mod_opal.o

OBJSM_LTE = $(DIR_OBJ)/mod_type.o \
            $(DIR_OBJ_OPAL)/mod_lte.o

########################################################################

#for creating atmospheric model as input to sc3d.eo
OBJS1_MODEL = $(DIR_OBJ)/math.o                         \
              $(DIR_OBJ)/eispackEV.o                    \
              $(DIR_OBJ)/model_laws.o                   \
              $(DIR_OBJ)/info_region.o                  \
              $(DIR_OBJ_MODEL)/model1d.o                \
              $(DIR_OBJ_MODEL)/model2d.o                \
              $(DIR_OBJ_MODEL)/model3d.o                \
              $(DIR_OBJ_MODEL)/model_output.o
OBJS2_MODEL = $(OBJS1_MODEL) $(OBJSM_MODEL) $(DIR_OBJ_MODEL)/model.o
#
#for main program sc3d.eo
OBJS1_SC3D = $(DIR_OBJ)/math.o                          \
             $(DIR_OBJ)/eispackEV.o                     \
             $(DIR_OBJ)/formal_sc.o                     \
             $(DIR_OBJ)/ng_extrapol.o                   \
             $(DIR_OBJ)/sobolev1d.o                     \
             $(DIR_OBJ)/sobolev3d.o                     \
             $(DIR_OBJ)/sparse.o                        \
             $(DIR_OBJ)/lebedev.o                       \
             $(DIR_OBJ)/info_region.o                   \
             $(DIR_OBJ)/model_laws.o                    \
             $(DIR_OBJ_SC3D)/anglenodes.o               \
             $(DIR_OBJ_SC3D)/formal_sc2d.o                   \
             $(DIR_OBJ_SC3D)/formal_sc3d.o                   \
             $(DIR_OBJ_SC3D)/formal_fvm2d.o                  \
             $(DIR_OBJ_SC3D)/formal_fvm3d.o                  \
             $(DIR_OBJ_SC3D)/output_sc3d.o                   \
             $(DIR_OBJ_SC3D)/mint_fvm3d.o                    \
             $(DIR_OBJ_SC3D)/mint_sc3d.o                     \
             $(DIR_OBJ_SC3D)/scont_new3d.o                   \
             $(DIR_OBJ_SC3D)/sline_new3d.o                   \
             $(DIR_OBJ_SC3D)/grid_spatial.o                  \
             $(DIR_OBJ_SC3D)/grid_delxyz.o                   \
             $(DIR_OBJ_SC3D)/model_sc3d.o                    \
             $(DIR_OBJ_SC3D)/benchmark01.o                   \
             $(DIR_OBJ_SC3D)/benchmark02.o                   \
             $(DIR_OBJ_SC3D)/benchmark03.o                   \
             $(DIR_OBJ_SC3D)/benchmark04.o                   \
             $(DIR_OBJ_SC3D)/benchmark05.o                   \
             $(DIR_OBJ_SC3D)/benchmark06.o                   \
             $(DIR_OBJ_SC3D)/benchmark07.o                   \
             $(DIR_OBJ_SC3D)/benchmark08.o                   \
             $(DIR_OBJ_SC3D)/benchmark09.o                   \
             $(DIR_OBJ_SC3D)/benchmark10.o                   \
             $(DIR_OBJ_SC3D)/benchmark11.o                   \
             $(DIR_OBJ_SC3D)/benchmark12.o                   \
             $(DIR_OBJ_SC3D)/benchmark13.o                   \
             $(DIR_OBJ_SC3D)/benchmark14.o                   \
             $(DIR_OBJ_SC3D)/tests.o                         
OBJS2_SC3D = $(OBJS1_SC3D) $(OBJSM_SC3D) $(DIR_OBJ_SC3D)/sc3d.o


OBJS1_PHOTPROF =  $(DIR_OBJ)/math.o                    \
                  $(DIR_OBJ)/eispackEV.o

OBJS2_PHOTPROF = $(OBJS1_PHOTPROF) $(OBJSM_PHOTPROF) $(DIR_OBJ_PHOTPROF)/phot_profile.o

#formal integral single star
OBJS1_SPEC =  $(DIR_OBJ)/info_region.o             \
              $(DIR_OBJ)/math.o                    \
              $(DIR_OBJ)/eispackEV.o               \
              $(DIR_OBJ)/model_laws.o              \
              $(DIR_OBJ_PHOTPROF)/phot_profile.o            \
              $(DIR_OBJ_SPEC)/spec_output.o             \
              $(DIR_OBJ_SPEC)/spec_input.o              \
              $(DIR_OBJ_SPEC)/spec_tests.o              \
              $(DIR_OBJ_SPEC)/formal_lc.o

OBJS2_SPEC = $(OBJS1_SPEC) $(OBJSM_SPEC) $(DIR_OBJ_SPEC)/spec.o


#for creating atmospheric model that can be read in by spec.eo
OBJS1_MODELSPEC =  $(DIR_OBJ)/model_laws.o            \
                   $(DIR_OBJ)/sobolev1d.o             \
                   $(DIR_OBJ)/info_region.o            \
                   $(DIR_OBJ)/math.o                    \
                   $(DIR_OBJ)/eispackEV.o              \
                   $(DIR_OBJ_OPAL)/mod_opal.o

OBJS2_MODELSPEC = $(OBJS1_MODELSPEC) $(OBJSM_MODELSPEC) $(DIR_OBJ_MODELSPEC)/modelspec.o


#formal integral binary system
OBJS1_SPECVBIN =  $(DIR_OBJ)/info_region.o             \
                  $(DIR_OBJ)/math.o                    \
                  $(DIR_OBJ)/eispackEV.o               \
                  $(DIR_OBJ)/geompack2.o               \
                  $(DIR_OBJ)/model_laws.o              \
                  $(DIR_OBJ_PHOTPROF)/phot_profile.o            \
                  $(DIR_OBJ_SPECVBIN)/spec_triangulation.o      \
                  $(DIR_OBJ_SPECVBIN)/spec_transformations.o     \
                  $(DIR_OBJ_SPECVBIN)/spec_output.o             \
                  $(DIR_OBJ_SPECVBIN)/spec_input.o              \
                  $(DIR_OBJ_SPECVBIN)/spec_tests.o              \
                  $(DIR_OBJ_SPECVBIN)/spec_setupray3d.o         \
                  $(DIR_OBJ_SPECVBIN)/formal_lc.o

OBJS2_SPECVBIN = $(OBJS1_SPECVBIN) $(OBJSM_SPECVBIN) $(DIR_OBJ_SPECVBIN)/spec.o

#for creating atmospheric model that can be read in by spec_vbin.eo
OBJS1_MODELSPECVBIN =  $(DIR_OBJ)/model_laws.o            \
                       $(DIR_OBJ)/sobolev1d.o             \
                       $(DIR_OBJ)/info_region.o            \
                       $(DIR_OBJ)/math.o                    \
                       $(DIR_OBJ)/eispackEV.o

OBJS2_MODELSPECVBIN = $(OBJS1_MODELSPECVBIN) $(OBJSM_MODELSPECVBIN) $(DIR_OBJ_MODELSPECVBIN)/modelspec.o



########################################################################

all: model sc3d spec modelspec spec_vbin modelspec_vbin


CPDIR = ../florian3d

copy:
	mkdir $(CPDIR)
	mkdir $(CPDIR)/outputFILES
	mkdir $(CPDIR)/outputFILES_TEST
	mkdir $(CPDIR)/outputFILES_TEMP
	mkdir $(CPDIR)/inputFILES
	mkdir $(CPDIR)/plotFILES
	mkdir $(CPDIR)/modules
	mkdir $(CPDIR)/modules/modules
	mkdir $(CPDIR)/modules/modules_model
	mkdir $(CPDIR)/modules/modules_modelspec
	mkdir $(CPDIR)/modules/modules_modelspec_vbin
	mkdir $(CPDIR)/modules/modules_spec
	mkdir $(CPDIR)/modules/modules_spec_vbin
	mkdir $(CPDIR)/modules/modules_opal
	mkdir $(CPDIR)/modules/modules_lte
	mkdir $(CPDIR)/modules/modules_photprof
	mkdir $(CPDIR)/modules/modules_sc3d
	mkdir $(CPDIR)/objects
	mkdir $(CPDIR)/objects/objects
	mkdir $(CPDIR)/objects/objects_model
	mkdir $(CPDIR)/objects/objects_modelspec
	mkdir $(CPDIR)/objects/objects_modelspec_vbin
	mkdir $(CPDIR)/objects/objects_spec
	mkdir $(CPDIR)/objects/objects_spec_vbin
	mkdir $(CPDIR)/objects/objects_opal
	mkdir $(CPDIR)/objects/objects_lte
	mkdir $(CPDIR)/objects/objects_photprof
	mkdir $(CPDIR)/objects/objects_sc3d

	cp -r outputFILES_TEST/*.pro $(CPDIR)/outputFILES_TEST/
	cp -r outputFILES_TEST/*.py $(CPDIR)/outputFILES_TEST/
	cp -r plotFILES/*.pro $(CPDIR)/plotFILES/
	cp -r plotFILES/*.py $(CPDIR)/plotFILES/
#	cp -r phot_flux $(CPDIR)/
	cp -r opal_tables $(CPDIR)/
	cp -r lte_tables $(CPDIR)/
	cp in_sc3d $(CPDIR)/
	cp in_spec $(CPDIR)/
	cp in_spec_vbin $(CPDIR)/
	cp in_modelspec $(CPDIR)/
	cp in_modelspec_vbin $(CPDIR)/
	cp indat_sc3d.nml $(CPDIR)/
	cp indat_spec.nml $(CPDIR)/
	cp indat_spec_vbin.nml $(CPDIR)/
	cp indat_modelspec.nml $(CPDIR)/
	cp indat_modelspec_vbin.nml $(CPDIR)/
	cp Makefile $(CPDIR)/
	cp -r src* $(CPDIR)/
	cp *.txt $(CPDIR)/
#
########################################################################
#
model: $(OBJS2_MODEL)
#	@echo $(OBJS2_MODEL)
	 $(LD) $(OBJS2_MODEL) $(OMP_FLAG) -L $(DIR_LIB_HDF5) -o model.eo

sc3d: $(OBJS2_SC3D)
	 $(LD) $(OBJS2_SC3D) $(OMP_FLAG) -L $(DIR_LIB_HDF5) -o sc3d.eo

spec: $(OBJS2_SPEC)
	 $(LD) $(OBJS2_SPEC) $(OMP_FLAG) -L $(DIR_LIB_HDF5) -o spec.eo

modelspec: $(OBJS2_MODELSPEC)
	 $(LD) $(OBJS2_MODELSPEC) $(OMP_FLAG) -L $(DIR_LIB_HDF5) -o modelspec.eo

spec_vbin: $(OBJS2_SPECVBIN)
	 $(LD) $(OBJS2_SPECVBIN) $(OMP_FLAG) -L $(DIR_LIB_HDF5) -o spec_vbin.eo

modelspec_vbin: $(OBJS2_MODELSPECVBIN)
	 $(LD) $(OBJS2_MODELSPECVBIN) $(OMP_FLAG) -L $(DIR_LIB_HDF5) -o modelspec_vbin.eo
#
########################################################################
#
$(DIR_OBJ)/mod_interp1d.o: $(DIR_SRC)/mod_interp1d.f90 \
                           $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $(MFLAGS) $(DIR_MOD) $< -o $(DIR_OBJ)/mod_interp1d.o

$(DIR_OBJ)/mod_interp2d.o: $(DIR_SRC)/mod_interp2d.f90 \
                           $(DIR_OBJ)/mod_interp1d.o   \
                           $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $(MFLAGS) $(DIR_MOD) $< -o $(DIR_OBJ)/mod_interp2d.o

$(DIR_OBJ)/mod_interp3d.o: $(DIR_SRC)/mod_interp3d.f90  \
                           $(DIR_OBJ)/mod_interp1d.o \
                           $(DIR_OBJ)/info_region.o \
                           $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $(MFLAGS) $(DIR_MOD) $< -o $(DIR_OBJ)/mod_interp3d.o

$(DIR_OBJ)/mod_opacities.o: $(DIR_SRC)/mod_opacities.f90  \
                            $(DIR_OBJ)/mod_interp1d.o \
                            $(DIR_OBJ)/mod_interp2d.o \
                            $(DIR_OBJ)/info_region.o \
                            $(DIR_OBJ)/mod_iline.o \
                            $(DIR_OBJ_LTE)/mod_lte.o \
                            $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) -I $(DIR_MOD_LTE) $(MFLAGS) $(DIR_MOD) $< -o $(DIR_OBJ)/mod_opacities.o

$(DIR_OBJ)/info_region.o: $(DIR_SRC)/info_region.f90  \
                          $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $< -o $(DIR_OBJ)/info_region.o

$(DIR_OBJ)/mod_integ1d.o: $(DIR_SRC)/mod_integ1d.f90 \
                      $(DIR_OBJ)/mod_interp2d.o \
                      $(DIR_OBJ)/math.o \
                      $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $(MFLAGS) $(DIR_MOD) $< -o $(DIR_OBJ)/mod_integ1d.o

$(DIR_OBJ)/eispackEV.o: $(DIR_SRC)/eispackEV.f90 \
                        $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $< -o $(DIR_OBJ)/eispackEV.o

$(DIR_OBJ)/geompack2.o: $(DIR_SRC)/geompack2.f90 \
                        $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $< -o $(DIR_OBJ)/geompack2.o

$(DIR_OBJ)/lebedev.o: $(DIR_SRC)/lebedev.f90 \
                      $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $< -o $(DIR_OBJ)/lebedev.o

$(DIR_OBJ)/lapackINV.o: $(DIR_SRC)/lapackINV.f
		$(F90) $(CFLAGS) $< -o $(DIR_OBJ)/lapackINV.o

$(DIR_OBJ)/math.o: $(DIR_SRC)/math.f90 \
                   $(DIR_OBJ)/eispackEV.o \
                   $(DIR_OBJ)/mod_sort.o \
                   $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $< -o $(DIR_OBJ)/math.o

$(DIR_OBJ)/mod_sort.o: $(DIR_SRC)/mod_sort.f90 \
                       $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD)  $(MFLAGS) $(DIR_MOD) $< -o $(DIR_OBJ)/mod_sort.o

$(DIR_OBJ)/mod_grid.o: $(DIR_SRC)/mod_grid.f90 \
                       $(DIR_OBJ)/mod_sort.o \
                       $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD)  $(MFLAGS) $(DIR_MOD) $< -o $(DIR_OBJ)/mod_grid.o

$(DIR_OBJ)/mod_amrvac_reader.o: $(DIR_SRC)/mod_amrvac_reader.f90 \
                                $(DIR_OBJ)/mod_sort.o \
                                $(DIR_OBJ)/mod_grid.o \
                                $(DIR_OBJ)/mod_interp1d.o \
                                $(DIR_OBJ)/mod_interp3d.o \
                                $(OBJSM_TYPE)
		$(F90) $(CFLAGS) $(MFLAGS) $(DIR_MOD) -I $(DIR_MOD) -I $(DIR_MOD_HDF5) $< -o $(DIR_OBJ)/mod_amrvac_reader.o

$(DIR_OBJ)/sobolev1d.o: $(DIR_SRC)/sobolev1d.f90 \
                        $(DIR_OBJ)/mod_integ1d.o \
                        $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $< -o $(DIR_OBJ)/sobolev1d.o

$(DIR_OBJ)/sobolev3d.o: $(DIR_SRC)/sobolev3d.f90 \
                        $(DIR_OBJ)/sobolev1d.o \
                        $(DIR_OBJ)/math.o \
                        $(DIR_OBJ)/mod_integ1d.o \
                        $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $< -o $(DIR_OBJ)/sobolev3d.o

$(DIR_OBJ)/ng_extrapol.o: $(DIR_SRC)/ng_extrapol.f90 \
                          $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $< -o $(DIR_OBJ)/ng_extrapol.o

$(DIR_OBJ)/formal_sc.o: $(DIR_SRC)/formal_sc.f90 \
                        $(DIR_OBJ)/mod_integ1d.o \
                        $(DIR_OBJ)/mod_interp1d.o \
                        $(DIR_OBJ)/math.o \
                        $(DIR_OBJ)/mod_type.o \
                        $(DIR_OBJ)/mod_interp2d.o \
                        $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $< -o $(DIR_OBJ)/formal_sc.o

$(DIR_OBJ)/model_laws.o: $(DIR_SRC)/model_laws.f90  \
                         $(DIR_OBJ)/mod_interp1d.o \
                         $(DIR_OBJ)/mod_interp2d.o \
                         $(DIR_OBJ)/mod_opacities.o \
                         $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $< -o $(DIR_OBJ)/model_laws.o

$(DIR_OBJ)/mod_iline.o: $(DIR_SRC)/mod_iline.f90 \
                        $(OBJSM_TYPE)
		$(F90) $(CFLAGS) $(MFLAGS) $(DIR_MOD) -I $(DIR_MOD) $< -o $(DIR_OBJ)/mod_iline.o

$(DIR_OBJ)/sparse.o: $(DIR_SRC)/sparse.f90 \
                     $(DIR_OBJ)/ng_extrapol.o \
                     $(DIR_OBJ)/math.o \
                     $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $< -o $(DIR_OBJ)/sparse.o

$(DIR_OBJ)/mod_type.o: $(DIR_SRC)/mod_type.f90
		$(F90) $(CFLAGS) -I $(DIR_MOD) $(MFLAGS) $(DIR_MOD) $< -o $(DIR_OBJ)/mod_type.o

########################################################################

$(DIR_OBJ_OPAL)/mod_opal.o: $(DIR_SRC_OPAL)/mod_opal.f90 \
                            $(OBJSM_TYPE)
		$(F90) $(CFLAGS) $(MFLAGS) $(DIR_MOD_OPAL) -I $(DIR_MOD) $< -o $(DIR_OBJ_OPAL)/mod_opal.o

########################################################################

$(DIR_OBJ_LTE)/mod_lte.o: $(DIR_SRC_LTE)/mod_lte.f90 \
                          $(DIR_OBJ/mod_iline.o \
                          $(OBJSM_TYPE)
		$(F90) $(CFLAGS) $(MFLAGS) $(DIR_MOD_LTE) -I $(DIR_MOD) $< -o $(DIR_OBJ_LTE)/mod_lte.o


########################################################################

$(DIR_OBJ_MODEL)/mod_model.o: $(DIR_SRC_MODEL)/mod_model.f90 \
                              $(OBJSM_TYPE)
		$(F90) $(CFLAGS) $(MFLAGS) $(DIR_MOD_MODEL) -I $(DIR_MOD) -I $(DIR_MOD_HDF5) -I $(DIR_MOD_LTE) $< -o $(DIR_OBJ_MODEL)/mod_model.o

$(DIR_OBJ_MODEL)/model1d.o: $(DIR_SRC_MODEL)/model1d.f90 \
                            $(DIR_OBJ)/model_laws.o \
                            $(OBJSM_MODEL)
		$(F90) $(CFLAGS) -I $(DIR_MOD_MODEL) -I $(DIR_MOD) $< -o $(DIR_OBJ_MODEL)/model1d.o

$(DIR_OBJ_MODEL)/model2d.o: $(DIR_SRC_MODEL)/model2d.f90 \
                            $(DIR_OBJ)/model_laws.o \
                            $(OBJSM_MODEL)
		$(F90) $(CFLAGS) -I $(DIR_MOD_MODEL) -I $(DIR_MOD) -I $(DIR_MOD_HDF5) $< -o $(DIR_OBJ_MODEL)/model2d.o

$(DIR_OBJ_MODEL)/model3d.o: $(DIR_SRC_MODEL)/model3d.f90 \
                            $(OBJSM_MODEL)
		$(F90) $(CFLAGS) -I $(DIR_MOD_MODEL) -I $(DIR_MOD) -I $(DIR_MOD_HDF5) $< -o $(DIR_OBJ_MODEL)/model3d.o

$(DIR_OBJ_MODEL)/model_output.o: $(DIR_SRC_MODEL)/model_output.f90 \
                                 $(OBJSM_MODEL)
		$(F90) $(CFLAGS) -I $(DIR_MOD_MODEL) -I $(DIR_MOD) -I $(DIR_MOD_HDF5) $< -o $(DIR_OBJ_MODEL)/model_output.o

$(DIR_OBJ_MODEL)/model.o: $(DIR_SRC_MODEL)/model.f90 \
                          $(OBJSM_MODEL) \
                          $(OBJS1_MODEL)
		$(F90) $(CFLAGS) -I $(DIR_MOD_MODEL) -I $(DIR_MOD) $< -o $(DIR_OBJ_MODEL)/model.o

########################################################################

$(DIR_OBJ_SC3D)/mod_sc3d.o: $(DIR_SRC_SC3D)/mod_sc3d.f90 \
                            $(OBJSM_TYPE)
		$(F90) $(CFLAGS) $(MFLAGS) $(DIR_MOD_SC3D) -I $(DIR_MOD) -I $(DIR_MOD_LTE) $< -o $(DIR_OBJ_SC3D)/mod_sc3d.o

$(DIR_OBJ_SC3D)/anglenodes.o: $(DIR_SRC_SC3D)/anglenodes.f90 \
                              $(DIR_OBJ)/math.o \
                              $(DIR_OBJ)/lebedev.o \
                              $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/anglenodes.o

$(DIR_OBJ_SC3D)/formal_sc2d.o: $(DIR_SRC_SC3D)/formal_sc2d.f90 \
                               $(DIR_OBJ)/math.o \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/formal_sc2d.o

$(DIR_OBJ_SC3D)/formal_sc3d.o: $(DIR_SRC_SC3D)/formal_sc3d.f90 \
                               $(DIR_OBJ)/math.o \
                               $(DIR_OBJ)/model_laws.o \
                               $(DIR_OBJ)/mod_opacities.o \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/formal_sc3d.o

$(DIR_OBJ_SC3D)/formal_fvm2d.o: $(DIR_SRC_SC3D)/formal_fvm2d.f90 \
                                $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/formal_fvm2d.o

$(DIR_OBJ_SC3D)/formal_fvm3d.o: $(DIR_SRC_SC3D)/formal_fvm3d.f90 \
                                $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/formal_fvm3d.o

$(DIR_OBJ_SC3D)/output_sc3d.o: $(DIR_SRC_SC3D)/output_sc3d.f90 \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) -I $(DIR_MOD_HDF5) $< -o $(DIR_OBJ_SC3D)/output_sc3d.o

$(DIR_OBJ_SC3D)/mint_sc3d.o: $(DIR_SRC_SC3D)/mint_sc3d.f90 \
                             $(DIR_OBJ_SC3D)/formal_sc3d.o \
                             $(DIR_OBJ_SC3D)/scont_new3d.o \
                             $(DIR_OBJ_SC3D)/sline_new3d.o \
                             $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/mint_sc3d.o

$(DIR_OBJ_SC3D)/mint_fvm3d.o: $(DIR_SRC_SC3D)/mint_fvm3d.f90 \
                              $(DIR_OBJ_SC3D)/formal_fvm3d.o \
                              $(DIR_OBJ_SC3D)/scont_new3d.o \
                              $(DIR_OBJ_SC3D)/sline_new3d.o \
                              $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/mint_fvm3d.o

$(DIR_OBJ_SC3D)/scont_new3d.o: $(DIR_SRC_SC3D)/scont_new3d.f90 \
                               $(DIR_OBJ)/math.o \
                               $(DIR_OBJ)/sparse.o \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/scont_new3d.o

$(DIR_OBJ_SC3D)/sline_new3d.o: $(DIR_SRC_SC3D)/sline_new3d.f90 \
                               $(DIR_OBJ)/math.o \
                               $(DIR_OBJ)/sparse.o \
                               $(DIR_OBJ)/mod_interp2d.o \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/sline_new3d.o

$(DIR_OBJ_SC3D)/grid_delxyz.o: $(DIR_SRC_SC3D)/grid_delxyz.f90 \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/grid_delxyz.o

$(DIR_OBJ_SC3D)/grid_spatial.o: $(DIR_SRC_SC3D)/grid_spatial.f90 \
                                $(DIR_OBJ)/mod_sort.o \
                                $(DIR_OBJ)/mod_grid.o \
                                $(DIR_OBJ)/mod_interp1d.o \
                                $(DIR_OBJ)/model_laws.o \
                                $(DIR_OBJ_SC3D)/model_sc3d.o \
                                $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) -I $(DIR_MOD_HDF5) $< -o $(DIR_OBJ_SC3D)/grid_spatial.o

$(DIR_OBJ_SC3D)/model_sc3d.o: $(DIR_SRC_SC3D)/model_sc3d.f90 \
                              $(DIR_OBJ)/math.o \
                              $(DIR_OBJ)/mod_interp2d.o \
                              $(DIR_OBJ)/mod_opacities.o \
                              $(DIR_OBJ)/info_region.o \
                              $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) -I $(DIR_MOD_HDF5) $< -o $(DIR_OBJ_SC3D)/model_sc3d.o

$(DIR_OBJ_SC3D)/tests.o: $(DIR_SRC_SC3D)/tests.f90 \
                         $(DIR_OBJ)/math.o \
                         $(DIR_OBJ)/mod_interp1d.o \
                         $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) -I $(DIR_MOD_HDF5) $< -o $(DIR_OBJ_SC3D)/tests.o

$(DIR_OBJ_SC3D)/benchmark01.o: $(DIR_SRC_SC3D)/benchmark01.f90 \
                               $(DIR_OBJ)/math.o \
                               $(DIR_OBJ)/mod_interp2d.o \
                               $(DIR_OBJ)/formal_sc.o \
                               $(DIR_OBJ_SC3D)/formal_sc2d.o \
                               $(DIR_OBJ_SC3D)/formal_fvm2d.o \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/benchmark01.o

$(DIR_OBJ_SC3D)/benchmark02.o: $(DIR_SRC_SC3D)/benchmark02.f90 \
                               $(DIR_OBJ)/formal_sc.o \
                               $(DIR_OBJ_SC3D)/formal_sc2d.o \
                               $(DIR_OBJ_SC3D)/formal_fvm2d.o \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/benchmark02.o

$(DIR_OBJ_SC3D)/benchmark03.o: $(DIR_SRC_SC3D)/benchmark03.f90 \
                               $(DIR_OBJ)/formal_sc.o \
                               $(DIR_OBJ_SC3D)/formal_sc2d.o \
                               $(DIR_OBJ_SC3D)/formal_fvm2d.o \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/benchmark03.o

$(DIR_OBJ_SC3D)/benchmark04.o: $(DIR_SRC_SC3D)/benchmark04.f90 \
                               $(DIR_OBJ)/formal_sc.o \
                               $(DIR_OBJ_SC3D)/formal_sc2d.o \
                               $(DIR_OBJ_SC3D)/formal_fvm2d.o \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/benchmark04.o

$(DIR_OBJ_SC3D)/benchmark05.o: $(DIR_SRC_SC3D)/benchmark05.f90 \
                               $(DIR_OBJ)/math.o \
                               $(DIR_OBJ)/mod_interp1d.o \
                               $(DIR_OBJ)/formal_sc.o \
                               $(DIR_OBJ_SC3D)/formal_sc2d.o \
                               $(DIR_OBJ_SC3D)/formal_fvm2d.o \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/benchmark05.o

$(DIR_OBJ_SC3D)/benchmark06.o: $(DIR_SRC_SC3D)/benchmark06.f90 \
                               $(DIR_OBJ)/mod_interp1d.o \
                               $(DIR_OBJ)/formal_sc.o \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/benchmark06.o

$(DIR_OBJ_SC3D)/benchmark07.o: $(DIR_SRC_SC3D)/benchmark07.f90 \
                               $(DIR_OBJ)/math.o \
                               $(DIR_OBJ)/mod_interp1d.o \
                               $(DIR_OBJ)/sobolev1d.o \
                               $(DIR_OBJ)/sobolev3d.o \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/benchmark07.o

$(DIR_OBJ_SC3D)/benchmark08.o: $(DIR_SRC_SC3D)/benchmark08.f90 \
                               $(DIR_OBJ)/mod_interp1d.o \
                               $(DIR_OBJ)/formal_sc.o \
                               $(DIR_OBJ)/sobolev1d.o \
                               $(DIR_OBJ)/sobolev3d.o \
                               $(DIR_OBJ_SC3D)/formal_sc2d.o \
                               $(DIR_OBJ_SC3D)/formal_fvm2d.o \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/benchmark08.o

$(DIR_OBJ_SC3D)/benchmark09.o: $(DIR_SRC_SC3D)/benchmark09.f90 \
                               $(DIR_OBJ)/sobolev3d.o \
                               $(DIR_OBJ_SC3D)/formal_sc2d.o \
                               $(DIR_OBJ_SC3D)/formal_fvm2d.o \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/benchmark09.o

$(DIR_OBJ_SC3D)/benchmark10.o: $(DIR_SRC_SC3D)/benchmark10.f90 \
                               $(DIR_OBJ)/mod_interp2d.o \
                               $(DIR_OBJ)/formal_sc.o \
                               $(DIR_OBJ_SC3D)/formal_sc3d.o \
                               $(DIR_OBJ_SC3D)/formal_fvm3d.o \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/benchmark10.o

$(DIR_OBJ_SC3D)/benchmark11.o: $(DIR_SRC_SC3D)/benchmark11.f90 \
                               $(DIR_OBJ_SC3D)/mint_sc3d.o \
                               $(DIR_OBJ_SC3D)/mint_fvm3d.o \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/benchmark11.o

$(DIR_OBJ_SC3D)/benchmark12.o: $(DIR_SRC_SC3D)/benchmark12.f90 \
                               $(DIR_OBJ)/math.o \
                               $(DIR_OBJ)/ng_extrapol.o \
                               $(DIR_OBJ_SC3D)/scont_new3d.o \
                               $(DIR_OBJ_SC3D)/mint_sc3d.o \
                               $(DIR_OBJ_SC3D)/mint_fvm3d.o \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/benchmark12.o

$(DIR_OBJ_SC3D)/benchmark13.o: $(DIR_SRC_SC3D)/benchmark13.f90 \
                               $(DIR_OBJ)/math.o \
                               $(DIR_OBJ)/ng_extrapol.o \
                               $(DIR_OBJ_SC3D)/sline_new3d.o \
                               $(DIR_OBJ_SC3D)/mint_sc3d.o \
                               $(DIR_OBJ_SC3D)/mint_fvm3d.o \
                               $(DIR_OBJ)/sobolev3d.o \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/benchmark13.o

$(DIR_OBJ_SC3D)/benchmark14.o: $(DIR_SRC_SC3D)/benchmark14.f90 \
                               $(DIR_OBJ)/math.o \
                               $(DIR_OBJ)/ng_extrapol.o \
                               $(DIR_OBJ_SC3D)/sline_new3d.o \
                               $(DIR_OBJ_SC3D)/mint_sc3d.o \
                               $(DIR_OBJ_SC3D)/mint_fvm3d.o \
                               $(DIR_OBJ)/sobolev3d.o \
                               $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/benchmark14.o

$(DIR_OBJ_SC3D)/sc3d.o: $(DIR_SRC_SC3D)/sc3d.f90 \
                        $(OBJSM_SC3D) \
                        $(OBJS1_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/sc3d.o


########################################################################

$(DIR_OBJ_PHOTPROF)/mod_photprof.o: $(DIR_SRC_PHOTPROF)/mod_photprof.f90 
		$(F90) $(CFLAGS) $(MFLAGS) $(DIR_MOD_PHOTPROF) -I $(DIR_MOD) $< -o $(DIR_OBJ_PHOTPROF)/mod_photprof.o

$(DIR_OBJ_PHOTPROF)/phot_profile.o: $(DIR_SRC_PHOTPROF)/phot_profile.f90 \
                                    $(OBJSM_PHOTPROF) \
                                    $(OBJS1_PHOTPROF)
		$(F90) $(CFLAGS) -I $(DIR_MOD_PHOTPROF) -I $(DIR_MOD) $< -o $(DIR_OBJ_PHOTPROF)/phot_profile.o

########################################################################

$(DIR_OBJ_SPEC)/formal_lc.o: $(DIR_SRC_SPEC)/formal_lc.f90 \
                              $(OBJSM_SPEC) 
		$(F90) $(CFLAGS) -I $(DIR_MOD_SPEC) -I $(DIR_MOD) $< -o $(DIR_OBJ_SPEC)/formal_lc.o

$(DIR_OBJ_SPEC)/spec_tests.o: $(DIR_SRC_SPEC)/spec_tests.f90 \
                              $(DIR_OBJ_SPEC)/formal_lc.o \
                              $(OBJSM_SPEC) 
		$(F90) $(CFLAGS) -I $(DIR_MOD_SPEC) -I $(DIR_MOD) $< -o $(DIR_OBJ_SPEC)/spec_tests.o

$(DIR_OBJ_SPEC)/spec_input.o: $(DIR_SRC_SPEC)/spec_input.f90 \
                              $(DIR_OBJ)/model_laws.o \
                              $(DIR_OBJ)/math.o \
                              $(DIR_OBJ)/info_region.o \
                              $(OBJSM_SPEC)
		$(F90) $(CFLAGS) -I $(DIR_MOD) -I $(DIR_MOD_SPEC) -I $(DIR_MOD_HDF5) $< -o $(DIR_OBJ_SPEC)/spec_input.o

$(DIR_OBJ_SPEC)/spec_output.o: $(DIR_SRC_SPEC)/spec_output.f90 \
                               $(OBJSM_SPEC)
		$(F90) $(CFLAGS) -I $(DIR_MOD) -I $(DIR_MOD_SPEC) -I $(DIR_MOD_HDF5) $< -o $(DIR_OBJ_SPEC)/spec_output.o

$(DIR_OBJ_SPEC)/spec.o: $(DIR_SRC_SPEC)/spec.f90 \
                        $(OBJSM_SPEC) \
                        $(OBJS1_SPEC)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SPEC) -I $(DIR_MOD) $< -o $(DIR_OBJ_SPEC)/spec.o

$(DIR_OBJ_SPEC)/mod_spec.o: $(DIR_SRC_SPEC)/mod_spec.f90 \
                            $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $(MFLAGS) $(DIR_MOD_SPEC) $< -o $(DIR_OBJ_SPEC)/mod_spec.o

########################################################################

$(DIR_OBJ_MODELSPEC)/modelspec.o: $(DIR_SRC_MODELSPEC)/modelspec.f90 \
                                  $(OBJSM_MODELSPEC) \
                                  $(OBJS1_MODELSPEC)
		$(F90) $(CFLAGS) -I $(DIR_MOD) -I $(DIR_MOD_MODELSPEC) -I $(DIR_MOD_HDF5) -I $(DIR_MOD_OPAL) -I $(DIR_MOD_LTE) $< -o $(DIR_OBJ_MODELSPEC)/modelspec.o

$(DIR_OBJ_MODELSPEC)/mod_modelspec.o: $(DIR_SRC_MODELSPEC)/mod_modelspec.f90 \
                                      $(OBJSM_TYPE)
		$(F90) $(CFLAGS) $(MFLAGS) $(DIR_MOD_MODELSPEC) -I $(DIR_MOD) $< -o $(DIR_OBJ_MODELSPEC)/mod_modelspec.o

########################################################################

$(DIR_OBJ_SPECVBIN)/formal_lc.o: $(DIR_SRC_SPECVBIN)/formal_lc.f90 \
                                 $(OBJSM_SPECVBIN) 
		$(F90) $(CFLAGS) -I $(DIR_MOD_SPECVBIN) -I $(DIR_MOD) $< -o $(DIR_OBJ_SPECVBIN)/formal_lc.o

$(DIR_OBJ_SPECVBIN)/spec_tests.o: $(DIR_SRC_SPECVBIN)/spec_tests.f90 \
                                  $(DIR_OBJ_SPECVBIN)/formal_lc.o \
                                  $(OBJSM_SPECVBIN) 
		$(F90) $(CFLAGS) -I $(DIR_MOD_SPECVBIN) -I $(DIR_MOD) $< -o $(DIR_OBJ_SPECVBIN)/spec_tests.o

$(DIR_OBJ_SPECVBIN)/spec_triangulation.o: $(DIR_SRC_SPECVBIN)/spec_triangulation.f90 \
                                          $(DIR_OBJ)/geompack2.o \
                                          $(OBJSM_SPECVBIN) 
		$(F90) $(CFLAGS) -I $(DIR_MOD_SPECVBIN) -I $(DIR_MOD) $< -o $(DIR_OBJ_SPECVBIN)/spec_triangulation.o

$(DIR_OBJ_SPECVBIN)/spec_transformations.o: $(DIR_SRC_SPECVBIN)/spec_transformations.f90 \
                                            $(OBJSM_SPECVBIN) 
		$(F90) $(CFLAGS) -I $(DIR_MOD_SPECVBIN) -I $(DIR_MOD) $< -o $(DIR_OBJ_SPECVBIN)/spec_transformations.o

$(DIR_OBJ_SPECVBIN)/spec_setupray3d.o: $(DIR_SRC_SPECVBIN)/spec_setupray3d.f90   \
                                       $(DIR_OBJ_SPECVBIN)/spec_transformations.o   \
                                       $(DIR_OBJ)/info_region.o                  \
                                       $(DIR_OBJ)/math.o                         \
                                       $(OBJSM_SPECVBIN) 
		$(F90) $(CFLAGS) -I $(DIR_MOD_SPECVBIN) -I $(DIR_MOD) $< -o $(DIR_OBJ_SPECVBIN)/spec_setupray3d.o

$(DIR_OBJ_SPECVBIN)/spec_input.o: $(DIR_SRC_SPECVBIN)/spec_input.f90 \
                                  $(DIR_OBJ)/model_laws.o $(DIR_OBJ)/mod_interp1d.o \
                                  $(DIR_OBJ)/math.o \
                                  $(DIR_OBJ)/info_region.o \
                                  $(OBJSM_SPECVBIN)
		$(F90) $(CFLAGS) -I $(DIR_MOD) -I $(DIR_MOD_SPECVBIN) -I $(DIR_MOD_HDF5) $< -o $(DIR_OBJ_SPECVBIN)/spec_input.o

$(DIR_OBJ_SPECVBIN)/spec_output.o: $(DIR_SRC_SPECVBIN)/spec_output.f90 \
                                   $(OBJSM_SPECVBIN)
		$(F90) $(CFLAGS) -I $(DIR_MOD) -I $(DIR_MOD_SPECVBIN) -I $(DIR_MOD_HDF5) $< -o $(DIR_OBJ_SPECVBIN)/spec_output.o

$(DIR_OBJ_SPECVBIN)/spec.o: $(DIR_SRC_SPECVBIN)/spec.f90 \
                            $(OBJSM_SPECVBIN) \
                            $(OBJS1_SPECVBIN)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SPECVBIN) -I $(DIR_MOD) $< -o $(DIR_OBJ_SPECVBIN)/spec.o

$(DIR_OBJ_SPECVBIN)/mod_spec.o: $(DIR_SRC_SPECVBIN)/mod_spec.f90 \
                                $(OBJSM_TYPE)
		$(F90) $(CFLAGS) $(MFLAGS) $(DIR_MOD_SPECVBIN) -I $(DIR_MOD) $< -o $(DIR_OBJ_SPECVBIN)/mod_spec.o

########################################################################

$(DIR_OBJ_MODELSPECVBIN)/modelspec.o: $(DIR_SRC_MODELSPECVBIN)/modelspec.f90 \
                                      $(OBJSM_MODELSPECVBIN) \
                                      $(OBJS1_MODELSPECVBIN)
		$(F90) $(CFLAGS) -I $(DIR_MOD) -I $(DIR_MOD_MODELSPECVBIN) -I $(DIR_MOD_HDF5) -I $(DIR_MOD_LTE) $< -o $(DIR_OBJ_MODELSPECVBIN)/modelspec.o

$(DIR_OBJ_MODELSPECVBIN)/mod_modelspec.o: $(DIR_SRC_MODELSPECVBIN)/mod_modelspec.f90 \
                                          $(OBJSM_TYPE)
		$(F90) $(CFLAGS) $(MFLAGS) $(DIR_MOD_MODELSPECVBIN) -I $(DIR_MOD) $< -o $(DIR_OBJ_MODELSPECVBIN)/mod_modelspec.o

########################################################################

# inputparam.mod gridstratrev.mod comvelo.mod comtluc.mod
.PHONY: clean cleanall

clean:	
		rm -f *.o *~ *.mod \#*
		rm -f $(DIR_OBJ)/*.o $(DIR_OBJ)/*~ $(DIR_OBJ)/*.mod
		rm -f $(DIR_OBJ_MODEL)/*.o $(DIR_OBJ_MODEL)/*~ $(DIR_OBJ_MODEL)/*.mod
		rm -f $(DIR_OBJ_SC3D)/*.o $(DIR_OBJ_SC3D)/*~ $(DIR_OBJ_SC3D)/*.mod $(DIR_MOD_SC3D)/*.mod
		rm -f $(DIR_OBJ_PHOTPROF)/*.o $(DIR_OBJ_PHOTPROF)/*~ $(DIR_OBJ_PHOTPROF)/*.mod $(DIR_MOD_PHOTPROF)/*.mod
		rm -f $(DIR_OBJ_SPEC)/*.o $(DIR_OBJ_SPEC)/*~ $(DIR_OBJ_SPEC)/*.mod $(DIR_MOD_SPEC)/*.mod
		rm -f $(DIR_OBJ_MODELSPEC)/*.o $(DIR_OBJ_MODELSPEC)/*~ $(DIR_OBJ_MODELSPEC)/*.mod $(DIR_MOD_MODELSPEC)/*.mod
		rm -f $(DIR_OBJ_SPECVBIN)/*.o $(DIR_OBJ_SPECVBIN)/*~ $(DIR_OBJ_SPECVBIN)/*.mod $(DIR_MOD_SPECVBIN)/*.mod
		rm -f $(DIR_OBJ_MODELSPECVBIN)/*.o $(DIR_OBJ_MODELSPECVBIN)/*~ $(DIR_OBJ_MODELSPECVBIN)/*.mod $(DIR_MOD_MODELSPECVBIN)/*.mod
		rm -f $(DIR_OBJ_OPAL)/*.o $(DIR_OBJ_OPAL)/*~ $(DIR_OBJ_OPAL)/*.mod $(DIR_MOD_OPAL)/*.mod
		rm -f $(DIR_OBJ_LTE)/*.o $(DIR_OBJ_LTE)/*~ $(DIR_OBJ_LTE)/*.mod $(DIR_MOD_LTE)/*.mod

		rm -f $(DIR_MOD)/*.o $(DIR_MOD)/*~ $(DIR_MOD)/*.mod
		rm -f $(DIR_MOD_MODEL)/*.o $(DIR_MOD_MODEL)/*~ $(DIR_MOD_MODEL)/*.mod
		rm -f $(DIR_MOD_SC3D)/*.o $(DIR_MOD_SC3D)/*~ $(DIR_MOD_SC3D)/*.mod
		rm -f $(DIR_MOD_PHOTPROF)/*.o $(DIR_MOD_PHOTPROF)/*~ $(DIR_MOD_PHOTPROF)/*.mod
		rm -f $(DIR_MOD_SPEC)/*.o $(DIR_MOD_SPEC)/*~ $(DIR_MOD_SPEC)/*.mod
		rm -f $(DIR_MOD_MODELSPEC)/*.o $(DIR_MOD_MODELSPEC)/*~ $(DIR_MOD_MODELSPEC)/*.mod
		rm -f $(DIR_MOD_SPECVBIN)/*.o $(DIR_MOD_SPECVBIN)/*~ $(DIR_MOD_SPECVBIN)/*.mod
		rm -f $(DIR_MOD_MODELSPECVBIN)/*.o $(DIR_MOD_MODELSPECVBIN)/*~ $(DIR_MOD_MODELSPECVBIN)/*.mod
		rm -f $(DIR_MOD_OPAL)/*.o $(DIR_MOD_OPAL)/*~ $(DIR_MOD_OPAL)/*.mod
		rm -f $(DIR_MOD_LTE)/*.o $(DIR_MOD_LTE)/*~ $(DIR_MOD_LTE)/*.mod


		rm -f $(DIR_SRC)/*~ $(DIR_SRC)/#*
		rm -f $(DIR_SRC_MODEL)/*~ $(DIR_SRC_MODEL)/#*
		rm -f $(DIR_SRC_SC3D)/*~ $(DIR_SRC_SC3D)/#*
		rm -f $(DIR_SRC_PHOTPROF)/*~ $(DIR_SRC_PHOTPROF)/#*
		rm -f $(DIR_SRC_SPEC)/*~ $(DIR_SRC_SPEC)/#*
		rm -f $(DIR_SRC_MODELSPEC)/*~ $(DIR_SRC_MODELSPEC)/#*
		rm -f $(DIR_SRC_SPECVBIN)/*~ $(DIR_SRC_SPECVBIN)/#*
		rm -f $(DIR_SRC_MODELSPECVBIN)/*~ $(DIR_SRC_MODELSPECVBIN)/#*
		rm -f $(DIR_SRC_OPAL)/*~ $(DIR_SRC_OPAL)/#*
		rm -f $(DIR_SRC_LTE)/*~ $(DIR_SRC_LTE)/#*


cleanall: 
		rm -f *.eo *.o *~ *.mod \#*

