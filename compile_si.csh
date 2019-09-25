#!/bin/csh

set MPAS_HOME = /global/u2/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR
set LAPACK_HOME = /global/u2/h/hgkang/my_programs/lapack-3.8.0/build/lapack
set LAPACK_LIB = ${LAPACK_HOME}/lib64
#set  LAPACK_HOME = $CRAY_LIBSCI_PREFIX_DIR
#set LAPACK_LIB = ${LAPACK_HOME}/lib
set PIO_HOME  = /global/u2/h/hgkang/my_programs/ParallelIO_1.9.23/build/pio
 
set exe = ocean_model

###################
cd ${MPAS_HOME}/src/core_ocean/mode_forward
rm -f mpas_ocn_forward_mode.o  mpas_ocn_forward_mode.mod

ftn -DUNDERSCORE -D_MPI   -DCORE_OCEAN -DMPAS_NAMELIST_SUFFIX=ocean -DMPAS_EXE_NAME=ocean_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=unknown -O3 -m64 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8 -c mpas_ocn_forward_mode.F -I/global/u2/h/hgkang/my_programs/ParallelIO_1.9.23/build/pio -I/global/u2/h/hgkang/my_programs/ParallelIO_1.9.23/build/pio -I${MPAS_HOME}/src/core_ocean/driver  -I${MPAS_HOME}/src/core_ocean/mode_forward -I${MPAS_HOME}/src/core_ocean/mode_analysis -I${MPAS_HOME}/src/core_ocean/mode_init -I${MPAS_HOME}/src/core_ocean/shared -I${MPAS_HOME}/src/core_ocean/analysis_members -I${MPAS_HOME}/src/core_ocean/cvmix -I${MPAS_HOME}/src/core_ocean/BGC -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/../framework -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/../external/esmf_time_f90 -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/../operators -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/BGC -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/shared -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/analysis_members -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/cvmix -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/mode_forward -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/mode_analysis -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/mode_init -L${LAPACK_LIB} -lblas -llapack

###################

###################
cd ${MPAS_HOME}/src/core_ocean/mode_forward
rm -f mpas_ocn_time_integration_si.o mpas_ocn_time_integration_si.mod

ftn -DUNDERSCORE -D_MPI   -DCORE_OCEAN -DMPAS_NAMELIST_SUFFIX=ocean -DMPAS_EXE_NAME=ocean_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=unknown -O3 -m64 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8 -c mpas_ocn_time_integration_si.F -I/global/u2/h/hgkang/my_programs/ParallelIO_1.9.23/build/pio -I/global/u2/h/hgkang/my_programs/ParallelIO_1.9.23/build/pio -I${MPAS_HOME}/src/core_ocean/driver  -I${MPAS_HOME}/src/core_ocean/mode_forward -I${MPAS_HOME}/src/core_ocean/mode_analysis -I${MPAS_HOME}/src/core_ocean/mode_init -I${MPAS_HOME}/src/core_ocean/shared -I${MPAS_HOME}/src/core_ocean/analysis_members -I${MPAS_HOME}/src/core_ocean/cvmix -I${MPAS_HOME}/src/core_ocean/BGC -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/../framework -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/../external/esmf_time_f90 -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/../operators -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/BGC -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/shared -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/analysis_members -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/cvmix -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/mode_forward -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/mode_analysis -I/global/homes/h/hgkang/MPAS/MPAS-Model-6.0_test_imp_pipeCR/src/core_ocean/mode_init -L${LAPACK_LIB} -lblas -llapack
###################

cd ${MPAS_HOME}/src/core_ocean

if ( -e libdycore.a ) then 
rm -f libdycore.a 
endif
echo "AR-RU .............."
ar -ru libdycore.a `find . -type f -name "*.o"`

###################
cd ${MPAS_HOME}/src
ln -sf core_ocean/libdycore.a libdycore.a
echo "mpas.o .............."

cd ${MPAS_HOME}/src/driver
rm -f mpas.o mpas.mod
ftn -DUNDERSCORE -D_MPI   -DCORE_OCEAN -DMPAS_NAMELIST_SUFFIX=ocean -DMPAS_EXE_NAME=ocean_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=unknown -O3 -m64 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8 -c mpas.F -I/global/u2/h/hgkang/my_programs/ParallelIO_1.9.23/build/pio -I/global/u2/h/hgkang/my_programs/ParallelIO_1.9.23/build/pio -I${MPAS_HOME}/src/core_ocean/driver  -I${MPAS_HOME}/src/core_ocean/mode_forward -I${MPAS_HOME}/src/core_ocean/mode_analysis -I${MPAS_HOME}/src/core_ocean/mode_init -I${MPAS_HOME}/src/core_ocean/shared -I${MPAS_HOME}/src/core_ocean/analysis_members -I${MPAS_HOME}/src/core_ocean/cvmix -I${MPAS_HOME}/src/core_ocean/BGC -I../framework -I../operators -I../core_ocean -I../external/esmf_time_f90

cd ${MPAS_HOME}/src

echo "EXE .............."
ftn -O3 -o $exe driver/*.o -L. -ldycore -lops -lframework -L${PIO_HOME} -lpio -I./external/esmf_time_f90 -L./external/esmf_time_f90 -lesmf_time
#ftn -O3 -o ocean_model_cg driver/*.o -L. -ldycore -lops -lframework -L/autofs/nccs-svm1_home1/hgkang/install_programs/ParallelIO_1.9.23_pgi2/build/pio -lpio -I./external/esmf_time_f90 -L./external/esmf_time_f90 -lesmf_time

cd ${MPAS_HOME}

if ( -e src/${exe} ) then
    mv src/${exe} .
endif
