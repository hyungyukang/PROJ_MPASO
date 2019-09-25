.SUFFIXES: .F .c .o

ifneq "$(ESM)" ""

include Makefile.in.$(ESM)

else

ifeq "$(AUTOCLEAN)" "true"
AUTOCLEAN_DEPS=clean_shared
else
AUTOCLEAN_DEPS=
endif

#setenv CRAYPE_LINK_TYPE dynamic
#setenv XTPE_LINK_TYPE dynamic

FDIR=/global/u2/h/hgkang/my_programs/Trilinos_MPI/build/Trilinos_install/lib
FDIR_inc=/global/u2/h/hgkang/my_programs/Trilinos_MPI/build/Trilinos_install/include

#libsrcs:=$(FDIR)/$(wildcard *.a)

all: mpas

mpas: $(AUTOCLEAN_DEPS) externals frame ops dycore drver
#	$(LINKER)  -g  $(LDFLAGS) -o $(EXE_NAME) driver/*.o -L. -ldycore -lops -lframework  $(LIBS) -I./external/esmf_time_f90 -L./external/esmf_time_f90 -lesmf_time $(FDIR)/libfortrilinos.a $(FDIR)/libforteuchos.a $(FDIR)/libfortpetra.a $(FDIR)/libforerror.a $(FDIR)/libforutils.a $(FDIR)/libforbelos.a $(FDIR)/libifpack2.a $(FDIR)/libbelos.a $(FDIR)/libnoxlapack.a $(FDIR)/libnox.a $(FDIR)/libthyratpetra.a $(FDIR)/libteuchoscomm.a $(FDIR)/libteuchosparameterlist.a $(FDIR)/libteuchoscore.a $(FDIR)/libtpetra.a $(FDIR)/libloca.a $(FDIR)/libtpetraclassic.a $(FDIR)/libnoxlapack.a $(FDIR)/libkokkoskernels.a $(FDIR)/libkokkostsqr.a $(FDIR)/libkokkoscore.a $(FDIR)/libkokkosalgorithms.a $(FDIR)/libkokkoscontainers.a $(FDIR)/libteuchoskokkoscompat.a  $(FDIR)/libthyracore.a $(FDIR)/liblocathyra.a $(FDIR)/libstratimikos.a  $(FDIR)/libanasazi.a $(FDIR)/libbelostpetra.a $(FDIR)/libloca.a $(FDIR)/libxpetra.a $(FDIR)/libteuchoscomm.a   /opt/cray/pe/libsci/19.02.1/GNU/8.1/x86_64/lib/libsci_gnu.a /usr/lib64/libdl.a -lmvec -lnetcdf_c++4 -lmpichcxx_gnu_82 -lmvec -lstdc++ 
#	$(LINKER)  -g  $(LDFLAGS) -o $(EXE_NAME) driver/*.o -L. -ldycore -lops -lframework  $(LIBS) -I./external/esmf_time_f90 -L./external/esmf_time_f90 -lesmf_time $(libsrcs) -lmvec -lnetcdf_c++4 -lmpichcxx_gnu_82 -lmvec -lstdc++ 
#	$(LINKER)  -g  $(LDFLAGS) -o $(EXE_NAME) driver/*.o -L. -ldycore -lops -lframework  $(LIBS) -I./external/esmf_time_f90 -L./external/esmf_time_f90 -lesmf_time -static -L$(FDIR)/libamesos2.a -L$(FDIR)/libanasazi.a -L$(FDIR)/libanasazitpetra.a -L$(FDIR)/libbelos.a -L$(FDIR)/libbelostpetra.a -L$(FDIR)/libbelosxpetra.a -L$(FDIR)/libforbelos.a -L$(FDIR)/libforerror.a -L$(FDIR)/libforteuchos.a -L$(FDIR)/libfortpetra.a -L$(FDIR)/libfortrilinos.a -L$(FDIR)/libforutils.a -L$(FDIR)/libgaleri-xpetra.a -L$(FDIR)/libifpack2-adapters.a -L$(FDIR)/libifpack2.a -L$(FDIR)/libkokkosalgorithms.a -L$(FDIR)/libkokkoscontainers.a -L$(FDIR)/libkokkoscore.a -L$(FDIR)/libkokkoskernels.a -L$(FDIR)/libkokkostsqr.a -L$(FDIR)/libloca.a -L$(FDIR)/liblocalapack.a -L$(FDIR)/liblocathyra.a -L$(FDIR)/libmuelu-adapters.a -L$(FDIR)/libmuelu-interface.a -L$(FDIR)/libmuelu.a -L$(FDIR)/libmuelu_lgn.a -L$(FDIR)/libnox.a -L$(FDIR)/libnoxlapack.a -L$(FDIR)/librtop.a -L$(FDIR)/libstratimikos.a -L$(FDIR)/libstratimikosamesos2.a -L$(FDIR)/libstratimikosbelos.a -L$(FDIR)/libteuchoscomm.a -L$(FDIR)/libteuchoscore.a -L$(FDIR)/libteuchoskokkoscomm.a -L$(FDIR)/libteuchoskokkoscompat.a -L$(FDIR)/libteuchosnumerics.a  -L$(FDIR)/libteuchosparameterlist.a -L$(FDIR)/libteuchosparser.a -L$(FDIR)/libteuchosremainder.a -L$(FDIR)/libthyracore.a -L$(FDIR)/libthyratpetra.a -L$(FDIR)/libtpetra.a -L$(FDIR)/libtpetraclassic.a -L$(FDIR)/libtpetraclassiclinalg.a -L$(FDIR)/libtpetraclassicnodeapi.a -L$(FDIR)/libtpetraext.a -L$(FDIR)/libtpetrainout.a -L$(FDIR)/libtrilinosss.a -L$(FDIR)/libxpetra-sup.a -L$(FDIR)/libxpetra.a  -lmvec -lmpichcxx_gnu_82 -lmvec -lstdc++ 

	$(LINKER)  -g  $(LDFLAGS) -o $(EXE_NAME) driver/*.o -L. -ldycore -lops -lframework  $(LIBS) -I./external/esmf_time_f90 -L./external/esmf_time_f90 -lesmf_time -Wl,-rpath,/global/homes/h/hgkang/my_programs/Trilinos_MPI/build/ForTrilinos/src/interface/src:/global/homes/h/hgkang/my_programs/Trilinos_MPI/build/packages/ifpack2/adapters:/global/homes/h/hgkang/my_programs/Trilinos_MPI/build/packages/ifpack2/src:/global/homes/h/hgkang/my_programs/Trilinos_MPI/build/ForTrilinos/src/tpetra/src:/global/homes/h/hgkang/my_programs/Trilinos_MPI/build/ForTrilinos/src/teuchos/src:/global/homes/h/hgkang/my_programs/Trilinos_MPI/build/ForTrilinos/src/utils/src:/global/homes/h/hgkang/my_programs/Trilinos_MPI/build/packages/nox/src:/global/homes/h/hgkang/my_programs/Trilinos_MPI/build/packages/thyra/adapters/tpetra/src:$(FDIR) -dynamic -L$(FDIR) -lfortrilinos -lforteuchos -lfortpetra -lnox -lforerror -lforutils -lstratimikos -lstratimikosbelos

#	$(LINKER)  -g  $(LDFLAGS) -o $(EXE_NAME) driver/*.o -L. -ldycore -lops -lframework  $(LIBS) -I./external/esmf_time_f90 -L./external/esmf_time_f90 -lesmf_time -L$(FDIR) -lfortrilinos
externals: $(AUTOCLEAN_DEPS)
	( cd external; $(MAKE) FC="$(FC)" SFC="$(SFC)" CC="$(CC)" SCC="$(SCC)" FFLAGS="$(FFLAGS)" CFLAGS="$(CFLAGS)" CPP="$(CPP)" NETCDF="$(NETCDF)" CORE="$(CORE)" all )

drver:  $(AUTOCLEAN_DEPS) externals frame ops dycore
	( cd driver; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" all ) 
endif

build_tools: externals
	(cd tools; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CC="$(SCC)" CFLAGS="$(CFLAGS)")

frame: $(AUTOCLEAN_DEPS) externals
	( cd framework; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" all ) 
	ln -sf framework/libframework.a libframework.a

ops: $(AUTOCLEAN_DEPS) externals frame
	( cd operators; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" all ) 
	ln -sf operators/libops.a libops.a

dycore: $(AUTOCLEAN_DEPS) build_tools externals frame ops
	( cd core_$(CORE); $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" REG_PARSE="$(PWD)/tools/registry/parse" gen_includes )
	( cd core_$(CORE); $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" NL_GEN="$(PWD)/tools/input_gen/namelist_gen" ST_GEN="$(PWD)/tools/input_gen/streams_gen" core_input_gen )
	( cd core_$(CORE); $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" all ) 
	ln -sf core_$(CORE)/libdycore.a libdycore.a


clean: clean_shared clean_core

clean_core:
	if [ -d core_$(CORE) ] ; then \
	   ( cd core_$(CORE); $(MAKE) clean ) \
	fi;

clean_shared:
ifeq "$(AUTOCLEAN)" "true"
	@echo ""
	@echo "*********************************************************************************************"
	@echo "The MPAS infrastructure is currently built for a core different from $(CORE)."
	@echo "The infrastructure will be cleaned and re-built for the $(CORE) core."
	@echo "*********************************************************************************************"
	@echo ""
endif
	$(RM) libframework.a libops.a libdycore.a lib$(CORE).a *.o
	( cd tools; $(MAKE) clean )
	( cd external; $(MAKE) clean )
	( cd framework; $(MAKE) clean )
	( cd operators; $(MAKE) clean )
	( cd driver; $(MAKE) clean )
