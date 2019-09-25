module fortrilinos_mod
#include "ForTrilinos.h"
  use forteuchos
  use fortpetra
  use fortrilinos
  use model_evaluator_mod
  use derived_type_mod, only : derived_type
  use iso_fortran_env
  implicit none

private

  ! Module data
  type(E3SMModelEvaluator), private, allocatable :: model_evaluator

  ! public interface
  public :: noxinit
  public :: noxsolve
  public :: noxfinish

contains

  subroutine noxinit(np, nlev, nelemd, state_vector, state_data, fcomm)
    use ,intrinsic :: iso_c_binding
    integer(c_int)               :: np, nlev, nelemd
    integer                      :: fcomm
    integer(size_type)           :: nState
    real(c_double), dimension(*) :: state_vector
    type(derived_type), pointer  :: state_data
    type(TpetraMap)              :: map
    type(TeuchosComm)            :: comm
    type(ParameterList)          :: params

#ifdef _PRIM
    ! PRIM
    nState = (3*np*np*nlev + np*np) * nelemd
#else
    ! SWIM
    nState = (3*np*np*nlev*nelemd)
#endif
    comm = TeuchosComm(fcomm); FORTRILINOS_CHECK_IERR()
    map = TpetraMap(TPETRA_GLOBAL_INVALID, nState, comm); FORTRILINOS_CHECK_IERR()

    ! TODO Read parameters from the input deck
    ! See homme/utils/trilinos/trilinosNoxSolver.cpp:253

    params = ParameterList()

    ! Allocate tpetra_me and initialize
    allocate(model_evaluator, source=E3SMModelEvaluator(map))
    call init_ForModelEvaluator(model_evaluator); FORTRILINOS_CHECK_IERR()
    call model_evaluator%setup(params); FORTRILINOS_CHECK_IERR()

    call params%release()
  end subroutine noxinit

  subroutine noxsolve(vector_size, state_vector, state_data, ierr)
    use, intrinsic :: iso_c_binding
    integer(c_int)               :: vector_size  ! length of state vector
    real(c_double), dimension(*) :: state_vector ! 1d state vector
    type(derived_type), pointer  :: state_data   ! derived type
    integer(c_int)               :: ierr         ! error flag

    type(NOXSolver)              :: nox_solver
    integer(kind(NOXStatusType)) :: status
    type(ParameterList)          :: params

    ierr = 0

    params = ParameterList()
    call load_from_xml(params, 'fortrilinosOptions.xml'); FORTRILINOS_CHECK_IERR()

    ! reset solver with new initial guess
    ! TODO?
    ! call model_evaluator%update_solution_vector(state_vec);

    ! reset interface with new state data
    call model_evaluator%update_state_data(state_data)
    call model_evaluator%setup(params); FORTRILINOS_CHECK_IERR()

    nox_solver = NOXSolver(model_evaluator)
    call nox_solver%setup(params); FORTRILINOS_CHECK_IERR()
    status = nox_solver%solve(); FORTRILINOS_CHECK_IERR()

    ! FIXME: where is the solution?
    call nox_solver%release(); FORTRILINOS_CHECK_IERR()

    if (status /= NOXConverged) ierr = 1

    call params%release()
  end subroutine noxsolve

  subroutine noxfinish()
    call model_evaluator%release(); FORTRILINOS_CHECK_IERR()
    deallocate(model_evaluator)
  end subroutine noxfinish

end module
