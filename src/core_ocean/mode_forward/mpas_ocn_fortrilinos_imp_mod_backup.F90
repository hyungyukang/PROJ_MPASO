module ocn_fortrilinos_imp_mod
!---------------------------
#include "ForTrilinosInterface_config.hpp"
#include "ForTrilinos.h"
  use forteuchos
  use fortpetra
  use fortrilinos
  use ocn_model_evaluator_mod
!---------------------------
  use, intrinsic :: iso_c_binding
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_constants
  use mpas_dmpar

! use mpas_vector_reconstruction
! use mpas_spline_interpolation
! use mpas_timer
! use mpas_threading
! use mpas_timekeeping
! use mpas_log

! use ocn_tendency
! use ocn_diagnostics
! use ocn_gm

! use ocn_equation_of_state
! use ocn_vmix
! use ocn_time_average_coupled

! use ocn_effective_density_in_land_ice
! use ocn_fortrilinos_mod
!---------------------------

  implicit none
  private  
  save


  public :: ocn_time_integration_imp_btrmode
! public :: noxinit ,noxsolve !,noxfinish
  public :: noxsolve !,noxfinish

!*********************************************************************
                                contains
!*********************************************************************

! subroutine ocn_time_integration_imp_btrmode(domain,nvec,xstate,fx)
  subroutine ocn_time_integration_imp_btrmode(domain,dt,xstate,nvec)
   
    implicit none
    !--------------
    type (TeuchosComm) :: comm
    type (ParameterList) :: params
    !--------------
    type (domain_type) :: domain
    real (kind=RKIND)  :: dt
!   real (c_double), dimension(nvec) :: xstate,fx
    real (c_double), dimension(nvec) :: xstate
    integer (c_int) :: ierr,nvec
    !--------------

    !--------------------------------------------------!
!   interface
!     subroutine noxsolve (comm,domain,dt,params,xstate,nvec,ierr) bind(C,name='noxsolve')
!       use,intrinsic :: iso_c_binding
!       use mpas_derived_types
!       type (TeuchosComm) :: comm
!       type (domain_type) :: domain
!       type (ParameterList) :: params
!       real (kind=RKIND)  :: dt
!       real (c_double), dimension(nvec) :: xstate
!       integer(c_int) :: ierr,nvec
!     end subroutine noxsolve
!   end interface
    !--------------------------------------------------!
    
!   print*, 'STOP in btrmode'
!   STOP
!   print*, minval(xstate),maxval(xstate)

    comm = TeuchosComm()
    params = ParameterList("ocnModelEvaluator")
    call load_from_xml(params,'nox_params.xml')

    call noxsolve(comm,domain,dt,params,xstate,nvec,ierr)

    call params%release()
    call comm%release()

  end subroutine ocn_time_integration_imp_btrmode

!------------------------------------------------------------------------

! subroutine noxinit(domain,fcomm)
  subroutine noxsolve(comm,domain,dt,params,nvec)
    use,intrinsic :: iso_c_binding
    implicit none

    integer :: fcomm
    type (domain_type) :: domain
    real (kind=RKIND)  :: dt
    integer(c_int) :: vector_size,ierr,nvec
    type(ParameterList) :: params
!    integer(size_type) :: nCellsNox
!    type(TpetraMap)   :: map
!    type(TeuchosComm) :: comm
!    type(ParameterList) :: params
!
!    !--------------------------------
!      
!    type (domain_type) :: domain
!!   real (kind=RKIND)  :: dt
!    type (block_type),pointer :: block
!    integer,pointer :: nCellsPtr,nEdgesPtr
!    integer,pointer,dimension(:) :: nCellsArray
!
!    !--------------------------------
!    
!!   comm = TeuchosComm(fcomm)
!    comm = TeuchosComm()
!
!    block => domain % blocklist
!    do while (associated(block))
!
!      call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
!
!      nCellsNox = nCellsArray(1)
!
!      map = TpetraMap(nCellsNox,comm)
!
!      params = ParameterList("ocnModelEvaluator")
!      call load_from_xml(params, 'fortrilinosOptions.xml')
!
!      allocate(model_evaluator, source=ocnModelEvaluator(map))
!      call init_ForModelEvaluator(model_evaluator)
!      call model_evaluator%setup(params) 
!      call params%release()
!
!      block => block % next
!    end do
!
  end subroutine noxinit
 
!------------------------------------------------------------------------

  subroutine noxsolve(comm,domain,dt,params,xstate,nvec,ierr)
    ! -------------------------------------------------------------------
    use,intrinsic :: iso_c_binding
    implicit none

    type (TeuchosComm) :: comm
    type (domain_type) :: domain
    real (kind=RKIND)  :: dt
    integer(c_int) :: vector_size,ierr,nvec
    real(c_double),dimension(nvec) :: xstate

    class(ForModelEvaluator),allocatable :: model_evaluator
    type(NOXSolver) :: nox_solver
    type(ocnMOdelEvaluator) :: self
    integer(kind(NOXStatusType)) :: status
    type(ParameterList) :: params
    integer(global_size_type) :: n_global_vec
    ! -------------------------------------------------------------------

    ierr = 0

    n_global_vec = nvec
    
     
    allocate(model_evaluator, source=ocnModelEvaluator(comm,domain,xstate,dt,n_global_vec))
    call init_ForModelEvaluator(model_evaluator)


    call model_evaluator%setup(params)
    nox_solver = NOXSolver(model_evaluator)
    call nox_solver%setup(params)
    status = nox_solver%solve()

    if (status /= NOXConverged) ierr = 1

    !--- Get the final solution
    call getFinalSolution(self,xstate,n_global_vec)

    call model_evaluator%release()
    deallocate(model_evaluator)
  
    call params%release()
    call nox_solver%release()

  end subroutine noxsolve

end module
