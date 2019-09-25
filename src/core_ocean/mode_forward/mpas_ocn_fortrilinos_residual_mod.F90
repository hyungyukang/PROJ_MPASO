module ocn_residual_mod

   use,intrinsic :: iso_c_binding

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer
   use mpas_threading
   use mpas_timekeeping
   use mpas_log

!  use ocn_tendency
!  use ocn_diagnostics
!  use ocn_gm

!  use ocn_equation_of_state
!  use ocn_vmix
!  use ocn_time_average_coupled

!  use ocn_effective_density_in_land_ice
!
!  use ocn_model_evaluator_mod, only : ocnModelEvaluator

   implicit none
   private
   save

!  public :: residual

   character (len=*), parameter :: subcycleGroupName = 'subcycleFields'

!  integer :: nBtrSubcycles
!  real (kind=RKIND), allocatable,dimension(:,:) :: block_prec_ivmat,block_prec_amat
!  real (kind=RKIND) :: total_num_cells,area_mean
!  integer :: istep=0

!***********************************************************************
                               contains
!***********************************************************************

  !=======================================================================

! subroutine residual(xstate,fx,nvec,dt)
! subroutine residual(xstate,fx,nvec,dt,nCells,nEdges,                             &
!                     nCellsArray,s1_nCellsArray,                                  &
!                     nEdgesArray,s1_nEdgesArray,                                  &
!                     nEdgesOnCell,s1_nEdgesOnCell,                                &
!                     dvEdge,s1_dvEdge,                                            & 
!                     dcEdge,s1_dcEdge,                                            & 
!                     areaCell,s1_areaCell,                                        &
!                     bottomDepth,s1_bottomDepth,                                  & 
!                     cellsOnEdge,s1_cellsOnEdge,s2_cellsOnEdge,                    &
!                     edgesOnCell,s1_edgesOnCell,s2_edgesOnCell,                   &
!                     edgeSignOnCell,s1_edgeSignOnCell,s2_edgeSignOnCell,          &
!                     sshCur, s1_sshCur,                                           &
!                     sshSubcycleCur,s1_sshSubcycleCur,                            &
!                     normalBarotropicVelocityCur, s1_normalBarotropicVelocityCur, &
!                     barotropicForcing, s1_barotropicForcing,                     &
!                     barotropicCoriolisterm, s1_barotropicCoriolisTerm,           &
!                     CGvec_r0,s1_CGvec_r0)
  subroutine residual(self,xstate,fx,nvec)

    !---------------------------------------------------------------------
    implicit none
    class(ocnModelEvaluator),intent(in) :: self
    real (c_double),dimension(nvec),intent(in)  :: xstate
    real (c_double),dimension(nvec),intent(inout) :: fx
    integer(c_int), intent(in), value :: nvec
    real (kind=RKIND) :: dt
    !---------------------------------------------------------------------
    type (domain_type) :: domain
    type (block_type),pointer :: block
    type (mpas_pool_type), pointer :: statePool
    type (mpas_pool_type), pointer :: tracersPool
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: diagnosticsPool
    type (mpas_pool_type), pointer :: tendPool
    type (mpas_pool_type), pointer :: tracersTendPool

    integer :: nCells, nEdges
    integer :: s1_nCellsArray,s1_nEdgesArray
    integer :: s1_nEdgesOnCell,s1_dvEdge,s1_dcEdge,s1_areaCell,s1_bottomDepth
    integer :: s1_cellsOnEdge,s2_cellsOnEdge,s1_edgesOnCell,s2_edgesOnCell
    integer :: s1_edgeSignOnCell,s2_edgeSignOnCell
    integer :: s1_sshCur,s1_sshSubcycleCur,s1_normalBarotropicVelocityCur
    integer :: s1_barotropicForcing, s1_barotropicCoriolisTerm,s1_CGvec_r0
    integer :: i,j,iCell,iEdge,cell1,cell2

!   integer :: nCellsPtr, nEdgesPtr    
!   integer, dimension(s1_nCellsArray) :: nCellsArray
!   integer, dimension(s1_nEdgesArray) :: nEdgesArray
!   integer, dimension(s1_nEdgesOnCell)   :: nEdgesOnCell    
!   integer, dimension(s1_cellsOnEdge,s2_cellsOnEdge)       :: cellsOnEdge
!   integer, dimension(s1_edgesOnCell,s2_edgesOnCell)       :: edgesOnCell
!   integer, dimension(s1_edgeSignOnCell,s2_edgeSignOnCell) :: edgeSignOnCell
!   real (kind=RKIND) :: sshTendb1,sshTendb2,sshTendAx,sshEdgeCur,thicknessSumCur,sshDiffCur
!   real (kind=RKIND) :: fluxb1,fluxb2,fluxAx,sshCurArea
!   real (kind=RKIND), dimension(s1_dcEdge)      :: dcEdge
!   real (kind=RKIND), dimension(s1_bottomDepth) :: bottomDepth
!   real (kind=RKIND), dimension(s1_dvEdge)      :: dvEdge
!   real (kind=RKIND), dimension(s1_areaCell)    :: areaCell

!   real (kind=RKIND), dimension(s1_sshCur) :: sshCur
!   real (kind=RKIND), dimension(s1_sshSubcycleCur) :: sshSubcycleCur
!   real (kind=RKIND), dimension(s1_normalBarotropicVelocityCur) :: normalBarotropicVelocityCur
!   real (kind=RKIND), dimension(s1_barotropicForcing) :: barotropicForcing
!   real (kind=RKIND), dimension(s1_barotropicCoriolisTerm) :: barotropicCoriolisTerm
!   real (kind=RKIND), dimension(s1_CGvec_r0) :: CGvec_r0

    integer,pointer :: nCellsPtr, nEdgesPtr    
    integer, dimension(:),pointer :: nCellsArray
    integer, dimension(:),pointer :: nEdgesArray
    integer, dimension(:),pointer :: nEdgesOnCell    
    integer, dimension(:,:),pointer :: cellsOnEdge
    integer, dimension(:,:),pointer :: edgesOnCell
    integer, dimension(:,:),pointer :: edgeSignOnCell
    real (kind=RKIND) :: sshTendb1,sshTendb2,sshTendAx,sshEdgeCur,thicknessSumCur,sshDiffCur
    real (kind=RKIND) :: fluxb1,fluxb2,fluxAx,sshCurArea
    real (kind=RKIND), dimension(:),pointer :: dcEdge
    real (kind=RKIND), dimension(:),pointer :: bottomDepth
    real (kind=RKIND), dimension(:),pointer :: dvEdge
    real (kind=RKIND), dimension(:),pointer :: areaCell

    real (kind=RKIND), dimension(:),pointer :: sshCur
    real (kind=RKIND), dimension(:),pointer :: sshSubcycleCur
    real (kind=RKIND), dimension(:),pointer :: normalBarotropicVelocityCur
    real (kind=RKIND), dimension(:),pointer :: barotropicForcing
    real (kind=RKIND), dimension(:),pointer :: barotropicCoriolisTerm
    real (kind=RKIND), dimension(:),pointer :: CGvec_r0
    !---------------------------------------------------------------------

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! Stage 2.2 : Compute : fx = A*x - b 
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!        call mpas_timer_start("si btr first r0")

         block => self%domain % blocklist
         do while (associated(block))
           call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
           call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
           call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
           call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

           call mpas_pool_get_subpool(block % structs, 'tend'       , tendPool       )
           call mpas_pool_get_subpool(tendPool       , 'tracersTend', tracersTendPool)
           call mpas_pool_get_subpool(block % structs, 'mesh'       , meshPool       )
           call mpas_pool_get_subpool(block % structs, 'state'      , statePool      )
           call mpas_pool_get_subpool(statePool      , 'tracers'    , tracersPool    )
           call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

           call mpas_pool_get_array(meshPool, 'nEdgesOnCell'  , nEdgesOnCell  )
           call mpas_pool_get_array(meshPool, 'edgesOnCell'   , edgesOnCell   )
           call mpas_pool_get_array(meshPool, 'cellsOnEdge'   , cellsOnEdge   )
           call mpas_pool_get_array(meshPool, 'dcEdge'        , dcEdge        )
           call mpas_pool_get_array(meshPool, 'bottomDepth'   , bottomDepth   )
           call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
           call mpas_pool_get_array(meshPool, 'dvEdge'        , dvEdge        )
           call mpas_pool_get_array(meshPool, 'areaCell'      , areaCell      )

           call mpas_pool_get_array(statePool, 'ssh', sshCur,1)
           call mpas_pool_get_array(statePool, 'ssh', sshSubcycleCur,1)
           call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur,1)

           call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
           call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm',barotropicCoriolisTerm)
           call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0', CGvec_r0)

           nCells = nCellsArray( 1 )
           nEdges = nEdgesArray( 2 )

           sshSubcycleCur(1:nCells) = xstate(1:nCells)

           call mpas_threading_barrier()
           call mpas_timer_start("si halo iter r0")
           call mpas_dmpar_exch_group_create(domain, subcycleGroupName)
           call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'sshSubcycleCur', 1)
           call mpas_threading_barrier()
           call mpas_dmpar_exch_group_full_halo_exch(domain, subcycleGroupName)
           call mpas_dmpar_exch_group_destroy(domain, subcycleGroupName)
           call mpas_timer_stop("si halo iter r0")

           print*, nCells,nEdges
           print*, 'STOP in residual routine  ******'
           STOP
       
           
           ! SpMV : A*x --------------------------------------------------------------------

           do iCell = 1, nCells
             sshTendb1 = 0.0_RKIND
             sshTendb2 = 0.0_RKIND
             sshTendAx = 0.0_RKIND
             do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)

               cell1 = cellsOnEdge(1, iEdge)
               cell2 = cellsOnEdge(2, iEdge)

               ! Interpolation sshEdge
               sshEdgeCur = 0.5_RKIND * (sshSubcycleCur(cell1) + sshSubcycleCur(cell2))

               ! method 1, matches method 0 without pbcs, works with pbcs.
               thicknessSumCur = sshEdgeCur + min(bottomDepth(cell1), bottomDepth(cell2))

               ! nabla (ssh^0)
               sshDiffCur = (  sshSubcycleCur(cell2) -   sshSubcycleCur(cell1)) / dcEdge(iEdge)

               !--------------------------------------------------------------!
               fluxb1 = thicknessSumCur * normalBarotropicVelocityCur(iEdge)
               fluxb2 = thicknessSumCur * (0.5_RKIND*gravity*sshDiffCur + (-barotropicCoriolisTerm(iEdge)-barotropicForcing(iEdge)))
               fluxAx = thicknessSumCur * sshDiffCur

               sshTendb1 = sshTendb1 + edgeSignOnCell(i, iCell) * fluxb1 * dvEdge(iEdge)
               sshTendb2 = sshTendb2 + edgeSignOnCell(i, iCell) * fluxb2 * dvEdge(iEdge)
               sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) * fluxAx * dvEdge(iEdge)
               !--------------------------------------------------------------!

             end do ! i

               sshTendb1 = (4.0_RKIND/(gravity*dt)) * sshTendb1
               sshTendb2 = (2.0_RKIND/(gravity   )) * sshTendb2

               sshCurArea = (4.0_RKIND/(gravity*dt**2.0)) *   sshCur(iCell) * areaCell(iCell)

               CGvec_r0(iCell) = -(-sshCurArea - sshTendb1 + sshTendb2)   &
                                 +(-sshCurArea - sshTendAx)

           end do ! iCell

           ! -------------------------------------------------------------------------------

           block => block % next
         end do  ! block

!        call mpas_threading_barrier()
!        call mpas_timer_start("si halo iter r0")
!        call mpas_dmpar_exch_group_create(domain, subcycleGroupName)
!        call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'CGvec_r0', 1)
!        call mpas_threading_barrier()
!        call mpas_dmpar_exch_group_full_halo_exch(domain, subcycleGroupName)
!        call mpas_dmpar_exch_group_destroy(domain, subcycleGroupName)
!        call mpas_timer_stop("si halo iter r0")

         fx(1:nCells) = CGvec_r0(1:nCells) 

!        call mpas_timer_stop ("si btr first r0")

  end subroutine residual

end module ocn_residual_mod

! vim: foldmethod=marker
