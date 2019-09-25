! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_si
!
!> \brief MPAS ocean split explicit time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the routine for the split explicit
!>  time integration scheme
!
!-----------------------------------------------------------------------


module ocn_time_integration_si

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

   use ocn_tendency
   use ocn_diagnostics
   use ocn_gm

   use ocn_equation_of_state
   use ocn_vmix
   use ocn_time_average_coupled

   use ocn_effective_density_in_land_ice

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_time_integrator_si, ocn_time_integration_si_init

   character (len=*), parameter :: subcycleGroupName = 'subcycleFields'
   character (len=*), parameter :: finalBtrGroupName = 'finalBtrFields'
   integer :: nBtrSubcycles
   real (kind=RKIND) :: total_num_cells,area_mean
   integer,save :: istep=0

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_si
!
!> \brief MPAS ocean split explicit time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine integrates a master time step (dt) using a
!>  split explicit time integrator.
!
!-----------------------------------------------------------------------

    subroutine ocn_time_integrator_si_first(domain, dt)!{{{
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Advance model state forward in time by the specified time step using
    !   Split_Explicit timestepping scheme
    !
    ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:))
    !                 plus mesh meta-data
    ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains
    !                  model state advanced forward in time by dt seconds
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt

      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: verticalMeshPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: tendPool
      type (mpas_pool_type), pointer :: tracersTendPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: swForcingPool

      type (dm_info) :: dminfo
      integer :: iCell, i,k,j, iEdge, cell1, cell2, split_explicit_step, split, &
                 eoe, oldBtrSubcycleTime, newBtrSubcycleTime, uPerpTime, BtrCorIter, &
                 stage1_tend_time,jj, iter,iter_max
      integer, dimension(:), allocatable :: n_bcl_iter
      type (block_type), pointer :: block
      real (kind=RKIND) :: normalThicknessFluxSum, thicknessSum,thicknessSumCur,thicknessSumLag,thicknessSumMid,flux,flux1,flux2, sshEdge,sshEdgeCur,sshEdgeLag,sshEdgeMid, hEdge1, &
                 CoriolisTerm, normalVelocityCorrection, temp, temp_h, coef, barotropicThicknessFlux_coeff, sshCell1, sshCell2, &
                 sshDiffCur,sshDiffNew,sshDiffLag,sshDiffMid
      real (kind=RKIND) :: fluxb1,fluxb2,fluxAx,sshTendb1,sshTendb2,sshTendAx,sshCurArea,sshLagArea

      real (kind=RKIND) :: CoriolisTermCur,CoriolisTermNew
      integer :: useVelocityCorrection, err
      real (kind=RKIND), dimension(:,:), pointer :: &
                 vertViscTopOfEdge, vertDiffTopOfCell
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroup
      real (kind=RKIND), dimension(:), allocatable:: uTemp
      real (kind=RKIND), dimension(:), pointer :: btrvel_temp
      type (field1DReal), pointer :: btrvel_tempField
      real (kind=RKIND), dimension(:,:), allocatable:: tracersTemp
      logical :: activeTracersOnly ! if true only compute tendencies for active tracers  
      integer :: tsIter
      integer :: edgeHaloComputeCounter, cellHaloComputeCounter

      ! Config options
      character (len=StrKIND), pointer :: config_time_integrator
      integer, pointer :: config_n_bcl_iter_mid, config_n_bcl_iter_beg, config_n_bcl_iter_end
      integer, pointer :: config_n_ts_iter, config_btr_subcycle_loop_factor
      integer, pointer :: config_n_btr_cor_iter, config_num_halos
      logical, pointer :: config_use_standardGM

      logical, pointer :: config_use_freq_filtered_thickness, config_btr_solve_SSH2, config_filter_btr_mode
      logical, pointer :: config_vel_correction, config_prescribe_velocity, config_prescribe_thickness
      logical, pointer :: config_use_cvmix_kpp
      logical, pointer :: config_use_tracerGroup
      logical, pointer :: config_compute_active_tracer_budgets

      character (len=StrKIND), pointer :: config_land_ice_flux_mode

      real (kind=RKIND), pointer :: config_mom_del4, config_btr_gam1_velWt1, config_btr_gam2_SSHWt1
      real (kind=RKIND), pointer :: config_btr_gam3_velWt2

      ! Dimensions
      integer :: nCells, nEdges
      integer, pointer :: nCellsPtr, nEdgesPtr, nVertLevels, num_tracersGroup, startIndex, endIndex
      integer, pointer :: nCellsSolve, nEdgesSolve
      integer, pointer :: indexTemperature, indexSalinity
      integer, pointer :: indexSurfaceVelocityZonal, indexSurfaceVelocityMeridional
      integer, pointer :: indexSSHGradientZonal, indexSSHGradientMeridional
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray

      ! Mesh array pointers
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop, nEdgesOnEdge, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgeMask, edgesOnEdge
      integer, dimension(:,:), pointer :: edgesOnCell, edgeSignOnCell

      real (kind=RKIND), dimension(:), pointer :: dcEdge, fEdge, bottomDepth, refBottomDepthTopOfCell
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge

      ! State Array Pointers
      real (kind=RKIND), dimension(:), pointer :: sshSubcycleCur, sshSubcycleNew,sshSubcycleLag
      real (kind=RKIND), dimension(:), pointer :: normalBarotropicVelocitySubcycleCur, normalBarotropicVelocitySubcycleNew
      real (kind=RKIND), dimension(:), pointer :: sshCur, sshNew
      real (kind=RKIND), dimension(:), pointer :: normalBarotropicVelocityCur, normalBarotropicVelocityNew
      real (kind=RKIND), dimension(:,:), pointer :: normalBaroclinicVelocityCur, normalBaroclinicVelocityNew
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityCur, normalVelocityNew
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessCur, layerThicknessNew
      real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessCur, highFreqThicknessNew
      real (kind=RKIND), dimension(:,:), pointer :: lowFreqDivergenceCur, lowFreqDivergenceNew
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroupCur, tracersGroupNew

      ! Tend Array Pointers
      real (kind=RKIND), dimension(:), pointer :: sshTend
      real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessTend
      real (kind=RKIND), dimension(:,:), pointer :: lowFreqDivergenceTend
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityTend, layerThicknessTend
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroupTend, activeTracersTend

      ! Diagnostics Array Pointers
      real (kind=RKIND), dimension(:), pointer :: barotropicForcing, barotropicThicknessFlux
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge, normalTransportVelocity, normalGMBolusVelocity
      real (kind=RKIND), dimension(:,:), pointer :: vertAleTransportTop
      real (kind=RKIND), dimension(:,:), pointer :: velocityX, velocityY, velocityZ
      real (kind=RKIND), dimension(:,:), pointer :: velocityZonal, velocityMeridional
      real (kind=RKIND), dimension(:), pointer :: gradSSH
      real (kind=RKIND), dimension(:), pointer :: gradSSHX, gradSSHY, gradSSHZ
      real (kind=RKIND), dimension(:), pointer :: gradSSHZonal, gradSSHMeridional
      real (kind=RKIND), dimension(:,:), pointer :: surfaceVelocity, SSHGradient

      ! BiCGSTAB Array Pointers
      real (kind=RKIND), dimension(:), pointer :: CGvec_r0,CGvec_r1,CGvec_p0,CGvec_p1,CGvec_h,CGvec_s,CGvec_v0,CGvec_v1,CGvec_t,CGvec_rh
      real (kind=RKIND), dimension(:), pointer :: barotropicCoriolisTerm
      real (kind=RKIND), dimension(:), pointer :: CGvec_u0,CGvec_u1,CGvec_s0,CGvec_s1,CGvec_w0,CGvec_w1
      real (kind=RKIND), dimension(:), pointer :: CGvec_z0,CGvec_z1
      ! BiCGSTAB variables 
      real (kind=RKIND) :: CGcst_rrold,CGcst_rrnew,CGcst_tssum,CGcst_ttsum,CGcst_alpha,CGcst_beta,CGcst_rvsum,CGcst_w0,CGcst_w1
      real (kind=RKIND) :: CGcst_alpha0,CGcst_alpha1,CGcst_beta0,CGcst_beta1
      real (kind=RKIND) :: CGcst_rrold_global,CGcst_rrnew_global,CGcst_tssum_global,CGcst_ttsum_global,CGcst_rvsum_global
      real (kind=RKIND) :: CGcst_wuold,CGcst_wuold_global,CGcst_wunew,CGcst_wunew_global
      real (kind=RKIND) :: CGcst_ruold,CGcst_ruold_global,CGcst_runew,CGcst_runew_global
      real (kind=RKIND), dimension(3) :: CGcst_allreduce,CGcst_allreduce_global
      real (kind=RKIND) :: sshTendA,sshTendB,sshTendC,temp1,temp2,temp3,temp4,crit_out,crit_main,kappa_out,kappa_main,resid

      ! Diagnostics Field Pointers
      type (field2DReal), pointer :: normalizedRelativeVorticityEdgeField, divergenceField, relativeVorticityField
      type (field1DReal), pointer :: barotropicThicknessFluxField, boundaryLayerDepthField, effectiveDensityField
      ! tracer tendencies brought in here to normalize by new layer thickness
      real (kind=RKIND), dimension(:,:,:), pointer :: &
        activeTracerHorizontalAdvectionTendency,      &
        activeTracerVerticalAdvectionTendency,        &
        activeTracerSurfaceFluxTendency,              &
        activeTracerNonLocalTendency

      real (kind=RKIND), dimension(:,:), pointer :: &
        temperatureShortWaveTendency
      ! State/Tend Field Pointers
      type (field1DReal), pointer :: normalBarotropicVelocitySubcycleField, sshSubcycleField
      type (field2DReal), pointer :: highFreqThicknessField, lowFreqDivergenceField
      type (field2DReal), pointer :: normalBaroclinicVelocityField, layerThicknessField
      type (field2DReal), pointer :: normalVelocityField
      type (field3DReal), pointer :: tracersGroupField

      ! tracer iterators
      type (mpas_pool_iterator_type) :: groupItr
      character (len=StrKIND) :: modifiedGroupName
      character (len=StrKIND) :: configName
      integer :: threadNum

      integer :: temp_mask

      dminfo = domain % dminfo

      call mpas_timer_start("se timestep")

      call mpas_pool_get_config(domain % configs, 'config_n_bcl_iter_beg', config_n_bcl_iter_beg)
      call mpas_pool_get_config(domain % configs, 'config_n_bcl_iter_mid', config_n_bcl_iter_mid)
      call mpas_pool_get_config(domain % configs, 'config_n_bcl_iter_end', config_n_bcl_iter_end)
      call mpas_pool_get_config(domain % configs, 'config_n_ts_iter', config_n_ts_iter)
      call mpas_pool_get_config(domain % configs, 'config_btr_subcycle_loop_factor', config_btr_subcycle_loop_factor)
      call mpas_pool_get_config(domain % configs, 'config_btr_gam1_velWt1', config_btr_gam1_velWt1)
      call mpas_pool_get_config(domain % configs, 'config_btr_gam3_velWt2', config_btr_gam3_velWt2)
      call mpas_pool_get_config(domain % configs, 'config_btr_solve_SSH2', config_btr_solve_SSH2)
      call mpas_pool_get_config(domain % configs, 'config_n_btr_cor_iter', config_n_btr_cor_iter)
      call mpas_pool_get_config(domain % configs, 'config_btr_gam2_SSHWt1', config_btr_gam2_SSHWt1)
      call mpas_pool_get_config(domain % configs, 'config_filter_btr_mode', config_filter_btr_mode)

      call mpas_pool_get_config(domain % configs, 'config_mom_del4', config_mom_del4)
      call mpas_pool_get_config(domain % configs, 'config_use_freq_filtered_thickness', config_use_freq_filtered_thickness)
      call mpas_pool_get_config(domain % configs, 'config_time_integrator', config_time_integrator)
      call mpas_pool_get_config(domain % configs, 'config_vel_correction', config_vel_correction)

      call mpas_pool_get_config(domain % configs, 'config_prescribe_velocity', config_prescribe_velocity)
      call mpas_pool_get_config(domain % configs, 'config_prescribe_thickness', config_prescribe_thickness)

      call mpas_pool_get_config(domain % configs, 'config_prescribe_velocity', config_prescribe_velocity)
      call mpas_pool_get_config(domain % configs, 'config_prescribe_thickness', config_prescribe_thickness)

      call mpas_pool_get_config(domain % configs, 'config_use_standardGM', config_use_standardGM)
      call mpas_pool_get_config(domain % configs, 'config_use_cvmix_kpp', config_use_cvmix_kpp)
      call mpas_pool_get_config(domain % configs, 'config_land_ice_flux_mode', config_land_ice_flux_mode)

      call mpas_pool_get_config(domain % configs, 'config_num_halos', config_num_halos)

      call mpas_pool_get_config(domain % configs, 'config_compute_active_tracer_budgets', config_compute_active_tracer_budgets)
      allocate(n_bcl_iter(config_n_ts_iter))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !  Prep variables before first iteration
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call mpas_timer_start("se prep")

      istep = istep + 1

      block => domain % blocklist


      do while (associated(block))
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
         call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
         call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

         call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)

         call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityNew, 2)
         call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew, 2)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)

         call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
         call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)

         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

         call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)

         call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceCur, 1)
         call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceNew, 2)

         call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

         nCells = nCellsPtr
         nEdges = nEdgesPtr
          
         ! Initialize * variables that are used to compute baroclinic tendencies below.

         !$omp do schedule(runtime) private(k)
         do iEdge = 1, nEdges
            do k = 1, nVertLevels !maxLevelEdgeTop % array(iEdge)

               ! The baroclinic velocity needs be recomputed at the beginning of a
               ! timestep because the implicit vertical mixing is conducted on the
               ! total u.  We keep normalBarotropicVelocity from the previous timestep.
               ! Note that normalBaroclinicVelocity may now include a barotropic component, because the
               ! weights layerThickness have changed.  That is OK, because the barotropicForcing variable
               ! subtracts out the barotropic component from the baroclinic.
               normalBaroclinicVelocityCur(k,iEdge) = normalVelocityCur(k,iEdge) - normalBarotropicVelocityCur(iEdge)

               normalVelocityNew(k,iEdge) = normalVelocityCur(k,iEdge)

               normalBaroclinicVelocityNew(k,iEdge) = normalBaroclinicVelocityCur(k,iEdge)
            end do
         end do
         !$omp end do

         !$omp do schedule(runtime) private(k)
         do iCell = 1, nCells
            sshNew(iCell) = sshCur(iCell)
            do k = 1, maxLevelCell(iCell)
               layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell)
               ! set vertAleTransportTop to zero for stage 1 velocity tendency, first time through.
               vertAleTransportTop(k,iCell) = 0.0_RKIND
            end do
         end do
         !$omp end do

         call mpas_pool_begin_iteration(tracersPool)
         do while ( mpas_pool_get_next_member(tracersPool, groupItr))
            if ( groupItr % memberType == MPAS_POOL_FIELD ) then
               call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupCur, 1)
               call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupNew, 2)

               if ( associated(tracersGroupCur) .and. associated(tracersGroupNew) ) then
                  !$omp do schedule(runtime) private(k)
                  do iCell = 1, nCells
                     do k = 1, maxLevelCell(iCell)
                        tracersGroupNew(:,k,iCell) = tracersGroupCur(:,k,iCell)
                     end do
                  end do
                  !$omp end do
               end if
            end if
         end do


         if (associated(highFreqThicknessNew)) then
            !$omp do schedule(runtime)
            do iCell = 1, nCells
               highFreqThicknessNew(:, iCell) = highFreqThicknessCur(:, iCell)
            end do
            !$omp end do
         end if

         if (associated(lowFreqDivergenceNew)) then
            !$omp do schedule(runtime)
            do iCell = 1, nCells
               lowFreqDivergenceNew(:, iCell) = lowFreqDivergenceCur(:, iCell)
            end do
            !$omp end do
         endif

         block => block % next
      end do


      call mpas_timer_stop("se prep")
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN large iteration loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      n_bcl_iter = config_n_bcl_iter_mid
      n_bcl_iter(1) = config_n_bcl_iter_beg
      n_bcl_iter(config_n_ts_iter) = config_n_bcl_iter_end

      do split_explicit_step = 1, config_n_ts_iter
         call mpas_timer_start('se loop')

         stage1_tend_time = min(split_explicit_step,2)

         call mpas_pool_get_subpool(domain % blocklist % structs, 'diagnostics', diagnosticsPool)

         call mpas_threading_barrier()
         ! ---  update halos for diagnostic ocean boundayr layer depth
         if (config_use_cvmix_kpp) then
            call mpas_timer_start("se halo diag obd")
            call mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
            call mpas_timer_stop("se halo diag obd")
         end if

         ! ---  update halos for diagnostic variables
         call mpas_timer_start("se halo diag")

         call mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
         if (config_mom_del4 > 0.0_RKIND) then
           call mpas_dmpar_field_halo_exch(domain, 'divergence')
           call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
         end if
         call mpas_timer_stop("se halo diag")

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         !  Stage 1: Baroclinic velocity (3D) prediction, explicit with long timestep
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (config_use_freq_filtered_thickness) then
            call mpas_timer_start("se freq-filtered-thick computations")

            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
               call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
               call mpas_pool_get_subpool(block % structs, 'state', statepool)
               call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
               call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

               call ocn_tend_freq_filtered_thickness(tendPool, statePool, diagnosticsPool, meshPool, stage1_tend_time)
               block => block % next
            end do
            call mpas_timer_stop("se freq-filtered-thick computations")

            call mpas_threading_barrier()

            call mpas_timer_start("se freq-filtered-thick halo update")

            call mpas_dmpar_field_halo_exch(domain, 'tendHighFreqThickness')
            call mpas_dmpar_field_halo_exch(domain, 'tendLowFreqDivergence')

            call mpas_timer_stop("se freq-filtered-thick halo update")
            call mpas_threading_barrier()

            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
               call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)

               call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
               call mpas_pool_get_subpool(block % structs, 'state', statePool)
               call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
               call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

               call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

               call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
               call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)

               call mpas_pool_get_array(tendPool, 'highFreqThickness', highFreqThicknessTend)

               nCells = nCellsPtr

               !$omp do schedule(runtime) private(k)
               do iCell = 1, nCells
                  do k = 1, maxLevelCell(iCell)
                     ! this is h^{hf}_{n+1}
                     highFreqThicknessNew(k,iCell) = highFreqThicknessCur(k,iCell) + dt * highFreqThicknessTend(k,iCell)
                  end do
               end do
               !$omp end do

               block => block % next
            end do

         endif

         ! compute velocity tendencies, T(u*,w*,p*)
         call mpas_timer_start("se bcl vel")

         call mpas_timer_start('se bcl vel tend')
         block => domain % blocklist
         do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
           call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
           call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
           call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
           call mpas_pool_get_subpool(block % structs, 'state', statePool)
           call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
           call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
           call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
           call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)

           call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
           call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, stage1_tend_time)
           call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)

           call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)

           call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)

           call ocn_tend_vel(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, scratchPool, stage1_tend_time)

           block => block % next
         end do
         call mpas_timer_stop('se bcl vel tend')

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! BEGIN baroclinic iterations on linear Coriolis term
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do j=1,n_bcl_iter(split_explicit_step)

            ! Use this G coefficient to avoid an if statement within the iEdge loop.
            if (trim(config_time_integrator) == 'unsplit_explicit') then
               split = 0
            elseif (trim(config_time_integrator) == 'split_explicit') then
               split = 1
            endif

            call mpas_timer_start('bcl iters on linear Coriolis')
            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
               call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

               call mpas_pool_get_subpool(block % structs, 'state', statePool)
               call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
               call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
               call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

               call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
               call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
               call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

               call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
               call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityCur, 1)
               call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityNew, 2)
               call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)

               call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)

               call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
               call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)

               ! Only need to loop over the 1 halo, since there is a halo exchange immediately after this computation.
               nEdges = nEdgesArray( 1 )

               ! Put f*normalBaroclinicVelocity^{perp} in normalVelocityNew as a work variable
               call ocn_fuperp(statePool, meshPool, 2)

               allocate(uTemp(nVertLevels))

               !$omp do schedule(runtime) private(cell1, cell2, k, normalThicknessFluxSum, thicknessSum)
               do iEdge = 1, nEdges
                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)

                  uTemp = 0.0_RKIND  ! could put this after with uTemp(maxleveledgetop+1:nvertlevels)=0
                  do k = 1, maxLevelEdgeTop(iEdge)

                     ! normalBaroclinicVelocityNew = normalBaroclinicVelocityOld + dt*(-f*normalBaroclinicVelocityPerp
                     !                             + T(u*,w*,p*) + g*grad(SSH*) )
                     ! Here uNew is a work variable containing -fEdge(iEdge)*normalBaroclinicVelocityPerp(k,iEdge)
                      uTemp(k) = normalBaroclinicVelocityCur(k,iEdge) &
                         + dt * (normalVelocityTend(k,iEdge) &
                         + normalVelocityNew(k,iEdge) &  ! this is f*normalBaroclinicVelocity^{perp}
                         + split * gravity * (  sshNew(cell2) - sshNew(cell1) ) &
                          / dcEdge(iEdge) )
                  enddo

                  ! thicknessSum is initialized outside the loop because on land boundaries
                  ! maxLevelEdgeTop=0, but I want to initialize thicknessSum with a
                  ! nonzero value to avoid a NaN.
                  normalThicknessFluxSum = layerThicknessEdge(1,iEdge) * uTemp(1)
                  thicknessSum  = layerThicknessEdge(1,iEdge)

                  do k = 2, maxLevelEdgeTop(iEdge)
                     normalThicknessFluxSum = normalThicknessFluxSum + layerThicknessEdge(k,iEdge) * uTemp(k)
                     thicknessSum  =  thicknessSum + layerThicknessEdge(k,iEdge)
                  enddo
                  barotropicForcing(iEdge) = split * normalThicknessFluxSum / thicknessSum / dt


                  do k = 1, maxLevelEdgeTop(iEdge)
                     ! These two steps are together here:
                     !{\bf u}'_{k,n+1} = {\bf u}'_{k,n} - \Delta t {\overline {\bf G}}
                     !{\bf u}'_{k,n+1/2} = \frac{1}{2}\left({\bf u}^{'}_{k,n} +{\bf u}'_{k,n+1}\right)
                     ! so that normalBaroclinicVelocityNew is at time n+1/2
                     normalBaroclinicVelocityNew(k,iEdge) = 0.5_RKIND*( &
                       normalBaroclinicVelocityCur(k,iEdge) + uTemp(k) - dt * barotropicForcing(iEdge))
                  enddo
               enddo ! iEdge
               !$omp end do

               deallocate(uTemp)

               block => block % next
            end do

            call mpas_threading_barrier()

            call mpas_timer_start("se halo normalBaroclinicVelocity")
            call mpas_dmpar_field_halo_exch(domain, 'normalBaroclinicVelocity', timeLevel=2)
            call mpas_timer_stop("se halo normalBaroclinicVelocity")

            call mpas_timer_stop('bcl iters on linear Coriolis')

         end do  ! do j=1,config_n_bcl_iter

         call mpas_timer_start('se halo barotropicForcing')
         call mpas_dmpar_field_halo_exch(domain, 'barotropicForcing')
         call mpas_timer_stop('se halo barotropicForcing')
         call mpas_timer_stop("se bcl vel")
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! END baroclinic iterations on linear Coriolis term
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!****************************************************************************************         
!****************************************************************************************         
!****************************************************************************************         
!****************************************************************************************         
!****************************************************************************************         


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         !  Stage 2: Barotropic velocity (2D) prediction, explicitly subcycled
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! Stage 2.1 : Preparation of varibales before the two outer iteration
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
         call mpas_timer_start("stage 2 timer")

         call mpas_timer_start("si btr prep")

         oldBtrSubcycleTime = 1
         newBtrSubcycleTime = 2

         cellHaloComputeCounter = config_num_halos
         edgeHaloComputeCounter = config_num_halos + 1

         kappa_out  = 1.0D-8
         kappa_main = 1.0D-15
         crit_out   = kappa_out  * area_mean
         crit_main  = kappa_main * area_mean

      
            ! Initialize variables for barotropic subcycling
            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
               call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
               call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

               call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
               call mpas_pool_get_subpool(block % structs, 'state', statePool)
               call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
               call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

                call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
                call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
                call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
                call mpas_pool_get_array(meshPool, 'fEdge' , fEdge)
                call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
                call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

               call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
               call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', barotropicThicknessFlux)
               call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm',barotropicCoriolisTerm)

               call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
               call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, 1)
               call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, 2)
               call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, 1)
               call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, 2)
               call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur, 1)
               call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew, 2)

	
               nCells = nCellsPtr
               nEdges = nEdgesPtr

               if (config_filter_btr_mode) then
                  !$omp do schedule(runtime)
                  do iEdge = 1, nEdges
                     barotropicForcing(iEdge) = 0.0_RKIND
                  end do
                  !$omp end do
               endif

               !$omp do schedule(runtime)
               do iCell = 1, nCells
                  ! sshSubcycleOld = sshOld
                  sshSubcycleCur(iCell) = sshCur(iCell)
               end do
               !$omp end do

               !$omp do schedule(runtime)
               do iEdge = 1, nEdges

                  ! normalBarotropicVelocitySubcycleOld = normalBarotropicVelocityOld
                  normalBarotropicVelocitySubcycleCur(iEdge) = normalBarotropicVelocityCur(iEdge)

                  ! normalBarotropicVelocitySubcycleNew = normalBarotropicVelocityOld
!                 normalBarotropicVelocitySubcycleNew(iEdge) = normalBarotropicVelocityCur(iEdge)

                  ! normalBarotropicVelocityNew = BtrOld  This is the first for the summation
!                         normalBarotropicVelocityNew(iEdge) = normalBarotropicVelocityCur(iEdge)

                  ! barotropicThicknessFlux = 0
                              barotropicThicknessFlux(iEdge) = 0.0_RKIND
               end do
               !$omp end do

               do iEdge = 1, nEdges

                          ! Compute the barotropic Coriolis term, -f*uPerp
                          barotropicCoriolisTerm(iEdge) = 0.0_RKIND
                          do i = 1, nEdgesOnEdge(iEdge)
                             eoe = edgesOnEdge(i,iEdge)
                             barotropicCoriolisTerm(iEdge) = barotropicCoriolisTerm(iEdge) + weightsOnEdge(i,iEdge) &
                                                           * normalBarotropicVelocityCur(eoe) * fEdge(eoe)
                          end do ! i
                end do ! iEdge

               block => block % next
            end do  ! block

            call mpas_timer_stop("si btr prep")

            call mpas_threading_barrier()

            call mpas_timer_start("si halo btrCoriolis")
            call mpas_dmpar_exch_group_create(domain, subcycleGroupName)
            call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'barotropicCoriolisTerm', 1 )
            call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'normalBarotropicVelocitySubcycle', 1 )
            call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'sshSubcycle', 1 )
            call mpas_threading_barrier()
            call mpas_dmpar_exch_group_full_halo_exch(domain, subcycleGroupName)
            call mpas_dmpar_exch_group_destroy(domain, subcycleGroupName)
            call mpas_timer_stop("si halo btrCoriolis")


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! Stage 2.2 : Initial guess, compute r0 = b - A * x0
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! Barotropic iteration : Compute r0 = b - Ax0
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              call mpas_timer_start("si btr first r0")

              block => domain % blocklist
              do while (associated(block))
                call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
                call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
                call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
                call mpas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)
                call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
                call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

                call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
                call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
                call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                call mpas_pool_get_subpool(block % structs, 'state', statePool)
                call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
                call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

                call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
                call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
                call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
                call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
                call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
                call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
                call mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell', refBottomDepthTopOfCell)
                call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
                call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
                call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

                call mpas_pool_get_array(statePool, 'ssh', sshCur,1)
                call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, 1)
                call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, 2)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur,1)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, 1)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, 2)

                call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
                call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm',barotropicCoriolisTerm)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0', CGvec_r0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_p0', CGvec_p0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_v0', CGvec_v0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_u0', CGvec_u0)

                call mpas_pool_get_array(diagnosticsPool, 'CGvec_z0', CGvec_z0)

                nCells = nCellsPtr
                nEdges = nEdgesPtr

                nCells = nCellsArray( 1 )
                nEdges = nEdgesArray( 2 )

!--- Make preconditioner at the first time stepping only ---------------------------!
              if ( istep == 1 ) then

                do iCell = 1, nCells

                  CGvec_z0(iCell) = 0.0_RKIND

                  do i = 1, nEdgesOnCell(iCell)
                    iEdge = edgesOnCell(i, iCell)

                    cell1 = cellsOnEdge(1, iEdge)
                    cell2 = cellsOnEdge(2, iEdge)

                    ! method 1, matches method 0 without pbcs, works with pbcs.
                    thicknessSumCur = min(bottomDepth(cell1), bottomDepth(cell2))

                    !--------------------------------------------------------------!
                    fluxAx = edgeSignOnCell(i,iCell)*dvEdge(iEdge)*thicknessSumCur / dcEdge(iEdge)

                    if (cell1 == iCell) then
                      CGvec_z0(iCell) = CGvec_z0(iCell) + fluxAx  ! reversed sign
                    elseif ( cell2 == iCell) then
                      CGvec_z0(iCell) = CGvec_z0(iCell) - fluxAx  ! reversed sign
                    endif
                    !--------------------------------------------------------------!

                  end do ! i

                  CGvec_z0(iCell) = CGvec_z0(iCell) - (4.0_RKIND/(gravity*dt**2.0))*areaCell(iCell)

                end do ! iCell

              endif ! istep
!-----------------------------------------------------------------------------------!

                CGcst_rrold = 0.0_RKIND
                CGcst_ruold = 0.0_RKIND

                do iCell = 1, nCells
                  sshTendb1 = 0.0_RKIND
                  sshTendb2 = 0.0_RKIND
                  sshTendAx = 0.0_RKIND
                  do i = 1, nEdgesOnCell(iCell)
                    iEdge = edgesOnCell(i, iCell)

                    cell1 = cellsOnEdge(1, iEdge)
                    cell2 = cellsOnEdge(2, iEdge)

                    ! Interpolation sshEdge
                    sshEdgeCur = 0.5_RKIND * (sshCur(cell1) + sshCur(cell2))

                    ! method 1, matches method 0 without pbcs, works with pbcs.
                    thicknessSumCur = sshEdgeCur + min(bottomDepth(cell1), bottomDepth(cell2))

                    ! nabla (ssh^0)
                    sshDiffCur = (  sshCur(cell2) -   sshCur(cell1)) / dcEdge(iEdge)

                    !--------------------------------------------------------------!
                       fluxb1 = thicknessSumCur * normalBarotropicVelocityCur(iEdge)
                       fluxb2 = thicknessSumCur * (0.5_RKIND*gravity*sshDiffCur + (-barotropicCoriolisTerm(iEdge)-barotropicForcing(iEdge)) )
                       fluxAx = thicknessSumCur * sshDiffCur
       
                    sshTendb1 = sshTendb1 + edgeSignOnCell(i, iCell) * fluxb1 * dvEdge(iEdge)
                    sshTendb2 = sshTendb2 + edgeSignOnCell(i, iCell) * fluxb2 * dvEdge(iEdge)
                    sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) * fluxAx * dvEdge(iEdge)
                    !--------------------------------------------------------------!

                  end do ! i

                    sshTendb1 = (4.0_RKIND/(gravity*dt)) * sshTendb1
                    sshTendb2 = (2.0_RKIND/(gravity   )) * sshTendb2

                    sshCurArea = (4.0_RKIND/(gravity*dt**2.0)) *   sshCur(iCell) * areaCell(iCell)

                    CGvec_r0(iCell) = (-sshCurArea - sshTendb1 + sshTendb2)   &
                                     -(-sshCurArea - sshTendAx) 

                  ! Initialize p0
                    CGvec_p0(iCell) = CGvec_r0(iCell)
                    CGvec_u0(iCell) = CGvec_r0(iCell) / CGvec_z0(iCell)

                    CGcst_ruold = CGcst_ruold + CGvec_r0(iCell) * CGvec_u0(iCell)

                end do ! iCell

                block => block % next
              end do  ! block

              
              ! Global sum across CPUs
              call mpas_dmpar_sum_real(dminfo,CGcst_ruold,CGcst_ruold_global)

              call mpas_threading_barrier()
              call mpas_timer_start("si halo iter r0")
              call mpas_dmpar_exch_group_create(domain, subcycleGroupName)
              call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'CGvec_u0', 1)
              call mpas_threading_barrier()
              call mpas_dmpar_exch_group_full_halo_exch(domain, subcycleGroupName)
              call mpas_dmpar_exch_group_destroy(domain, subcycleGroupName)
              call mpas_timer_stop("si halo iter r0")


              !--- A * u0 = A * r0---------------------------------------------------------------------------!

              block => domain % blocklist
              do while (associated(block))
                call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
                call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
                call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
                call mpas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)
                call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
                call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

                call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
                call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
                call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                call mpas_pool_get_subpool(block % structs, 'state', statePool)
                call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
                call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

                call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
                call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
                call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
                call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
                call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
                call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
                call mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell', refBottomDepthTopOfCell)
                call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
                call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
                call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
                     call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

                call mpas_pool_get_array(statePool, 'ssh', sshCur,1)
                call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, 1)
                call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, 2)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur,1)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, 1)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, 2)

                call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
                call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm',barotropicCoriolisTerm)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0', CGvec_r0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_p0', CGvec_p0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_v0', CGvec_v0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_u0', CGvec_u0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_w0', CGvec_w0)


                nCells = nCellsPtr
                nEdges = nEdgesPtr

                nCells = nCellsArray( 1 )
                nEdges = nEdgesArray( 2 )

                CGcst_wuold = 0.0_RKIND

                 
                do iCell = 1, nCells

                  sshTendAx = 0.0_RKIND

                  do i = 1, nEdgesOnCell(iCell)
                    iEdge = edgesOnCell(i, iCell)

                    cell1 = cellsOnEdge(1, iEdge)
                    cell2 = cellsOnEdge(2, iEdge)

                    ! Interpolation sshEdge
                    sshEdgeCur = 0.5_RKIND * (sshCur(cell1) + sshCur(cell2))

                    ! method 1, matches method 0 without pbcs, works with pbcs.
                    thicknessSumCur = sshEdgeCur + min(bottomDepth(cell1), bottomDepth(cell2))

                    ! nabla (ssh^0)
                    sshDiffCur = (CGvec_u0(cell2)- CGvec_u0(cell1)) / dcEdge(iEdge)

                    !--------------------------------------------------------------!
                       fluxAx = thicknessSumCur * sshDiffCur
       
                    sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) * fluxAx * dvEdge(iEdge)
                    !--------------------------------------------------------------!

                  end do ! i

                    sshCurArea = (4.0_RKIND/(gravity*dt**2.0)) * CGvec_u0(iCell) * areaCell(iCell)
                     
                    CGvec_w0(iCell) = -sshCurArea - sshTendAx 

                    CGcst_wuold = CGcst_wuold + CGvec_w0(iCell) * CGvec_u0(iCell)

                end do ! iCell


              ! Global sum across CPUs
              call mpas_dmpar_sum_real(dminfo,CGcst_wuold,CGcst_wuold_global)

                CGcst_alpha0 = CGcst_ruold_global / CGcst_wuold_global

                CGcst_beta0  = 0.0_RKIND

                block => block % next
              end do  ! block

              call mpas_timer_stop("si btr first r0")



         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! Stage 2.3 : Outer iterations - lagged values are sufficiently up to date
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


              call mpas_timer_start("si btr outer iteration")

              iter = 0
              resid = (crit_out+100.0)**2.0

              !**************************************************************!
              do while ( dsqrt(resid) > crit_out ) 
              !**************************************************************!


              block => domain % blocklist
              do while (associated(block))
                call mpas_pool_get_dimension(block % dimensions, 'nCells'     , nCellsPtr  )
                call mpas_pool_get_dimension(block % dimensions, 'nEdges'     , nEdgesPtr  )
                call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
                call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

                call mpas_pool_get_subpool(block % structs, 'tend'       , tendPool       )
                call mpas_pool_get_subpool(block % structs, 'mesh'       , meshPool       )
                call mpas_pool_get_subpool(block % structs, 'state'      , statePool      )
                call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

                call mpas_pool_get_subpool(statePool, 'tracers'    , tracersPool    )

                call mpas_pool_get_subpool(tendPool , 'tracersTend', tracersTendPool)

                call mpas_pool_get_array(meshPool, 'nEdgesOnCell',            nEdgesOnCell           )
                call mpas_pool_get_array(meshPool, 'edgesOnCell',             edgesOnCell            )
                call mpas_pool_get_array(meshPool, 'cellsOnEdge',             cellsOnEdge            )
                call mpas_pool_get_array(meshPool, 'dcEdge',                  dcEdge                 )
                call mpas_pool_get_array(meshPool, 'bottomDepth',             bottomDepth            )
                call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop',         maxLevelEdgeTop        )
                call mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell', refBottomDepthTopOfCell)
                call mpas_pool_get_array(meshPool, 'edgeSignOnCell',          edgeSignOnCell         )
                call mpas_pool_get_array(meshPool, 'dvEdge',                  dvEdge                 )
                call mpas_pool_get_array(meshPool, 'areaCell',                areaCell               )
                call mpas_pool_get_array(meshPool, 'nEdgesOnEdge',            nEdgesOnEdge           )
                call mpas_pool_get_array(meshPool, 'edgesOnEdge',             edgesOnEdge            )
                call mpas_pool_get_array(meshPool, 'weightsOnEdge',           weightsOnEdge          )
                call mpas_pool_get_array(meshPool, 'fEdge',                   fEdge                  )
                call mpas_pool_get_array(meshPool, 'edgeMask',                edgeMask               )

                call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
                call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)
                call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, 1)
                call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, 2)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur,1)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew,2)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, 1)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, 2)

                call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing',      barotropicForcing     )
                call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm', barotropicCoriolisTerm)

                call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0', CGvec_r0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_p0', CGvec_p0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_v0', CGvec_v0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_r1', CGvec_r1)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_p1', CGvec_p1)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_t' , CGvec_t )
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_h' , CGvec_h )

                call mpas_pool_get_array(diagnosticsPool, 'CGvec_s0' , CGvec_s0 )
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_u0' , CGvec_u0 )
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_w0' , CGvec_w0 )
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_s1' , CGvec_s1 )
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_u1' , CGvec_u1 )
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_w1' , CGvec_w1 )

                call mpas_pool_get_array(diagnosticsPool, 'CGvec_z0', CGvec_z0)

                call mpas_pool_get_array(tendPool, 'ssh', sshTend)

                nCells = nCellsPtr
                nEdges = nEdgesPtr

                nCells = nCellsArray( 1 )
                nEdges = nEdgesArray( 2 )

                iter = iter + 1

! 1. p1 = u0 + beta0 * p0
                do iCell = 1, nCells
                  CGvec_p1(iCell) = CGvec_u0(iCell) + CGcst_beta0 * CGvec_p0(iCell)
                end do ! iCell

! 2. s1 = w0 + beta0 * s0
                do iCell = 1, nCells
                  CGvec_s1(iCell) = CGvec_w0(iCell) + CGcst_beta0 * CGvec_s0(iCell)
                end do ! iCell

! 3. x1 = x0 + alpha0 * p1
                do iCell = 1,nCells
 		  sshSubcycleNew(iCell) = sshSubcycleCur(iCell) + CGcst_alpha0 * CGvec_p1(iCell)
              	end do ! iCell

! 4. r1 = r0 - alpha * s1
                do iCell = 1,nCells
 		  CGvec_r1(iCell) = CGvec_r0(iCell) - CGcst_alpha0 * CGvec_s1(iCell) 
              	end do ! iCell

! 5. u1 = r1
                do iCell = 1,nCells
 		  CGvec_u1(iCell) = CGvec_r1(iCell) / CGvec_z0(iCell)
              	end do ! iCell

               call mpas_timer_start("si halo iter")
               call mpas_dmpar_exch_group_create(domain, subcycleGroupName)
               call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'CGvec_u1', 1)
               call mpas_threading_barrier()
               call mpas_dmpar_exch_group_full_halo_exch(domain, subcycleGroupName)
               call mpas_dmpar_exch_group_destroy(domain, subcycleGroupName)
               call mpas_timer_stop("si halo iter")

               call mpas_threading_barrier()

                !------------------------------------------------------------------------------!

! 6. w1 = A * u1


                do iCell = 1,nCells
                   sshTend(iCell) = 0.0_RKIND
                end do

                do iEdge = 1,nEdges
                  cell1 = cellsOnEdge(1, iEdge)
                  cell2 = cellsOnEdge(2, iEdge)

                    ! Interpolation sshEdge
                    sshEdgeLag = 0.5_RKIND * (sshSubcycleCur(cell1) + sshSubcycleCur(cell2))

                    thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1), bottomDepth(cell2))

                    ! nabla (ssh^0)
                    sshDiffNew = (CGvec_u1(cell2)-CGvec_u1(cell1)) / dcEdge(iEdge)

                        fluxAx = thicknessSumLag * sshDiffNew
       
                    sshTend(cell1) = sshTend(cell1) - fluxAx * dvEdge(iEdge)
                    sshTend(cell2) = sshTend(cell2) + fluxAx * dvEdge(iEdge)
                end do          

                CGcst_runew = 0.0_RKIND
                CGcst_wunew = 0.0_RKIND
                CGcst_rrnew = 0.0_RKIND

                do iCell = 1,nCells
                     sshLagArea = (4.0_RKIND/(gravity*dt**2.0)) * CGvec_u1(iCell) * areaCell(iCell)

                  CGvec_w1(iCell) = -sshLagArea - sshTend(iCell)

                  CGcst_runew = CGcst_runew + CGvec_r1(iCell) * CGvec_u1(iCell)
                  CGcst_wunew = CGcst_wunew + CGvec_w1(iCell) * CGvec_u1(iCell)
                  CGcst_rrnew = CGcst_rrnew + CGvec_r1(iCell) * CGvec_r1(iCell)
                end do

                !------------------------------------------------------------------------------!
       
                  CGcst_allreduce(1) = CGcst_runew
                  CGcst_allreduce(2) = CGcst_wunew
                  CGcst_allreduce(3) = CGcst_rrnew

                ! Global sum across CPUs
                call mpas_dmpar_sum_real_array(dminfo, 3, CGcst_allreduce, CGcst_allreduce_global)

                  CGcst_runew_global = CGcst_allreduce_global(1)
                  CGcst_wunew_global = CGcst_allreduce_global(2)
                               resid = CGcst_allreduce_global(3)

! 7. beta1 = runew / ruold
               CGcst_beta1 = CGcst_runew_global / CGcst_ruold_global

! 8. alpha1 = 1 / (wunew / runew - beta1 / alpha0 ) 
               CGcst_alpha1 = 1.0_RKIND / ( CGcst_wunew_global/CGcst_runew_global     &
                                           -CGcst_beta1/CGcst_alpha0)

! 14. Stepping for next iteration
                do iCell = 1,nCells
                  CGvec_r0(iCell) = CGvec_r1(iCell)               
                  CGvec_p0(iCell) = CGvec_p1(iCell)               
                  CGvec_s0(iCell) = CGvec_s1(iCell)               
                  CGvec_u0(iCell) = CGvec_u1(iCell)               
                  CGvec_w0(iCell) = CGvec_w1(iCell)               

                  temp1 = sshCur(iCell)
                  temp2 = sshSubcycleNew(iCell)

                  sshSubcycleCur(iCell) = sshSubcycleNew(iCell)               
                   
                  ! sshSubcycleNew => sshSubcycleNew^* for u update
                    CGvec_h(iCell) = 0.5_RKIND * (temp1+temp2)
              	end do ! iCell
		
  		  CGcst_rrold_global = CGcst_rrnew_global		
  		  CGcst_ruold_global = CGcst_runew_global		
  		  CGcst_wuold_global = CGcst_wunew_global		
  		  CGcst_beta0        = CGcst_beta1       		
  		  CGcst_alpha0       = CGcst_alpha1

                block => block % next
              end do  ! block

               !   boundary update on SSHnew
               call mpas_timer_start("si halo iter")
               call mpas_dmpar_exch_group_create(domain, subcycleGroupName)
               call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'sshSubcycle', timeLevel=1)
               call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'CGvec_h', 1)
               call mpas_threading_barrier()
               call mpas_dmpar_exch_group_full_halo_exch(domain, subcycleGroupName)
               call mpas_dmpar_exch_group_destroy(domain, subcycleGroupName)
               call mpas_timer_stop("si halo iter")


              ! Update normalBarotropicVelocitySubcycle
              
                  block => domain % blocklist
                  do while (associated(block))
                     call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
                     call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

                     call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                     call mpas_pool_get_subpool(block % structs, 'state', statePool)
                     call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
                     call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

                     call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
                     call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
                     call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
                     call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
                     call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
                     call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
                     call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

                     call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, 1)
                     call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, 2)
                     call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur,1)
                     call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew,2)

                     call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
                     call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm', barotropicCoriolisTerm)
                     call mpas_pool_get_array(diagnosticsPool, 'CGvec_h', CGvec_h)

                     nEdges = nEdgesPtr
                     nEdges = nEdgesArray( edgeHaloComputeCounter )
!                    nEdges = nEdgesArray( 1)

                     !$omp do schedule(runtime) private(cell1, cell2, CoriolisTerm, i, eoe)
                     do iEdge = 1, nEdges

                        temp_mask = edgeMask(1, iEdge)

                          cell1 = cellsOnEdge(1,iEdge)
                          cell2 = cellsOnEdge(2,iEdge)

                          ! Compute the barotropic Coriolis term, -f*uPerp
                          CoriolisTerm = 0.0_RKIND
                                 temp1 = 0.0_RKIND
                          do i = 1, nEdgesOnEdge(iEdge)
                             eoe = edgesOnEdge(i,iEdge)

                             CoriolisTerm = CoriolisTerm + weightsOnEdge(i,iEdge)  &
                                          * normalBarotropicVelocityCur(eoe) * fEdge(eoe)
                                    temp1 =        temp1 + weightsOnEdge(i,iEdge)  &
                                          * normalBarotropicVelocitySubcycleCur(eoe) * fEdge(eoe)
                          end do

                          barotropicCoriolisTerm(iEdge) = 0.5_RKIND * (CoriolisTerm+temp1)

                          normalBarotropicVelocitySubcycleNew(iEdge) &
                            = temp_mask &
                            * (normalBarotropicVelocityCur(iEdge) &
                            - dt * (-barotropicCoriolisTerm(iEdge) + gravity &
                            * (CGvec_h(cell2) - CGvec_h(cell1)) / dcEdge(iEdge) - barotropicForcing(iEdge)))

                     end do
                     !$omp end do

                     do iEdge = 1, nEdges
                          ! Cur = New for next iteration
                          normalBarotropicVelocitySubcycleCur(iEdge) = normalBarotropicVelocitySubcycleNew(iEdge) 
                     end do

                     block => block % next
                  end do  ! block


               !   boundary update on SSHnew
!              call mpas_timer_start("si halo iter")
!              call mpas_dmpar_exch_group_create(domain, subcycleGroupName)
!              call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'normalBarotropicVelocitySubcycle', 1)
!              call mpas_threading_barrier()
!              call mpas_dmpar_exch_group_full_halo_exch(domain, subcycleGroupName)
!              call mpas_dmpar_exch_group_destroy(domain, subcycleGroupName)
!              call mpas_timer_stop("si halo iter")


                if ( iter > total_num_cells ) then
                  print*,'******************************************************'
                  print*,'******************************************************'
                  print*,'Iteration number exceeds Max. #iteration: PROGRAM STOP'
                  print*,'******************************************************'
                  print*,'******************************************************'
                  stop
                endif

              !**************************************************************!
              end do ! do while iter
              !**************************************************************!

              call mpas_timer_stop("si btr outer iteration")

              ! Update normalBarotropicVelocitySubcycle
              
                  block => domain % blocklist
                  do while (associated(block))
                     call mpas_pool_get_dimension(block % dimensions, 'nCells'     , nCellsPtr  )
                     call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
                     call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
                     call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

                     call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                     call mpas_pool_get_subpool(block % structs, 'state', statePool)
                     call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
                     call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
                     call mpas_pool_get_subpool(block % structs, 'tend'       , tendPool       )

                     call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
                     call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
                     call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
                     call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
                     call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
                     call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
                     call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

                     call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, 1)
                     call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, 2)
                     call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur,1)
                     call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew,2)
                     call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, 1)
                     call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, 2)
                     call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
                     call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)

                     call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
                     call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm', barotropicCoriolisTerm)
                     call mpas_pool_get_array(diagnosticsPool, 'CGvec_h', CGvec_h)

                     nCells = nCellsPtr
                     nEdges = nEdgesPtr

!                    nCells = nCellsArray( 1)
!                    nEdges = nEdgesArray( 1)

!                 print*,'s_sqrt_outer',dsqrt(resid),iter,split_explicit_step
                     call mpas_log_write( 'iter $i', intArgs=(/ iter /) )

                  ! sshCur = h

                  do iCell = 1, nCells
                    sshSubcycleCur(iCell) = sshSubcycleNew(iCell)
                  end do ! iCell

                      !--- Update u_lag ------------------------------------------------------!

                     !$omp do schedule(runtime) private(cell1, cell2, CoriolisTerm, i, eoe)
                     do iEdge = 1, nEdges

                      temp_mask = edgeMask(1, iEdge)
                          cell1 = cellsOnEdge(1,iEdge)
                          cell2 = cellsOnEdge(2,iEdge)

                          ! Compute the barotropic Coriolis term, -f*uPerp
                          CoriolisTerm = 0.0_RKIND
                                 temp1 = 0.0_RKIND
                          do i = 1, nEdgesOnEdge(iEdge)
                             eoe = edgesOnEdge(i,iEdge)

                             CoriolisTerm = CoriolisTerm + weightsOnEdge(i,iEdge)  &
                                          * normalBarotropicVelocityCur(eoe) * fEdge(eoe)
                                    temp1 =        temp1 + weightsOnEdge(i,iEdge)  &
                                          * normalBarotropicVelocitySubcycleCur(eoe) * fEdge(eoe)
                          end do

                          barotropicCoriolisTerm(iEdge) = 0.5_RKIND * (CoriolisTerm+temp1)

                     end do
                     !$omp end do

                     do iEdge = 1, nEdges
                          ! Cur = New for next iteration
                          normalBarotropicVelocitySubcycleCur(iEdge) = normalBarotropicVelocitySubcycleNew(iEdge) 
                     end do

                     block => block % next
                  end do  ! block

               !   boundary update on SSHnew
               call mpas_timer_start("si halo iter")
               call mpas_dmpar_exch_group_create(domain, subcycleGroupName)
               call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'normalBarotropicVelocitySubcycle', 1)
               call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'barotropicCoriolisTerm', 1 )
               call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'sshSubcycle', 1 )
               call mpas_threading_barrier()
               call mpas_dmpar_exch_group_full_halo_exch(domain, subcycleGroupName)

               call mpas_dmpar_exch_group_destroy(domain, subcycleGroupName)
               call mpas_timer_stop("si halo iter")



         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! Stage 2.4 : Initial guess using lagged values
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! Barotropic iteration : Compute r0 = b - Ax0
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              call mpas_timer_start("si btr first r0")

              block => domain % blocklist
              do while (associated(block))
                call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
                call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
                call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
                call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

                call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
                call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
                call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                call mpas_pool_get_subpool(block % structs, 'state', statePool)
                call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
                call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

                call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
                call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
                call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

                call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

                call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
                call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
                call mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell', refBottomDepthTopOfCell)
                call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
                call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
                call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

                call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
                call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, 1)
                call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, 2)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur, 1)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, 1)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, 2)

                call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
                call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm',barotropicCoriolisTerm)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0', CGvec_r0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_p0', CGvec_p0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_v0', CGvec_v0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_h' , CGvec_h )
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_u0', CGvec_u0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_w0', CGvec_w0)

                call mpas_pool_get_array(diagnosticsPool, 'CGvec_z0', CGvec_z0)

                nCells = nCellsPtr
                nEdges = nEdgesPtr

!               nCells = nCellsArray(cellHaloComputeCounter)
!               nEdges = nEdgesArray(edgeHaloComputeCounter)

                nCells = nCellsArray( 1 )
                nEdges = nEdgesArray( 2 )

                CGcst_ruold = 0.0_RKIND





                do iCell = 1, nCells

                  sshTendb1 = 0.0_RKIND
                  sshTendb2 = 0.0_RKIND
                  sshTendAx = 0.0_RKIND

                  do i = 1, nEdgesOnCell(iCell)
                    iEdge = edgesOnCell(i, iCell)

                    cell1 = cellsOnEdge(1, iEdge)
                    cell2 = cellsOnEdge(2, iEdge)

                    ! Interpolation sshEdge
                    sshEdgeCur = 0.5_RKIND * (        sshCur(cell1) +         sshCur(cell2))
                    sshEdgeMid = 0.5_RKIND * (       CGvec_h(cell1) +        CGvec_h(cell2))
                    sshEdgeLag = 0.5_RKIND * (sshSubcycleCur(cell1) + sshSubcycleCur(cell2))

                    ! method 1, matches method 0 without pbcs, works with pbcs.
                    thicknessSumCur = sshEdgeCur + min(bottomDepth(cell1), bottomDepth(cell2))
                    thicknessSumMid = sshEdgeMid + min(bottomDepth(cell1), bottomDepth(cell2))
                    thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1), bottomDepth(cell2))

                    ! nabla (ssh^0)
                    sshDiffCur = (        sshCur(cell2)-        sshCur(cell1)) / dcEdge(iEdge)
                    sshDiffMid = (       CGvec_h(cell2)-       CGvec_h(cell1)) / dcEdge(iEdge)
                    sshDiffLag = (sshSubcycleCur(cell2)-sshSubcycleCur(cell1)) / dcEdge(iEdge)

                    !--------------------------------------------------------------!
                       fluxb1 = thicknessSumMid * normalBarotropicVelocityCur(iEdge)
                       fluxb2 = thicknessSumLag * (0.5_RKIND*gravity*sshDiffCur + (-barotropicCoriolisTerm(iEdge)-barotropicForcing(iEdge)) )
                       fluxAx = thicknessSumLag * sshDiffLag
       
                    sshTendb1 = sshTendb1 + edgeSignOnCell(i, iCell) * fluxb1 * dvEdge(iEdge)
                    sshTendb2 = sshTendb2 + edgeSignOnCell(i, iCell) * fluxb2 * dvEdge(iEdge)
                    sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) * fluxAx * dvEdge(iEdge)
                    !--------------------------------------------------------------!

                  end do ! i
                    
                    sshTendb1 = (4.0_RKIND/(gravity*dt)) * sshTendb1
                    sshTendb2 = (2.0_RKIND/(gravity   )) * sshTendb2

                    sshCurArea = (4.0_RKIND/(gravity*dt**2.0)) *         sshCur(iCell) * areaCell(iCell)
                    sshLagArea = (4.0_RKIND/(gravity*dt**2.0)) * sshSubcycleCur(iCell) * areaCell(iCell)

                    CGvec_r0(iCell) = (-sshCurArea - sshTendb1 + sshTendb2)   &
                                     -(-sshLagArea - sshTendAx) 

                  ! Initialize v0 and p0
                  CGvec_p0(iCell) = CGvec_r0(iCell)
                  CGvec_u0(iCell) = CGvec_r0(iCell) / CGvec_z0(iCell)

                  CGcst_ruold = CGcst_ruold + CGvec_r0(iCell) * CGvec_u0(iCell)

                end do ! iCell

                block => block % next
              end do  ! block

              ! Global sum across CPUs
              call mpas_dmpar_sum_real(dminfo,CGcst_ruold,CGcst_ruold_global)

              call mpas_timer_start("si halo iter r0")
              call mpas_dmpar_exch_group_create(domain, subcycleGroupName)
              call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'CGvec_u0', 1)
              call mpas_threading_barrier()
              call mpas_dmpar_exch_group_full_halo_exch(domain, subcycleGroupName)
              call mpas_dmpar_exch_group_destroy(domain, subcycleGroupName)
              call mpas_timer_stop("si halo iter r0")



              !--- A * u0 = A * r0---------------------------------------------------------------------------!

              block => domain % blocklist
              do while (associated(block))
                call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
                call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
                call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
                call mpas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)
                call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
                call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

                call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
                call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
                call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                call mpas_pool_get_subpool(block % structs, 'state', statePool)
                call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
                call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

                call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
                call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
                call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
                call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
                call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
                call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
                call mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell', refBottomDepthTopOfCell)
                call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
                call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
                call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

                call mpas_pool_get_array(statePool, 'ssh', sshCur,1)
                call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, 1)
                call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, 2)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur,1)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, 1)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, 2)

                call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
                call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm',barotropicCoriolisTerm)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0', CGvec_r0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_p0', CGvec_p0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_v0', CGvec_v0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_u0', CGvec_u0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_w0', CGvec_w0)

                nCells = nCellsPtr
                nEdges = nEdgesPtr

                nCells = nCellsArray( 1 )
                nEdges = nEdgesArray( 2 )

                CGcst_wuold = 0.0_RKIND

                do iCell = 1, nCells

                  sshTendAx = 0.0_RKIND

                  do i = 1, nEdgesOnCell(iCell)
                    iEdge = edgesOnCell(i, iCell)

                    cell1 = cellsOnEdge(1, iEdge)
                    cell2 = cellsOnEdge(2, iEdge)

                    ! Interpolation sshEdge
                    sshEdgeLag = 0.5_RKIND * (sshSubcycleCur(cell1) + sshSubcycleCur(cell2))

                    ! method 1, matches method 0 without pbcs, works with pbcs.
                    thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1), bottomDepth(cell2))

                    ! nabla (ssh^0)
                    sshDiffLag = (CGvec_u0(cell2) - CGvec_u0(cell1)) / dcEdge(iEdge)

                    !--------------------------------------------------------------!
                       fluxAx = thicknessSumLag * sshDiffLag
       
                    sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) * fluxAx * dvEdge(iEdge)
                    !--------------------------------------------------------------!

                  end do ! i
                    
                    sshLagArea = (4.0_RKIND/(gravity*dt**2.0)) * CGvec_u0(iCell) * areaCell(iCell)

                    CGvec_w0(iCell) = -sshLagArea - sshTendAx

                    CGcst_wuold = CGcst_wuold + CGvec_w0(iCell) * CGvec_u0(iCell)

                end do ! iCell

              ! Global sum across CPUs
              call mpas_dmpar_sum_real(dminfo,CGcst_wuold,CGcst_wuold_global)

                CGcst_alpha0 = CGcst_ruold_global / CGcst_wuold_global

                CGcst_beta0  = 0.0_RKIND

                block => block % next
              end do  ! block

              call mpas_timer_stop("si btr first r0")



         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! Stage 2.5 : Main iteration
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              call mpas_timer_start("si btr iteration")

              iter = 0
              resid = (crit_out+100.0)**2.0

              !**************************************************************!
              do while ( dsqrt(resid) > crit_main ) 
              !**************************************************************!

              block => domain % blocklist
              do while (associated(block))

                call mpas_pool_get_dimension(block % dimensions, 'nCells'     , nCellsPtr  )
                call mpas_pool_get_dimension(block % dimensions, 'nEdges'     , nEdgesPtr  )
                call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
                call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

                call mpas_pool_get_subpool(block % structs, 'tend'       , tendPool       )
                call mpas_pool_get_subpool(block % structs, 'mesh'       , meshPool       )
                call mpas_pool_get_subpool(block % structs, 'state'      , statePool      )
                call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

                call mpas_pool_get_subpool(statePool, 'tracers'    , tracersPool    )

                call mpas_pool_get_subpool(tendPool , 'tracersTend', tracersTendPool)

                call mpas_pool_get_array(meshPool, 'nEdgesOnCell',            nEdgesOnCell           )
                call mpas_pool_get_array(meshPool, 'edgesOnCell',             edgesOnCell            )
                call mpas_pool_get_array(meshPool, 'cellsOnEdge',             cellsOnEdge            )
                call mpas_pool_get_array(meshPool, 'dcEdge',                  dcEdge                 )
                call mpas_pool_get_array(meshPool, 'bottomDepth',             bottomDepth            )
                call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop',         maxLevelEdgeTop        )
                call mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell', refBottomDepthTopOfCell)
                call mpas_pool_get_array(meshPool, 'edgeSignOnCell',          edgeSignOnCell         )
                call mpas_pool_get_array(meshPool, 'dvEdge',                  dvEdge                 )
                call mpas_pool_get_array(meshPool, 'areaCell',                areaCell               )
                call mpas_pool_get_array(meshPool, 'nEdgesOnEdge',            nEdgesOnEdge           )
                call mpas_pool_get_array(meshPool, 'edgesOnEdge',             edgesOnEdge            )
                call mpas_pool_get_array(meshPool, 'weightsOnEdge',           weightsOnEdge          )
                call mpas_pool_get_array(meshPool, 'fEdge',                   fEdge                  )
                call mpas_pool_get_array(meshPool, 'edgeMask',                edgeMask               )

                call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
                call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)
                call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, 1)
                call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, 2)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur,1)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew,2)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, 1)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, 2)

                call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing',      barotropicForcing     )
                call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm', barotropicCoriolisTerm)

                call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0', CGvec_r0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_p0', CGvec_p0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_v0', CGvec_v0)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_r1', CGvec_r1)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_p1', CGvec_p1)
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_t' , CGvec_t )
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_h' , CGvec_h )

                call mpas_pool_get_array(diagnosticsPool, 'CGvec_s0' , CGvec_s0 )
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_u0' , CGvec_u0 )
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_w0' , CGvec_w0 )
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_s1' , CGvec_s1 )
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_u1' , CGvec_u1 )
                call mpas_pool_get_array(diagnosticsPool, 'CGvec_w1' , CGvec_w1 )

                call mpas_pool_get_array(diagnosticsPool, 'CGvec_z0', CGvec_z0)

                call mpas_pool_get_array(tendPool, 'ssh', sshTend)

                nCells = nCellsPtr
                nEdges = nEdgesPtr

                nCells = nCellsArray( 1 )
                nEdges = nEdgesArray( 2 )

                iter = iter + 1

! 1. p1 = u0 + beta0 * p0
                do iCell = 1, nCells
                  CGvec_p1(iCell) = CGvec_u0(iCell) + CGcst_beta0 * CGvec_p0(iCell)
                end do ! iCell

! 2. s1 = w0 + beta0 * s0
                do iCell = 1, nCells
                  CGvec_s1(iCell) = CGvec_w0(iCell) + CGcst_beta0 * CGvec_s0(iCell)
                end do ! iCell

! 3. x1 = x0 + alpha0 * p1
                do iCell = 1,nCells
 		  sshSubcycleNew(iCell) = sshSubcycleCur(iCell) + CGcst_alpha0 * CGvec_p1(iCell)
              	end do ! iCell

! 4. r1 = r0 - alpha * s1
                do iCell = 1,nCells
 		  CGvec_r1(iCell) = CGvec_r0(iCell) - CGcst_alpha0 * CGvec_s1(iCell) 
              	end do ! iCell

! 5. u1 = r1
                do iCell = 1,nCells
 		  CGvec_u1(iCell) = CGvec_r1(iCell) / CGvec_z0(iCell) 
              	end do ! iCell

               call mpas_timer_start("si halo iter")
               call mpas_dmpar_exch_group_create(domain, subcycleGroupName)
               call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'CGvec_u1', 1)
               call mpas_threading_barrier()
               call mpas_dmpar_exch_group_full_halo_exch(domain, subcycleGroupName)
               call mpas_dmpar_exch_group_destroy(domain, subcycleGroupName)
               call mpas_timer_stop("si halo iter")


!---------------------------------------------------------------------!
! 6. w1 = A * u1


                do iCell = 1,nCells
                   sshTend(iCell) = 0.0_RKIND
                end do

                do iEdge = 1,nEdges
                  cell1 = cellsOnEdge(1, iEdge)
                  cell2 = cellsOnEdge(2, iEdge)

                    ! Interpolation sshEdge
                    sshEdgeLag = 0.5_RKIND * (sshSubcycleCur(cell1) + sshSubcycleCur(cell2))

                    thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1), bottomDepth(cell2))

                    ! nabla (ssh^0)
                    sshDiffNew = (CGvec_u1(cell2)-CGvec_u1(cell1)) / dcEdge(iEdge)

                        fluxAx = thicknessSumLag * sshDiffNew
       
                    sshTend(cell1) = sshTend(cell1) - fluxAx * dvEdge(iEdge)
                    sshTend(cell2) = sshTend(cell2) + fluxAx * dvEdge(iEdge)
                end do          

                CGcst_runew = 0.0_RKIND
                CGcst_wunew = 0.0_RKIND
                CGcst_rrnew = 0.0_RKIND

                do iCell = 1,nCells
                     sshLagArea = (4.0_RKIND/(gravity*dt**2.0)) * CGvec_u1(iCell) * areaCell(iCell)

                  CGvec_w1(iCell) = -sshLagArea - sshTend(iCell)

                  CGcst_runew = CGcst_runew + CGvec_r1(iCell) * CGvec_u1(iCell)
                  CGcst_wunew = CGcst_wunew + CGvec_w1(iCell) * CGvec_u1(iCell)
                  CGcst_rrnew = CGcst_rrnew + CGvec_r1(iCell) * CGvec_r1(iCell)
                end do

!---------------------------------------------------------------------!
       
                  CGcst_allreduce(1) = CGcst_runew
                  CGcst_allreduce(2) = CGcst_wunew
                  CGcst_allreduce(3) = CGcst_rrnew

                ! Global sum across CPUs
               call mpas_dmpar_sum_real_array(dminfo, 3, CGcst_allreduce, CGcst_allreduce_global)

                  CGcst_runew_global = CGcst_allreduce_global(1)
                  CGcst_wunew_global = CGcst_allreduce_global(2)
                               resid = CGcst_allreduce_global(3)

! 6. beta1 = runew / ruold
               CGcst_beta1 = CGcst_runew_global / CGcst_ruold_global

! 7. alpha1 = 1 / (wunew / runew - beta1 / alpha0 ) 
               CGcst_alpha1 = 1.0_RKIND / ( CGcst_wunew_global/CGcst_runew_global - CGcst_beta1/CGcst_alpha0 )

! 14. Stepping for next iteration
                do iCell = 1,nCells
                  CGvec_r0(iCell) = CGvec_r1(iCell) 
                  CGvec_p0(iCell) = CGvec_p1(iCell)
                  CGvec_s0(iCell) = CGvec_s1(iCell)
                  CGvec_u0(iCell) = CGvec_u1(iCell)
                  CGvec_w0(iCell) = CGvec_w1(iCell)

                  temp1 = sshCur(iCell)
                  temp2 = sshSubcycleNew(iCell)

                  sshSubcycleCur(iCell) = sshSubcycleNew(iCell)
                   
                  ! sshSubcycleNew => sshSubcycleNew^* for u update
                    CGvec_h(iCell) = 0.5_RKIND * (temp1+temp2)
              	end do ! iCell
		
  		  CGcst_rrold_global = CGcst_rrnew_global
  		  CGcst_ruold_global = CGcst_runew_global
  		  CGcst_wuold_global = CGcst_wunew_global
  		  CGcst_beta0        = CGcst_beta1
  		  CGcst_alpha0       = CGcst_alpha1

                block => block % next
              end do  ! block


               !   boundary update on SSHnew
               call mpas_timer_start("si halo iter")
               call mpas_dmpar_exch_group_create(domain, subcycleGroupName)
               call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'sshSubcycle', 1)
               call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'CGvec_h', 1)
               call mpas_threading_barrier()
               call mpas_dmpar_exch_group_full_halo_exch(domain, subcycleGroupName)
               call mpas_dmpar_exch_group_destroy(domain, subcycleGroupName)
               call mpas_timer_stop("si halo iter")


              ! Update normalBarotropicVelocitySubcycle
              
                  block => domain % blocklist
                  do while (associated(block))
                     call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
                     call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

                     call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                     call mpas_pool_get_subpool(block % structs, 'state', statePool)
                     call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
                     call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

                     call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
                     call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
                     call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
                     call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
                     call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
                     call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
                     call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

                     call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, 1)
                     call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, 2)
                     call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur,1)
                     call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew,2)
!                    call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, 2)
                     call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)

                     call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
                     call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm',barotropicCoriolisTerm)
                     call mpas_pool_get_array(diagnosticsPool, 'CGvec_h', CGvec_h)

                     nEdges = nEdgesPtr
                     nEdges = nEdgesArray( edgeHaloComputeCounter )
!                    nEdges = nEdgesArray( 1)

                     !$omp do schedule(runtime) private(cell1, cell2, CoriolisTerm, i, eoe)
                     do iEdge = 1, nEdges

                        temp_mask = edgeMask(1, iEdge)

                          cell1 = cellsOnEdge(1,iEdge)
                          cell2 = cellsOnEdge(2,iEdge)

                          ! Compute the barotropic Coriolis term, -f*uPerp
                          CoriolisTerm = 0.0_RKIND
                                 temp1 = 0.0_RKIND
                          do i = 1, nEdgesOnEdge(iEdge)
                             eoe = edgesOnEdge(i,iEdge)

                             CoriolisTerm = CoriolisTerm + weightsOnEdge(i,iEdge)  &
                                          * normalBarotropicVelocityCur(eoe) * fEdge(eoe)
                                    temp1 =        temp1 + weightsOnEdge(i,iEdge)  &
                                          * normalBarotropicVelocitySubcycleCur(eoe) * fEdge(eoe)
                          end do
                             barotropicCoriolisTerm(iEdge) = 0.5_RKIND * (CoriolisTerm+temp1)

                          normalBarotropicVelocitySubcycleNew(iEdge) &
                            = temp_mask &
                            * (normalBarotropicVelocityCur(iEdge) &
                            - dt * (-barotropicCoriolisTerm(iEdge) + gravity &
                            * ( CGvec_h(cell2) - CGvec_h(cell1)) / dcEdge(iEdge) - barotropicForcing(iEdge)))
                     end do
                     !$omp end do

                     do iEdge = 1, nEdges
                          ! Cur = New for next iteration
                          normalBarotropicVelocitySubcycleCur(iEdge) = normalBarotropicVelocitySubcycleNew(iEdge) 
                     end do

                     block => block % next
                  end do  ! block


                  !   boundary update on SSHnew
!                 call mpas_timer_start("si halo iter")
!                 call mpas_dmpar_exch_group_create(domain, subcycleGroupName)
!                 call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'normalBarotropicVelocitySubcycle', 1)
!                 call mpas_threading_barrier()
!                 call mpas_dmpar_exch_group_full_halo_exch(domain, subcycleGroupName)
!                 call mpas_dmpar_exch_group_destroy(domain, subcycleGroupName)
!                 call mpas_timer_stop("si halo iter")


                if ( iter > total_num_cells ) then
                  print*,'******************************************************'
                  print*,'******************************************************'
                  print*,'Iteration number exceeds Max. #iteration: PROGRAM STOP'
                  print*,'******************************************************'
                  print*,'******************************************************'
                  stop
                endif

              !**************************************************************!
              end do ! do while iter
              !**************************************************************!

              call mpas_timer_stop("si btr iteration")



                  block => domain % blocklist
                  do while (associated(block))
                     call mpas_pool_get_dimension(block % dimensions, 'nCells'     , nCellsPtr  )
                     call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
                     call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
                     call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

                     call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                     call mpas_pool_get_subpool(block % structs, 'state', statePool)
                     call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
                     call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
                     call mpas_pool_get_subpool(block % structs, 'tend'       , tendPool       )

                     call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
                     call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
                     call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
                     call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
                     call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
                     call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
                     call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

                     call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, 1)
                     call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew,2)
                     call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, 1)
                     call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
                     call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)

                     nCells = nCellsPtr
                     nEdges = nEdgesPtr

!                    call mpas_log_write( 'iter $i', intArgs=(/ iter /) )
!                    print*,'s_sqrt_main ',dsqrt(resid),iter,split_explicit_step

                     do iCell = 1, nCells
                       sshNew(iCell) = sshSubcycleCur(iCell) ! Already stepped, so we use 'Cur'
                     end do ! iCell


                      !--- Update u_lag ------------------------------------------------------!

                     !$omp do schedule(runtime) private(cell1, cell2, CoriolisTerm, i, eoe)
                     do iEdge = 1, nEdges
                          normalBarotropicVelocityNew(iEdge) = normalBarotropicVelocitySubcycleCur(iEdge)
                     end do
                     !$omp end do

               block => block % next
            end do  ! block

                   call mpas_timer_start("se halo btr vel corr")
                   call mpas_dmpar_field_halo_exch(domain, 'normalBarotropicVelocity', timeLevel=2)
                   call mpas_dmpar_field_halo_exch(domain, 'ssh', timeLevel=2)
                   call mpas_timer_stop("se halo btr vel corr")


        !--- Update barotropicThicknessFlux ---------------------------------------------------!

              block => domain % blocklist
              do while (associated(block))
                call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
                call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
                call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
                call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

                call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
                call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
                call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                call mpas_pool_get_subpool(block % structs, 'state', statePool)
                call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
                call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

                call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
                call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
                call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
                call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
                call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
                call mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell', refBottomDepthTopOfCell)
                call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
                call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
                call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

                call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
                call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur,1)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew,2)

                call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', barotropicThicknessFlux)

                nCells = nCellsPtr
                nEdges = nEdgesPtr

                nCells = nCellsArray( cellHaloComputeCounter )
                nEdges = nEdgesArray( edgeHaloComputeCounter )

                ! config_btr_gam1_velWt1 sets the forward weighting of velocity in the SSH computation
                ! config_btr_gam1_velWt1=  1     flux = normalBarotropicVelocityNew*H
                ! config_btr_gam1_velWt1=0.5     flux = 1/2*(normalBarotropicVelocityNew+normalBarotropicVelocityOld)*H
                ! config_btr_gam1_velWt1=  0     flux = normalBarotropicVelocityOld*H

                !! asarje: changed to avoid redundant computations when config_btr_solve_SSH2 is true

                if (config_btr_solve_SSH2) then

                  ! If config_btr_solve_SSH2=.true.,
                  ! then do NOT accumulate barotropicThicknessFlux in this SSH predictor
                  ! section, because it will be accumulated in the SSH corrector section.
                  barotropicThicknessFlux_coeff = 0.0_RKIND

                  ! othing else to do

                else

                  ! otherwise, DO accumulate barotropicThicknessFlux in this SSH predictor section
                  barotropicThicknessFlux_coeff = 1.0_RKIND

                  !$omp do schedule(runtime)
                  do iEdge = 1, nEdges
                     cell1 = cellsOnEdge(1,iEdge)
                     cell2 = cellsOnEdge(2,iEdge)

!                    sshEdge = 0.5_RKIND * ( 0.5_RKIND * (sshCur(cell1) + sshCur(cell2)) &
!                                          + 0.5_RKIND * (sshNew(cell1) + sshNew(cell2)))

                     sshEdge = 0.25_RKIND * ( 0.5_RKIND * (sshCur(cell1) + sshCur(cell2)) &
                                            + 0.5_RKIND * (sshNew(cell1) + sshNew(cell2)))
                     
                     ! method 1, matches method 0 without pbcs, works with pbcs.
                     thicknessSum = sshEdge + min(bottomDepth(cell1), bottomDepth(cell2))

!                    barotropicThicknessFlux(iEdge) =0.5_RKIND*( normalBarotropicVelocityCur(iEdge) &
!                                                               +normalBarotropicVelocityNew(iEdge) )* thicknessSum

                     barotropicThicknessFlux(iEdge) =0.25_RKIND*( normalBarotropicVelocityCur(iEdge) &
                                                                +normalBarotropicVelocityNew(iEdge) )* thicknessSum

                  end do
                  !$omp end do

                endif

                !********************************************************************
                !* Here is ssh_(1/2) and u_(1/2) ************************************
                !********************************************************************
                do iCell = 1,nCells
                  sshCur(iCell) = 0.5_RKIND * (sshCur(iCell) + sshNew(iCell))
                end do
                do iEdge = 1,nEdges
                  normalBarotropicVelocityCur(iCell) = 0.5_RKIND * (normalBarotropicVelocityCur(iCell) + &
                                                                    normalBarotropicVelocityNew(iCell))
                end do

                block => block % next
              end do  ! block

                   call mpas_timer_start("se halo btr vel corr")
                   call mpas_dmpar_field_halo_exch(domain, 'normalBarotropicVelocity', timeLevel=1)
                   call mpas_dmpar_field_halo_exch(domain, 'ssh', timeLevel=1)

            ! boundary update on F
            call mpas_timer_start("se halo F and btr vel")
            call mpas_dmpar_exch_group_create(domain, finalBtrGroupName)
            call mpas_dmpar_exch_group_add_field(domain, finalBtrGroupName, 'barotropicThicknessFlux')
            call mpas_threading_barrier()
            call mpas_dmpar_exch_group_full_halo_exch(domain, finalBtrGroupName)
            call mpas_dmpar_exch_group_destroy(domain, finalBtrGroupName)
            call mpas_timer_stop("se halo F and btr vel")
            call mpas_threading_barrier()



	! u correct and transport velocity

            call mpas_timer_start('btr si u cor and trans')
            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
               call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

               call mpas_pool_get_subpool(block % structs, 'state', statePool)
               call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
               call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

               call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityNew, 2)
               call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew, 2)

               call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
               call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)
               call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
               call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', barotropicThicknessFlux)

               call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
               call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

               nEdges = nEdgesPtr

               nEdges = nEdgesArray( config_num_halos )

               allocate(uTemp(nVertLevels))

               ! Correction velocity    normalVelocityCorrection = (Flux - Sum(h u*))/H
               ! or, for the full latex version:
               !{\bf u}^{corr} = \left( {\overline {\bf F}}
               !  - \sum_{k=1}^{N^{edge}} h_{k,*}^{edge}  {\bf u}_k^{avg} \right)
               ! \left/ \sum_{k=1}^{N^{edge}} h_{k,*}^{edge}   \right.

               if (config_vel_correction) then
                  useVelocityCorrection = 1
               else
                  useVelocityCorrection = 0
               endif

               !$omp do schedule(runtime) private(k, normalThicknessFluxSum, thicknessSum, normalVelocityCorrection)
               do iEdge = 1, nEdges

                  ! velocity for normalVelocityCorrectionection is normalBarotropicVelocity + normalBaroclinicVelocity + uBolus
                  uTemp(:) = normalBarotropicVelocityNew(iEdge) + normalBaroclinicVelocityNew(:,iEdge) &
                           + normalGMBolusVelocity(:,iEdge)

                  ! thicknessSum is initialized outside the loop because on land boundaries
                  ! maxLevelEdgeTop=0, but I want to initialize thicknessSum with a
                  ! nonzero value to avoid a NaN.
                  normalThicknessFluxSum = layerThicknessEdge(1,iEdge) * uTemp(1)
                  thicknessSum  = layerThicknessEdge(1,iEdge)

                  do k = 2, maxLevelEdgeTop(iEdge)
                     normalThicknessFluxSum = normalThicknessFluxSum + layerThicknessEdge(k,iEdge) * uTemp(k)
                     thicknessSum  =  thicknessSum + layerThicknessEdge(k,iEdge)
                  enddo

                  normalVelocityCorrection = useVelocityCorrection * (( barotropicThicknessFlux(iEdge) - normalThicknessFluxSum) &
                                           / thicknessSum)

                  do k = 1, nVertLevels

                     ! normalTransportVelocity = normalBarotropicVelocity + normalBaroclinicVelocity + normalGMBolusVelocity
                     !                         + normalVelocityCorrection
                     ! This is u used in advective terms for layerThickness and tracers
                     ! in tendency calls in stage 3.
!mrp note: in QC version, there is an if (config_use_standardGM) on adding normalGMBolusVelocity
! I think it is not needed because normalGMBolusVelocity=0 when GM not on.
                     normalTransportVelocity(k,iEdge) &
                           = edgeMask(k,iEdge) &
                           *( normalBarotropicVelocityNew(iEdge) + normalBaroclinicVelocityNew(k,iEdge) &
                           + normalGMBolusVelocity(k,iEdge) + normalVelocityCorrection )
                  enddo

               end do ! iEdge
               !$omp end do

               deallocate(uTemp)

               block => block % next
            end do  ! block
            call mpas_timer_stop('btr si u cor and trans')


         call mpas_timer_stop("stage 2 timer")

!****************************************************************************************         
!****************************************************************************************         
!****************************************************************************************         
!****************************************************************************************         
!****************************************************************************************         

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         !  Stage 3: Tracer, density, pressure, vertical velocity prediction
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! only compute tendencies for active tracers on last large iteration
         if (split_explicit_step < config_n_ts_iter) then
            activeTracersOnly = .true.
         else
            activeTracersOnly = .false.
         endif

         ! Thickness tendency computations and thickness halo updates are completed before tracer
         ! tendency computations to allow monotonic advection.
         call mpas_timer_start('se thick tend')
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)

            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
            call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
            call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)

            call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
            call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
            call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)

            ! compute vertAleTransportTop.  Use normalTransportVelocity for advection of layerThickness and tracers.
            ! Use time level 1 values of layerThickness and layerThicknessEdge because
            ! layerThickness has not yet been computed for time level 2.
            call mpas_timer_start('thick vert trans vel top')
            if (associated(highFreqThicknessNew)) then
               call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
                 layerThicknessCur, layerThicknessEdge, normalTransportVelocity, &
                 sshCur, dt, vertAleTransportTop, err, highFreqThicknessNew)
            else
               call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
                 layerThicknessCur, layerThicknessEdge, normalTransportVelocity, &
                 sshCur, dt, vertAleTransportTop, err)
            endif
            call mpas_timer_stop('thick vert trans vel top')

            call ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool)

            block => block % next
         end do
         call mpas_timer_stop('se thick tend')

         call mpas_threading_barrier()

         ! update halo for thickness tendencies
         call mpas_timer_start("se halo thickness")

         call mpas_dmpar_field_halo_exch(domain, 'tendLayerThickness')

         call mpas_timer_stop("se halo thickness")

         call mpas_threading_barrier()

         call mpas_timer_start('se tracer tend', .false.)
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            call mpas_pool_get_subpool(block % structs, 'shortwave', swForcingPool)
            call ocn_tend_tracer(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, swForcingPool, scratchPool, &
                    dt, activeTracersOnly, 2)

            block => block % next
         end do
         call mpas_timer_stop('se tracer tend')

         call mpas_threading_barrier()

         ! update halo for tracer tendencies
         call mpas_timer_start("se halo tracers")
         call mpas_pool_get_subpool(domain % blocklist % structs, 'tend', tendPool)
         call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

         call mpas_pool_begin_iteration(tracersTendPool)
         do while ( mpas_pool_get_next_member(tracersTendPool, groupItr) )
            if ( groupItr % memberType == MPAS_POOL_FIELD ) then
               ! Only compute tendencies for active tracers if activeTracersOnly flag is true.
               if ( .not.activeTracersOnly .or. trim(groupItr % memberName)=='activeTracersTend') then
                  call mpas_dmpar_field_halo_exch(domain, groupItr % memberName)
               end if
            end if
         end do
         call mpas_timer_stop("se halo tracers")

         call mpas_threading_barrier()

         call mpas_timer_start('se loop fini')
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
            call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
            call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
            call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
            call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

            call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
            call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)
            call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

            call mpas_pool_get_array(tracersPool, 'activeTracers', tracersGroupCur, 1)
            call mpas_pool_get_array(tracersPool, 'activeTracers', tracersGroupNew, 2)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)
            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
            call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
            call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)
            call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceCur, 1)
            call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceNew, 2)
            call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur, 1)
            call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew, 2)
            call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityCur, 1)
            call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityNew, 2)

            call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)
            call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
            call mpas_pool_get_array(tendPool, 'highFreqThickness', highFreqThicknessTend)
            call mpas_pool_get_array(tendPool, 'lowFreqDivergence', lowFreqDivergenceTend)

            call mpas_pool_get_array(tracersTendPool, 'activeTracersTend', activeTracersTend)

            nCells = nCellsPtr
            nEdges = nEdgesPtr

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !  If iterating, reset variables for next iteration
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (split_explicit_step < config_n_ts_iter) then

               ! Get indices for dynamic tracers (Includes T&S).
               call mpas_pool_get_dimension(tracersPool, 'activeGRP_start', startIndex)
               call mpas_pool_get_dimension(tracersPool, 'activeGRP_end', endIndex)

               ! Only need T & S for earlier iterations,
               ! then all the tracers needed the last time through.

               !$omp do schedule(runtime) private(k, temp_h, temp, i)
               do iCell = 1, nCells
                  ! sshNew is a pointer, defined above.
                  do k = 1, maxLevelCell(iCell)

                     ! this is h_{n+1}
                     temp_h = layerThicknessCur(k,iCell) + dt * layerThicknessTend(k,iCell)

                     ! this is h_{n+1/2}
                     layerThicknessNew(k,iCell) = 0.5*( layerThicknessCur(k,iCell) + temp_h)

                     do i = startIndex, endIndex
                        ! This is Phi at n+1
                        temp = ( tracersGroupCur(i,k,iCell) * layerThicknessCur(k,iCell) + dt * activeTracersTend(i,k,iCell)) &
                             / temp_h

                        ! This is Phi at n+1/2
                        tracersGroupNew(i,k,iCell) = 0.5_RKIND * ( tracersGroupCur(i,k,iCell) + temp )
                     end do
                  end do
               end do ! iCell
               !$omp end do

               if (config_use_freq_filtered_thickness) then
                  !$omp do schedule(runtime) private(k, temp)
                  do iCell = 1, nCells
                     do k = 1, maxLevelCell(iCell)

                        ! h^{hf}_{n+1} was computed in Stage 1

                        ! this is h^{hf}_{n+1/2}
                        highFreqThicknessnew(k,iCell) = 0.5_RKIND * (highFreqThicknessCur(k,iCell) + highFreqThicknessNew(k,iCell))

                        ! this is D^{lf}_{n+1}
                        temp = lowFreqDivergenceCur(k,iCell) &
                         + dt * lowFreqDivergenceTend(k,iCell)

                        ! this is D^{lf}_{n+1/2}
                        lowFreqDivergenceNew(k,iCell) = 0.5_RKIND * (lowFreqDivergenceCur(k,iCell) + temp)
                     end do
                  end do
                  !$omp end do
               end if

               !$omp do schedule(runtime) private(k)
               do iEdge = 1, nEdges

                  do k = 1, nVertLevels

                     ! u = normalBarotropicVelocity + normalBaroclinicVelocity
                     ! here normalBaroclinicVelocity is at time n+1/2
                     ! This is u used in next iteration or step
                     normalVelocityNew(k,iEdge) = edgeMask(k,iEdge) * ( normalBarotropicVelocityNew(iEdge) &
                                                + normalBaroclinicVelocityNew(k,iEdge) )

                  enddo

               end do ! iEdge
               !$omp end do

               ! Efficiency note: We really only need this to compute layerThicknessEdge, density, pressure, and SSH
               ! in this diagnostics solve.
               call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 2)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !  If large iteration complete, compute all variables at time n+1
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            elseif (split_explicit_step == config_n_ts_iter) then

               !$omp do schedule(runtime) private(k)
               do iCell = 1, nCells
                  do k = 1, maxLevelCell(iCell)
                     ! this is h_{n+1}
                     layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell) + dt * layerThicknessTend(k,iCell)
                  end do
               end do
               !$omp end do

               if (config_compute_active_tracer_budgets) then
                  call mpas_pool_get_array(diagnosticsPool,'activeTracerHorizontalAdvectionTendency', &
                          activeTracerHorizontalAdvectionTendency)
                  call mpas_pool_get_array(diagnosticspool,'activeTracerVerticalAdvectionTendency', &
                          activeTracerVerticalAdvectionTendency)
                  call mpas_pool_get_array(diagnosticsPool,'activeTracerSurfaceFluxTendency',activeTracerSurfaceFluxTendency)
                  call mpas_pool_get_array(diagnosticsPool,'temperatureShortWaveTendency',temperatureShortWaveTendency)
                  call mpas_pool_get_array(diagnosticsPool,'activeTracerNonLocalTendency',activeTracerNonLocalTendency)
               !$omp do schedule(runtime) private(k)
                  do iCell = 1, nCells
                     do k= 1, maxLevelCell(iCell)
                        activeTracerHorizontalAdvectionTendency(:,k,iCell) = &
                           activeTracerHorizontalAdvectionTendency(:,k,iCell) / &
                           layerThicknessNew(k,iCell)
        
                        activeTracerVerticalAdvectionTendency(:,k,iCell) = &
                           activeTracerVerticalAdvectionTendency(:,k,iCell) / &
                           layerThicknessNew(k,iCell)

                        activeTracerSurfaceFluxTendency(:,k,iCell) = &
                           activeTracerSurfaceFluxTendency(:,k,iCell) / &
                           layerThicknessNew(k,iCell)

                        temperatureShortWaveTendency(k,iCell) = &
                           temperatureShortWaveTendency(k,iCell) / &
                           layerThicknessNew(k,iCell)

                        activeTracerNonLocalTendency(:,k,iCell) = &
                           activeTracerNonLocalTendency(:,k,iCell) / &
                           layerThicknessNew(k,iCell)
                     end do
                  end do 
                  !$omp end do
               endif

               call mpas_pool_begin_iteration(tracersPool)
               do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
                  if ( groupItr % memberType == MPAS_POOL_FIELD ) then
                     configName = 'config_use_' // trim(groupItr % memberName)
                     call mpas_pool_get_config(domain % configs, configName, config_use_tracerGroup)

                     if ( config_use_tracerGroup ) then
                        call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupCur, 1)
                        call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupNew, 2)

                        modifiedGroupName = trim(groupItr % memberName) // 'Tend'
                        call mpas_pool_get_array(tracersTendPool, modifiedGroupName, tracersGroupTend)

                        !$omp do schedule(runtime) private(k)
                        do iCell = 1, nCells
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupNew(:,k,iCell) = (tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) + dt &
                                                         * tracersGroupTend(:,k,iCell) ) / layerThicknessNew(k,iCell)
                           end do
                        end do
                        !$omp end do

                        ! limit salinity in separate loop
                        if ( trim(groupItr % memberName) == 'activeTracers' ) then
                           !$omp do schedule(runtime) private(k)
                           do iCell = 1, nCells
                              do k = 1, maxLevelCell(iCell)
                                 tracersGroupNew(indexSalinity,k,iCell) = max(0.001_RKIND, tracersGroupNew(indexSalinity,k,iCell))
                              end do
                           end do
                           !$omp end do
                        end if

                     end if
                  end if
               end do

               if (config_use_freq_filtered_thickness) then
                  !$omp do schedule(runtime) private(k)
                  do iCell = 1, nCells
                     do k = 1, maxLevelCell(iCell)

                        ! h^{hf}_{n+1} was computed in Stage 1

                        ! this is D^{lf}_{n+1}
                        lowFreqDivergenceNew(k,iCell) = lowFreqDivergenceCur(k,iCell) + dt * lowFreqDivergenceTend(k,iCell)
                     end do
                  end do
                  !$omp end do
               end if

               ! Recompute final u to go on to next step.
               ! u_{n+1} = normalBarotropicVelocity_{n+1} + normalBaroclinicVelocity_{n+1}
               ! Right now normalBaroclinicVelocityNew is at time n+1/2, so back compute to get normalBaroclinicVelocity
               !   at time n+1 using normalBaroclinicVelocity_{n+1/2} = 1/2*(normalBaroclinicVelocity_n + u_Bcl_{n+1})
               ! so the following lines are
               ! u_{n+1} = normalBarotropicVelocity_{n+1} + 2*normalBaroclinicVelocity_{n+1/2} - normalBaroclinicVelocity_n
               ! note that normalBaroclinicVelocity is recomputed at the beginning of the next timestep due to Imp Vert mixing,
               ! so normalBaroclinicVelocity does not have to be recomputed here.

               !$omp do schedule(runtime) private(k)
               do iEdge = 1, nEdges
                  do k = 1, maxLevelEdgeTop(iEdge)
                     normalVelocityNew(k,iEdge) = normalBarotropicVelocityNew(iEdge) + 2 * normalBaroclinicVelocityNew(k,iEdge) &
                                                - normalBaroclinicVelocityCur(k,iEdge)
                  end do
               end do ! iEdges
               !$omp end do

            endif ! split_explicit_step

            block => block % next
         end do

         call mpas_timer_stop('se loop fini')
         call mpas_timer_stop('se loop')

      end do  ! split_explicit_step = 1, config_n_ts_iter
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END large iteration loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call mpas_timer_start("se implicit vert mix")

      block => domain % blocklist
      do while(associated(block))
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
        call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
        call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

        ! Call ocean diagnostic solve in preparation for vertical mixing.  Note
        ! it is called again after vertical mixing, because u and tracers change.
        ! For Richardson vertical mixing, only density, layerThicknessEdge, and kineticEnergyCell need to
        ! be computed.  For kpp, more variables may be needed.  Either way, this
        ! could be made more efficient by only computing what is needed for the
        ! implicit vmix routine that follows.
        call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 2)

        ! Compute normalGMBolusVelocity; it will be added to the baroclinic modes in Stage 2 above.
        if (config_use_standardGM) then
           call ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, scratchPool)
        end if
        call ocn_vmix_implicit(dt, meshPool, diagnosticsPool, statePool, forcingPool, scratchPool, err, 2)

        block => block % next
      end do

      call mpas_threading_barrier()

      ! Update halo on u and tracers, which were just updated for implicit vertical mixing.  If not done,
      ! this leads to lack of volume conservation.  It is required because halo updates in stage 3 are only
      ! conducted on tendencies, not on the velocity and tracer fields.  So this update is required to
      ! communicate the change due to implicit vertical mixing across the boundary.
      call mpas_timer_start('se vmix halos')
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)


      call mpas_timer_start('se vmix halos normalVelFld')
      call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=2)
      call mpas_timer_stop('se vmix halos normalVelFld')

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_dmpar_field_halo_exch(domain, groupItr % memberName, timeLevel=2)
         end if
      end do
      call mpas_timer_stop('se vmix halos')

      call mpas_timer_stop("se implicit vert mix")

      call mpas_threading_barrier()

      call mpas_timer_start('se fini')
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
         call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
         call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

         call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'velocityX', velocityX)
         call mpas_pool_get_array(diagnosticsPool, 'velocityY', velocityY)
         call mpas_pool_get_array(diagnosticsPool, 'velocityZ', velocityZ)
         call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velocityZonal)
         call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', velocityMeridional)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSH', gradSSH)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHX', gradSSHX)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHY', gradSSHY)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHZ', gradSSHZ)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHZonal', gradSSHZonal)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHMeridional', gradSSHMeridional)

         call mpas_pool_get_array(diagnosticsPool, 'surfaceVelocity', surfaceVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'SSHGradient', SSHGradient)

         call mpas_pool_get_dimension(diagnosticsPool, 'index_surfaceVelocityZonal', indexSurfaceVelocityZonal)
         call mpas_pool_get_dimension(diagnosticsPool, 'index_surfaceVelocityMeridional', indexSurfaceVelocityMeridional)
         call mpas_pool_get_dimension(diagnosticsPool, 'index_SSHGradientZonal', indexSSHGradientZonal)
         call mpas_pool_get_dimension(diagnosticsPool, 'index_SSHGradientMeridional', indexSSHGradientMeridional)

         nCells = nCellsPtr
         nEdges = nEdgesPtr

         if (config_prescribe_velocity) then
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               normalVelocityNew(:, iEdge) = normalVelocityCur(:, iEdge)
            end do
            !$omp end do
         end if

         if (config_prescribe_thickness) then
            !$omp do schedule(runtime)
            do iCell = 1, nCells
               layerThicknessNew(:, iCell) = layerThicknessCur(:, iCell)
            end do
            !$omp end do
         end if

         call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 2)

         ! Update the effective desnity in land ice if we're coupling to land ice
         call ocn_effective_density_in_land_ice_update(meshPool, forcingPool, statePool, scratchPool, err)

         ! Compute normalGMBolusVelocity; it will be added to normalVelocity in Stage 2 of the next cycle.
         if (config_use_standardGM) then
            call ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, scratchPool)
         end if

         call mpas_timer_start('se final mpas reconstruct', .false.)

         call mpas_reconstruct(meshPool, normalVelocityNew,  &
                          velocityX, velocityY, velocityZ,   &
                          velocityZonal, velocityMeridional, &
                          includeHalos = .true.)

         call mpas_reconstruct(meshPool, gradSSH,          &
                          gradSSHX, gradSSHY, gradSSHZ,    &
                          gradSSHZonal, gradSSHMeridional, &
                          includeHalos = .true.)

         call mpas_timer_stop('se final mpas reconstruct')

         !$omp do schedule(runtime)
         do iCell = 1, nCells
            surfaceVelocity(indexSurfaceVelocityZonal, iCell) = velocityZonal(1, iCell)
            surfaceVelocity(indexSurfaceVelocityMeridional, iCell) = velocityMeridional(1, iCell)

            SSHGradient(indexSSHGradientZonal, iCell) = gradSSHZonal(iCell)
            SSHGradient(indexSSHGradientMeridional, iCell) = gradSSHMeridional(iCell)
         end do
         !$omp end do

         call ocn_time_average_coupled_accumulate(diagnosticsPool, statePool, forcingPool, 2)
         call mpas_threading_barrier()

         if (config_use_standardGM) then
            call ocn_reconstruct_gm_vectors(diagnosticsPool, meshPool)
         end if

         block => block % next
      end do

      if (trim(config_land_ice_flux_mode) == 'coupled') then
         call mpas_timer_start("se effective density halo")
         call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
         call mpas_pool_get_field(statePool, 'effectiveDensityInLandIce', effectiveDensityField, 2)
         call mpas_dmpar_exch_halo_field(effectiveDensityField)
         call mpas_timer_stop("se effective density halo")
      end if

      call mpas_timer_stop('se fini')
      call mpas_timer_stop("se timestep")

      deallocate(n_bcl_iter)

   end subroutine ocn_time_integrator_si!}}}

!***********************************************************************
!
!  routine ocn_time_integration_si_init
!
!> \brief   Initialize split-explicit time stepping within MPAS-Ocean core
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine initializes variables required for the split-explicit time
!>  stepper.
!
!-----------------------------------------------------------------------
end module ocn_time_integration_si

! vim: foldmethod=marker
