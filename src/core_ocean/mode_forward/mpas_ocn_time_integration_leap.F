! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_leap
!
!> \brief MPAS ocean Leap-frog Time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   March 2019
!> \details
!>  This module contains the forward time integration routine.
!
!-----------------------------------------------------------------------

module ocn_time_integration_leap

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_threading
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer

   use ocn_constants
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

   public :: ocn_time_integrator_leap

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integrator_leap
!
!> \brief MPAS ocean Leap-frog Time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine integrates one timestep (dt) using an Leap-frog time integrator.
!
!-----------------------------------------------------------------------

   subroutine ocn_time_integrator_leap(domain, dt)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Advance model state forward in time by the specified time step using Leap-frog scheme
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:))
   !                 plus mesh meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information
      real (kind=RKIND), intent(in) :: dt !< Input: timestep
      real (kind=RKIND)             :: dth ! Half of dt

      integer :: iCell, iEdge, k, i, err
      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: tendPool
      type (mpas_pool_type), pointer :: tracersTendPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: provisStatePool
      type (mpas_pool_type), pointer :: provisTracersPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: verticalMeshPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: swForcingPool


      type (mpas_pool_type), pointer :: nextProvisPool, prevProvisPool

      real (kind=RKIND), dimension(4) :: rk_weights, rk_substep_weights

      real (kind=RKIND) :: coef
      real (kind=RKIND), dimension(:,:), pointer :: &
        vertViscTopOfEdge, vertDiffTopOfCell
      real (kind=RKIND), dimension(:), allocatable:: A,C,uTemp
      real (kind=RKIND), dimension(:,:), allocatable:: tracersTemp

      ! Dimensions
      integer, pointer :: nCells, nEdges, nVertLevels, num_tracers

      ! Config options
      logical, pointer :: config_prescribe_velocity, config_prescribe_thickness
      logical, pointer :: config_filter_btr_mode, config_use_freq_filtered_thickness
      logical, pointer :: config_use_standardGM
      logical, pointer :: config_use_cvmix_kpp
      logical, pointer :: config_use_tracerGroup
      real (kind=RKIND), pointer :: config_mom_del4
      character (len=StrKIND), pointer :: config_land_ice_flux_mode

      ! State indices
      integer, pointer :: indexTemperature
      integer, pointer :: indexSalinity

      ! Diagnostics Indices
      integer, pointer :: indexSurfaceVelocityZonal, indexSurfaceVelocityMeridional
      integer, pointer :: indexSSHGradientZonal, indexSSHGradientMeridional

      ! Mesh array pointers
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop

      ! Provis Array Pointers
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityProvis, layerThicknessProvis
      real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessProvis
      real (kind=RKIND), dimension(:,:), pointer :: lowFreqDivergenceProvis
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroupProvis

      ! Tend Array Pointers
      real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessTend, lowFreqDivergenceTend, normalVelocityTend, &
                                                    layerThicknessTend
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroupTend

      ! Diagnostics Array Pointers
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge
      real (kind=RKIND), dimension(:,:), pointer :: vertAleTransportTop
      real (kind=RKIND), dimension(:,:), pointer :: normalTransportVelocity, normalGMBolusVelocity
      real (kind=RKIND), dimension(:,:), pointer :: velocityX, velocityY, velocityZ
      real (kind=RKIND), dimension(:,:), pointer :: velocityZonal, velocityMeridional
      real (kind=RKIND), dimension(:), pointer :: gradSSH
      real (kind=RKIND), dimension(:), pointer :: gradSSHX, gradSSHY, gradSSHZ
      real (kind=RKIND), dimension(:), pointer :: gradSSHZonal, gradSSHMeridional
      real (kind=RKIND), dimension(:,:), pointer :: surfaceVelocity, sshGradient

      ! State Array Pointers
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityCur, normalVelocityNew
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessCur, layerThicknessNew
      real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessCur, highFreqThicknessNew
      real (kind=RKIND), dimension(:,:), pointer :: lowFreqDivergenceCur, lowFreqDivergenceNew
      real (kind=RKIND), dimension(:), pointer :: sshCur, sshNew

      real (kind=RKIND), dimension(:,:,:), pointer :: tracerGroup, tracersCur, tracersNew

      ! Diagnostics Field Pointers
      type (field1DReal), pointer :: boundaryLayerDepthField, effectiveDensityField
      type (field2DReal), pointer :: normalizedRelativeVorticityEdgeField, divergenceField, relativeVorticityField

      ! State/Tend Field Pointers
      type (field2DReal), pointer :: normalVelocityField, layerThicknessField
      type (field3DReal), pointer :: tracersGroupField

      ! Tracer Group Iteartion
      type (mpas_pool_iterator_type) :: groupItr
      character (len=StrKIND) :: modifiedGroupName
      character (len=StrKIND) :: configName

      ! Half of time step size
      dth = 0.5 * dt

      ! Get config options
      call mpas_pool_get_config(domain % configs, 'config_mom_del4', config_mom_del4)
      call mpas_pool_get_config(domain % configs, 'config_filter_btr_mode', config_filter_btr_mode)
      call mpas_pool_get_config(domain % configs, 'config_prescribe_velocity', config_prescribe_velocity)
      call mpas_pool_get_config(domain % configs, 'config_prescribe_thickness', config_prescribe_thickness)
      call mpas_pool_get_config(domain % configs, 'config_use_freq_filtered_thickness', config_use_freq_filtered_thickness)
      call mpas_pool_get_config(domain % configs, 'config_use_standardGM', config_use_standardGM)
      call mpas_pool_get_config(domain % configs, 'config_use_cvmix_kpp', config_use_cvmix_kpp)
      call mpas_pool_get_config(domain % configs, 'config_land_ice_flux_mode', config_land_ice_flux_mode)

      !
      ! Initialize time_levs(2) with state at current time
      ! Couple tracers time_levs(2) with layerThickness in time-levels
      !
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

         call mpas_pool_create_pool(provisStatePool)

         call mpas_pool_clone_pool(statePool, provisStatePool, 1)
         call mpas_pool_add_subpool(block % structs, 'provis_state', provisStatePool)
         call mpas_threading_barrier()

         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

         call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)
         call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceCur, 1)
         call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceNew, 2)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

         !$omp do schedule(runtime) private(k)
         do iEdge = 1, nEdges
            do k = 1, maxLevelEdgeTop(iEdge)
               normalVelocityNew(k, iEdge) = normalVelocityCur(k, iEdge)
            end do
         end do
         !$omp end do

         !$omp do schedule(runtime) private(k)
         do iCell = 1, nCells
            do k = 1, maxLevelCell(iCell)
               layerThicknessNew(k, iCell) = layerThicknessCur(k, iCell)
            end do
         end do
         !$omp end do

         call mpas_pool_begin_iteration(tracersPool)
         do while ( mpas_pool_get_next_member(tracersPool, groupItr) )

            if ( groupItr % memberType == MPAS_POOL_FIELD ) then

               call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersCur, 1)
               call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersNew, 2)

               if ( associated(tracersCur) .and. associated(tracersNew) ) then
                  !$omp do schedule(runtime) private(k)
                  do iCell = 1, nCells  ! couple tracers to thickness
                     do k = 1, maxLevelCell(iCell)
                        tracersNew(:, k, iCell) = tracersCur(:, k, iCell) * layerThicknessCur(k, iCell)
                     end do
                  end do
                  !$omp end do
               end if
            end if
         end do

         if (associated(highFreqThicknessCur)) then
              !$omp do schedule(runtime)
              do iCell = 1, nCells
                 highFreqThicknessNew(:, iCell) = highFreqThicknessCur(:, iCell)
              end do
              !$omp end do
         end if

         if (associated(lowFreqDivergenceCur)) then
              !$omp do schedule(runtime)
              do iCell = 1, nCells
                 lowFreqDivergenceNew(:, iCell) = lowFreqDivergenceCur(:, iCell)
              end do
              !$omp end do
         end if

         block => block % next
      end do

      block => domain % blocklist
      do while(associated(block))
         if (associated(block % prev)) then
            call mpas_pool_get_subpool(block % prev % structs, 'provis_state', prevProvisPool)
         else
            nullify(prevProvisPool)
         end if

         if (associated(block % next)) then
            call mpas_pool_get_subpool(block % next % structs, 'provis_state', nextProvisPool)
         else
            nullify(nextProvisPool)
         end if

         call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

         if (associated(prevProvisPool) .and. associated(nextProvisPool)) then
            call mpas_pool_link_pools(provisStatePool, prevProvisPool, nextProvisPool)
         else if (associated(prevProvisPool)) then
            call mpas_pool_link_pools(provisStatePool, prevProvisPool)
         else if (associated(nextProvisPool)) then
            call mpas_pool_link_pools(provisStatePool, nextPool=nextProvisPool)
         else
            call mpas_pool_link_pools(provisStatePool)
         end if

         call mpas_pool_link_parinfo(block, provisStatePool)

         block => block % next
      end do

      call mpas_threading_barrier()


      call mpas_timer_start("LEAP-main ")


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN Leap-frog time integration
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-----------------------------------------------------------------------!     
! Step 1 | Compute tendency : u^(n) -> F(u^(n))
!-----------------------------------------------------------------------!     

        call mpas_pool_get_subpool(domain % blocklist % structs, 'diagnostics', diagnosticsPool)

        ! Update halos for diagnostic variables.
        if (config_use_cvmix_kpp) then
           call mpas_timer_start("LEAP-boundary layer depth halo update")
           call mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
           call mpas_timer_stop("LEAP-boundary layer depth halo update")
        end if

        call mpas_timer_start("LEAP-diagnostic halo update")

        call mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
        if (config_mom_del4 > 0.0_RKIND) then
           call mpas_dmpar_field_halo_exch(domain, 'divergence')
           call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
        end if
        call mpas_timer_stop("LEAP-diagnostic halo update")
        call mpas_threading_barrier()

        !----------------------------------------------------------------------------------------------------------------!

        ! Compute tendencies for high frequency thickness

        if (config_use_freq_filtered_thickness) then
           call mpas_timer_start("LEAP-tendency computations")

           block => domain % blocklist
           do while (associated(block))
              call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
              call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
              call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
              call mpas_pool_get_subpool(block % structs, 'state', statePool)
              call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
              call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

              call ocn_tend_freq_filtered_thickness(tendPool, provisStatePool, diagnosticsPool, meshPool, 1)
              call mpas_threading_barrier()
              block => block % next
           end do

           call mpas_timer_stop("LEAP-tendency computations")

           call mpas_timer_start("LEAP-prognostic halo update")

           call mpas_dmpar_field_halo_exch(domain, 'tendHighFreqThickness')
           call mpas_dmpar_field_halo_exch(domain, 'tendLowFreqDivergence')

           call mpas_timer_stop("LEAP-prognostic halo update")
           call mpas_threading_barrier()


           ! Compute next substep state for high frequency thickness.
           block => domain % blocklist
           do while (associated(block))
              call mpas_pool_get_subpool(block % structs, 'state', statePool)
              call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
              call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
              call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

              call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
              call mpas_pool_get_array(provisStatePool, 'highFreqThickness', highFreqThicknessProvis, 1)
              call mpas_pool_get_array(tendPool, 'highFreqThickness', highFreqThicknessTend)

              call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)

              !$omp do schedule(runtime)
              do iCell = 1, nCells
                 highFreqThicknessProvis(:, iCell) = highFreqThicknessCur(:, iCell) + dth * highFreqThicknessTend(:, iCell)
              end do
              !$omp end do
              call mpas_threading_barrier()
              block => block % next
           end do

        endif


        ! Compute tendencies for velocity, thickness, and tracers.
        ! In RK4 notation, we are computing the right hand side f(t,y),
        ! which is the same as k_j / h.
        call mpas_timer_start("LEAP vel/thick tendency computations")

        block => domain % blocklist
        do while (associated(block))
           call ocn_time_integrator_leap_compute_vel_tends( block, dth, err )

           call ocn_time_integrator_leap_compute_thick_tends( block, dth, err )
           block => block % next
        end do

        call mpas_timer_stop("LEAP vel/thick tendency computations")

        ! Update halos for prognostic variables.

        call mpas_timer_start("LEAP vel/thick prognostic halo update")

        call mpas_dmpar_field_halo_exch(domain, 'tendNormalVelocity')
        call mpas_dmpar_field_halo_exch(domain, 'tendLayerThickness')

        call mpas_timer_stop("LEAP vel/thick prognostic halo update")

        call mpas_timer_start("LEAP tracer tendency computations")

        block => domain % blocklist
        do while (associated(block))
           call ocn_time_integrator_leap_compute_tracer_tends( block, dth, err )
           block => block % next
        end do

        call mpas_timer_stop("LEAP tracer tendency computations")

        call mpas_timer_start("LEAP tracer prognostic halo update")

        call mpas_pool_get_subpool(domain % blocklist % structs, 'tend', tendPool)
        call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

        call mpas_pool_begin_iteration(tracersTendPool)
        do while ( mpas_pool_get_next_member(tracersTendPool, groupItr) )
           if ( groupItr % memberType == MPAS_POOL_FIELD ) then
              call mpas_dmpar_field_halo_exch(domain, trim(groupItr % memberName))
           end if
        end do

        call mpas_timer_stop("LEAP tracer prognostic halo update")
        call mpas_threading_barrier()


!-----------------------------------------------------------------------------------------!     
! Step 2 | Forward half dt : u^(n+1/2) = u^(n)+(dt/2)*F[u^(n)]
!-----------------------------------------------------------------------------------------!     

        ! Diagnostic update

        call mpas_timer_start("LEAP-update diagnostic variables")

           block => domain % blocklist
           do while (associated(block))
              call ocn_time_integrator_leap_diagnostic_update(block, dth, err)
              block => block % next
           end do

        call mpas_timer_stop("LEAP-update diagnostic variables")
        call mpas_threading_barrier()


!-----------------------------------------------------------------------------------------!     
! Step 3 | Compute tendency at (n+1/2) : u^(n+1/2) -> F[u^(n+1/2)]
!-----------------------------------------------------------------------------------------!     

        call mpas_pool_get_subpool(domain % blocklist % structs, 'diagnostics', diagnosticsPool)

        ! Update halos for diagnostic variables.
        if (config_use_cvmix_kpp) then
           call mpas_timer_start("LEAP-boundary layer depth halo update")
           call mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
           call mpas_timer_stop("LEAP-boundary layer depth halo update")
        end if

        call mpas_timer_start("LEAP-diagnostic halo update")

        call mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
        if (config_mom_del4 > 0.0_RKIND) then
           call mpas_dmpar_field_halo_exch(domain, 'divergence')
           call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
        end if
        call mpas_timer_stop("LEAP-diagnostic halo update")
        call mpas_threading_barrier()

        !----------------------------------------------------------------------------------------------------------------!

        ! Compute tendencies for high frequency thickness

        if (config_use_freq_filtered_thickness) then
           call mpas_timer_start("LEAP-tendency computations")

           block => domain % blocklist
           do while (associated(block))
              call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
              call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
              call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
              call mpas_pool_get_subpool(block % structs, 'state', statePool)
              call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
              call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

              call ocn_tend_freq_filtered_thickness(tendPool, provisStatePool, diagnosticsPool, meshPool, 1)
              call mpas_threading_barrier()
              block => block % next
           end do

           call mpas_timer_stop("LEAP-tendency computations")

           call mpas_timer_start("LEAP-prognostic halo update")

           call mpas_dmpar_field_halo_exch(domain, 'tendHighFreqThickness')
           call mpas_dmpar_field_halo_exch(domain, 'tendLowFreqDivergence')

           call mpas_timer_stop("LEAP-prognostic halo update")
           call mpas_threading_barrier()

        endif

        ! Compute tendencies for velocity, thickness, and tracers.
        call mpas_timer_start("LEAP vel/thick tendency computations")

        block => domain % blocklist
        do while (associated(block))
           call ocn_time_integrator_leap_compute_vel_tends( block, dth, err )

           call ocn_time_integrator_leap_compute_thick_tends( block, dth, err )
           block => block % next
        end do

        call mpas_timer_stop("LEAP vel/thick tendency computations")

        ! Update halos for prognostic variables.

        call mpas_timer_start("LEAP vel/thick prognostic halo update")

        call mpas_dmpar_field_halo_exch(domain, 'tendNormalVelocity')
        call mpas_dmpar_field_halo_exch(domain, 'tendLayerThickness')

        call mpas_timer_stop("LEAP vel/thick prognostic halo update")

        call mpas_timer_start("LEAP tracer tendency computations")

        block => domain % blocklist
        do while (associated(block))
           call ocn_time_integrator_leap_compute_tracer_tends( block, dth, err )
           block => block % next
        end do

        call mpas_timer_stop("LEAP tracer tendency computations")

        call mpas_timer_start("LEAP tracer prognostic halo update")

        call mpas_pool_get_subpool(domain % blocklist % structs, 'tend', tendPool)
        call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

        call mpas_pool_begin_iteration(tracersTendPool)
        do while ( mpas_pool_get_next_member(tracersTendPool, groupItr) )
           if ( groupItr % memberType == MPAS_POOL_FIELD ) then
              call mpas_dmpar_field_halo_exch(domain, trim(groupItr % memberName))
           end if
        end do

        call mpas_timer_stop("LEAP tracer prognostic halo update")
        call mpas_threading_barrier()
 
!-----------------------------------------------------------------------------------------!     
! Step 4 | Leap-frog for dt : u(n+1) = u^(n) + F[u^(n+1/2)]
!-----------------------------------------------------------------------------------------!     

        ! Accumulate update.
        call mpas_timer_start("LEAP accumulate update")

        block => domain % blocklist
        do while (associated(block))
           call ocn_time_integrator_leap_accumulate_update(block, dt, err)
           block => block % next
        end do

        call mpas_timer_stop("LEAP accumulate update")
        call mpas_threading_barrier()

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END Leap-frog time integration
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call mpas_timer_stop("LEAP-main ")
      call mpas_threading_barrier()


      !
      !  A little clean up at the end: rescale tracer fields and compute diagnostics for new state
      !
      call mpas_timer_start("LEAP-cleanup phase")

      ! Rescale tracers
      block => domain % blocklist
      do while(associated(block))
        call ocn_time_integrator_leap_cleanup(block, dt, err)

        block => block % next
      end do

      call mpas_timer_stop("LEAP-cleanup phase")

      call mpas_timer_start("LEAP-implicit vert mix")
      ! Update halo on u and tracers, which were just updated for implicit vertical mixing.  If not done,
      ! this leads to lack of volume conservation.  It is required because halo updates in RK4 are only
      ! conducted on tendencies, not on the velocity and tracer fields.  So this update is required to
      ! communicate the change due to implicit vertical mixing across the boundary.
      call mpas_timer_start("LEAP-implicit vert mix halos")

      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=2)

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_dmpar_field_halo_exch(domain, groupItr % memberName, timeLevel=2)
         end if
      end do

      call mpas_timer_stop("LEAP-implicit vert mix halos")

      call mpas_timer_stop("LEAP-implicit vert mix")
      call mpas_threading_barrier()

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

         call mpas_pool_get_dimension(diagnosticsPool, 'index_surfaceVelocityZonal', indexSurfaceVelocityZonal)
         call mpas_pool_get_dimension(diagnosticsPool, 'index_surfaceVelocityMeridional', indexSurfaceVelocityMeridional)
         call mpas_pool_get_dimension(diagnosticsPool, 'index_SSHGradientZonal', indexSSHGradientZonal)
         call mpas_pool_get_dimension(diagnosticsPool, 'index_SSHGradientMeridional', indexSSHGradientMeridional)

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
         call mpas_threading_barrier()

         ! Update the effective desnity in land ice if we're coupling to land ice
         call ocn_effective_density_in_land_ice_update(meshPool, forcingPool, statePool, scratchPool, err)

         ! ------------------------------------------------------------------
         ! Accumulating various parameterizations of the transport velocity
         ! ------------------------------------------------------------------
         !$omp do schedule(runtime)
         do iEdge = 1, nEdges
            normalTransportVelocity(:, iEdge) = normalVelocityNew(:, iEdge)
         end do
         !$omp end do
         call mpas_threading_barrier()

         ! Compute normalGMBolusVelocity and the tracer transport velocity
         if (config_use_standardGM) then
             call ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, scratchPool)
         end if
         call mpas_threading_barrier()

         if (config_use_standardGM) then
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               normalTransportVelocity(:, iEdge) = normalTransportVelocity(:, iEdge) + normalGMBolusVelocity(:, iEdge)
            end do
            !$omp end do
         end if
         ! ------------------------------------------------------------------
         ! End: Accumulating various parameterizations of the transport velocity
         ! ------------------------------------------------------------------

         call mpas_reconstruct(meshPool,  normalVelocityNew, &
                          velocityX, velocityY, velocityZ,   &
                          velocityZonal, velocityMeridional, &
                          includeHalos = .true.)

         call mpas_reconstruct(meshPool, gradSSH,          &
                          gradSSHX, gradSSHY, gradSSHZ,    &
                          gradSSHZonal, gradSSHMeridional, &
                          includeHalos = .true.)
         call mpas_threading_barrier()

         !$omp do schedule(runtime)
         do iCell = 1, nCells
            surfaceVelocity(indexSurfaceVelocityZonal, iCell) = velocityZonal(1, iCell)
            surfaceVelocity(indexSurfaceVelocityMeridional, iCell) = velocityMeridional(1, iCell)

            SSHGradient(indexSSHGradientZonal, iCell) = gradSSHZonal(iCell)
            SSHGradient(indexSSHGradientMeridional, iCell) = gradSSHMeridional(iCell)
         end do
         !$omp end do

         call ocn_time_average_coupled_accumulate(diagnosticsPool, statePool, forcingPool, 2)

         if (config_use_standardGM) then
            call ocn_reconstruct_gm_vectors(diagnosticsPool, meshPool)
         end if
         call mpas_threading_barrier()

         block => block % next
      end do

      if (trim(config_land_ice_flux_mode) == 'coupled') then
         call mpas_timer_start("LEAP-effective density halo")
         call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
         call mpas_pool_get_field(statePool, 'effectiveDensityInLandIce', effectiveDensityField, 2)
         call mpas_dmpar_exch_halo_field(effectiveDensityField)
         call mpas_timer_stop("LEAP-effective density halo")
      end if


      call mpas_threading_barrier()

      block => domain % blocklist
      do while(associated(block))
         call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

         call mpas_pool_destroy_pool(provisStatePool)

         call mpas_pool_remove_subpool(block % structs, 'provis_state')
         block => block % next
      end do
      call mpas_threading_barrier()

   end subroutine ocn_time_integrator_leap!}}}

!***********************************************************************************!   
!***********************************************************************************! 
!***********************************************************************************! fend  

   subroutine ocn_time_integrator_leap_compute_vel_tends(block, dt, err)!{{{
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool
      type (mpas_pool_type), pointer :: statePool, diagnosticsPool, forcingPool
      type (mpas_pool_type), pointer :: scratchPool, tendPool, provisStatePool
      type (mpas_pool_type), pointer :: tracersPool

      real (kind=RKIND), dimension(:), pointer :: sshCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCur, normalVelocityCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessEdge, vertAleTransportTop
      real (kind=RKIND), dimension(:, :), pointer :: normalTransportVelocity
      real (kind=RKIND), dimension(:, :), pointer ::  normalVelocityProvis, highFreqThicknessProvis

      logical, pointer :: config_filter_btr_mode

      err = 0

      call mpas_pool_get_config(block % configs, 'config_filter_btr_mode', config_filter_btr_mode)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
      call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)

      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)
      call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)

      call mpas_pool_get_array(provisStatePool, 'normalVelocity', normalVelocityProvis, 1)
      call mpas_pool_get_array(provisStatePool, 'highFreqThickness', highFreqThicknessProvis, 1)

      ! advection of u uses u, while advection of layerThickness and tracers use normalTransportVelocity.
      if (associated(highFreqThicknessProvis)) then
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur,layerThicknessEdge, normalVelocityProvis, &
            sshCur, dt, &
            vertAleTransportTop, err, highFreqThicknessProvis)
      else
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur,layerThicknessEdge, normalVelocityProvis, &
            sshCur, dt, &
            vertAleTransportTop, err)
      endif
      call mpas_threading_barrier()

      call ocn_tend_vel(tendPool, provisStatePool, forcingPool, diagnosticsPool, meshPool, scratchPool, 1)
      call mpas_threading_barrier()

   end subroutine ocn_time_integrator_leap_compute_vel_tends!}}}

!***********************************************************************************!   
!***********************************************************************************!   

   subroutine ocn_time_integrator_leap_compute_thick_tends(block, dt, err)!{{{
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool
      type (mpas_pool_type), pointer :: statePool, diagnosticsPool, forcingPool
      type (mpas_pool_type), pointer :: scratchPool, tendPool, provisStatePool
      type (mpas_pool_type), pointer :: tracersPool

      real (kind=RKIND), dimension(:), pointer :: sshCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCur, normalVelocityCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessEdge, vertAleTransportTop
      real (kind=RKIND), dimension(:, :), pointer :: normalTransportVelocity
      real (kind=RKIND), dimension(:, :), pointer ::  normalVelocityProvis, highFreqThicknessProvis

      logical, pointer :: config_filter_btr_mode

      err = 0

      call mpas_pool_get_config(block % configs, 'config_filter_btr_mode', config_filter_btr_mode)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
      call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)

      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)
      call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)

      call mpas_pool_get_array(provisStatePool, 'normalVelocity', normalVelocityProvis, 1)
      call mpas_pool_get_array(provisStatePool, 'highFreqThickness', highFreqThicknessProvis, 1)

      ! advection of u uses u, while advection of layerThickness and tracers use normalTransportVelocity.
      if (associated(highFreqThicknessProvis)) then
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur, layerThicknessEdge, normalTransportVelocity, &
            sshCur, dt, &
            vertAleTransportTop, err, highFreqThicknessProvis)
      else
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur, layerThicknessEdge, normalTransportVelocity, &
            sshCur, dt, &
            vertAleTransportTop, err)
      endif
      call mpas_threading_barrier()

      call ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool)

      call mpas_threading_barrier()

   end subroutine ocn_time_integrator_leap_compute_thick_tends!}}}

!***********************************************************************************!   
!***********************************************************************************!   

   subroutine ocn_time_integrator_leap_compute_tracer_tends(block, dt, err)!{{{
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool
      type (mpas_pool_type), pointer :: statePool, diagnosticsPool, forcingPool
      type (mpas_pool_type), pointer :: scratchPool, tendPool, provisStatePool
      type (mpas_pool_type), pointer :: swForcingPool, tracersPool

      real (kind=RKIND), dimension(:), pointer :: sshCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCur, normalVelocityCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessEdge, vertAleTransportTop
      real (kind=RKIND), dimension(:, :), pointer :: normalTransportVelocity
      real (kind=RKIND), dimension(:, :), pointer ::  normalVelocityProvis, highFreqThicknessProvis

      logical, pointer :: config_filter_btr_mode

      err = 0

      call mpas_pool_get_config(block % configs, 'config_filter_btr_mode', config_filter_btr_mode)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)
      call mpas_pool_get_subpool(block % structs, 'shortwave', swForcingPool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
      call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)

      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)
      call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)

      call mpas_pool_get_array(provisStatePool, 'normalVelocity', normalVelocityProvis, 1)
      call mpas_pool_get_array(provisStatePool, 'highFreqThickness', highFreqThicknessProvis, 1)

      ! advection of u uses u, while advection of layerThickness and tracers use normalTransportVelocity.
      if (associated(highFreqThicknessProvis)) then
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur, layerThicknessEdge, normalTransportVelocity, &
            sshCur, dt, &
            vertAleTransportTop, err, highFreqThicknessProvis)
      else
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur, layerThicknessEdge, normalTransportVelocity, &
            sshCur, dt, &
            vertAleTransportTop, err)
      endif
      call mpas_threading_barrier()

      if (config_filter_btr_mode) then
          call ocn_filter_btr_mode_tend_vel(tendPool, provisStatePool, diagnosticsPool, meshPool, 1)
      endif
      call mpas_threading_barrier()

      call ocn_tend_tracer(tendPool, provisStatePool, forcingPool, diagnosticsPool, meshPool, swForcingPool, &
                           scratchPool, dt, activeTracersOnlyIn=.false., timeLevelIn=1)
      call mpas_threading_barrier()

   end subroutine ocn_time_integrator_leap_compute_tracer_tends!}}}

!***********************************************************************************!   
!***********************************************************************************!   

   subroutine ocn_time_integrator_leap_compute_tends_prov(block, dt, err)!{{{
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool
      type (mpas_pool_type), pointer :: statePool, diagnosticsPool, forcingPool
      type (mpas_pool_type), pointer :: scratchPool, tendPool, provisStatePool
      type (mpas_pool_type), pointer :: swForcingPool, tracersPool

      real (kind=RKIND), dimension(:), pointer :: sshCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCur, normalVelocityCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessEdge, vertAleTransportTop
      real (kind=RKIND), dimension(:, :), pointer :: normalTransportVelocity
      real (kind=RKIND), dimension(:, :), pointer ::  normalVelocityProvis, highFreqThicknessProvis

      logical, pointer :: config_filter_btr_mode

      err = 0

      call mpas_pool_get_config(block % configs, 'config_filter_btr_mode', config_filter_btr_mode)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)
      call mpas_pool_get_subpool(block % structs, 'shortwave', swForcingPool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
      call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)

      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)
      call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)

      call mpas_pool_get_array(provisStatePool, 'normalVelocity', normalVelocityProvis, 1)
      call mpas_pool_get_array(provisStatePool, 'highFreqThickness', highFreqThicknessProvis, 1)


      ! advection of u uses u, while advection of layerThickness and tracers use normalTransportVelocity.
      if (associated(highFreqThicknessProvis)) then
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur,layerThicknessEdge, normalVelocityProvis, &
            sshCur, dt, &
            vertAleTransportTop, err, highFreqThicknessProvis)
      else
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur,layerThicknessEdge, normalVelocityProvis, &
            sshCur, dt, &
            vertAleTransportTop, err)
      endif
      call mpas_threading_barrier()

      call ocn_tend_vel(tendPool, provisStatePool, forcingPool, diagnosticsPool, meshPool, scratchPool, 1)
      call mpas_threading_barrier()

      if (associated(highFreqThicknessProvis)) then
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur, layerThicknessEdge, normalTransportVelocity, &
            sshCur, dt, &
            vertAleTransportTop, err, highFreqThicknessProvis)
      else
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessCur, layerThicknessEdge, normalTransportVelocity, &
            sshCur, dt, &
            vertAleTransportTop, err)
      endif
      call mpas_threading_barrier()

      call ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool)

      if (config_filter_btr_mode) then
          call ocn_filter_btr_mode_tend_vel(tendPool, provisStatePool, diagnosticsPool, meshPool, 1)
      endif
      call mpas_threading_barrier()

      call ocn_tend_tracer(tendPool, provisStatePool, forcingPool, diagnosticsPool, meshPool, swForcingPool, &
                           scratchPool, dt, activeTracersOnlyIn=.false., timeLevelIn=1)
      call mpas_threading_barrier()

   end subroutine ocn_time_integrator_leap_compute_tends_prov!}}}

!***********************************************************************************!   
!***********************************************************************************!   

   subroutine ocn_time_integrator_leap_compute_tends_new(block, dt, err)!{{{
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool
      type (mpas_pool_type), pointer :: statePool, diagnosticsPool, forcingPool
      type (mpas_pool_type), pointer :: scratchPool, tendPool, provisStatePool
      type (mpas_pool_type), pointer :: swForcingPool, tracersPool

      real (kind=RKIND), dimension(:), pointer :: sshCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCur, normalVelocityCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessEdge, vertAleTransportTop
      real (kind=RKIND), dimension(:, :), pointer :: normalTransportVelocity
      real (kind=RKIND), dimension(:, :), pointer :: normalVelocityProvis, highFreqThicknessProvis
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessProvis

      logical, pointer :: config_filter_btr_mode

      err = 0

      call mpas_pool_get_config(block % configs, 'config_filter_btr_mode', config_filter_btr_mode)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)
      call mpas_pool_get_subpool(block % structs, 'shortwave', swForcingPool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
      call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)

      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)
      call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)

      call mpas_pool_get_array(provisStatePool, 'normalVelocity', normalVelocityProvis, 1)
      call mpas_pool_get_array(provisStatePool, 'highFreqThickness', highFreqThicknessProvis, 1)
      call mpas_pool_get_array(provisStatePool,    'layerThickness',    layerThicknessProvis, 1)


      ! advection of u uses u, while advection of layerThickness and tracers use normalTransportVelocity.
      if (associated(highFreqThicknessProvis)) then
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessProvis,layerThicknessEdge, normalVelocityProvis, &
            sshCur, dt, &
            vertAleTransportTop, err, highFreqThicknessProvis)
      else
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessProvis,layerThicknessEdge, normalVelocityProvis, &
            sshCur, dt, &
            vertAleTransportTop, err)
      endif
      call mpas_threading_barrier()

      call ocn_tend_vel(tendPool, provisStatePool, forcingPool, diagnosticsPool, meshPool, scratchPool, 1)
      call mpas_threading_barrier()

      if (associated(highFreqThicknessProvis)) then
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessProvis, layerThicknessEdge, normalTransportVelocity, &
            sshCur, dt, &
            vertAleTransportTop, err, highFreqThicknessProvis)
      else
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThicknessProvis, layerThicknessEdge, normalTransportVelocity, &
            sshCur, dt, &
            vertAleTransportTop, err)
      endif
      call mpas_threading_barrier()

      call ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool)

      if (config_filter_btr_mode) then
          call ocn_filter_btr_mode_tend_vel(tendPool, provisStatePool, diagnosticsPool, meshPool, 1)
      endif
      call mpas_threading_barrier()

      call ocn_tend_tracer(tendPool, provisStatePool, forcingPool, diagnosticsPool, meshPool, swForcingPool, &
                           scratchPool, dt, activeTracersOnlyIn=.false., timeLevelIn=1)
      call mpas_threading_barrier()

   end subroutine ocn_time_integrator_leap_compute_tends_new!}}}

!***********************************************************************************!   
!***********************************************************************************!   

   subroutine ocn_time_integrator_leap_diagnostic_update(block, dt, err)!{{{
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err

      logical, pointer :: config_prescribe_velocity, config_prescribe_thickness, config_use_standardGM

      integer, pointer :: nCells, nEdges
      integer :: iCell, iEdge, k

      type (mpas_pool_type), pointer :: statePool, tendPool, meshPool, scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool, provisStatePool, forcingPool
      type (mpas_pool_type), pointer :: tracersPool, tracersTendPool, provisTracersPool

      real (kind=RKIND), dimension(:, :), pointer :: normalVelocityCur, normalVelocityProvis, normalVelocityTend
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCur, layerThicknessProvis, layerThicknessTend
      real (kind=RKIND), dimension(:, :), pointer :: lowFreqDivergenceCur, lowFreqDivergenceProvis, lowFreqDivergenceTend
      real (kind=RKIND), dimension(:, :), pointer :: normalTransportVelocity, normalGMBolusVelocity

      real (kind=RKIND), dimension(:, :, :), pointer :: tracersGroupCur, tracersGroupProvis, tracersGroupTend

      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop

      logical, pointer :: config_use_tracerGroup
      type (mpas_pool_iterator_type) :: groupItr
      character (len=StrKIND) :: modifiedGroupName
      character (len=StrKIND) :: configName

      err = 0

      call mpas_pool_get_config(block % configs, 'config_prescribe_velocity', config_prescribe_velocity)
      call mpas_pool_get_config(block % configs, 'config_prescribe_thickness', config_prescribe_thickness)
      call mpas_pool_get_config(block % configs, 'config_use_standardGM', config_use_standardGM)

      call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

      call mpas_pool_get_subpool(provisStatePool, 'tracers', provisTracersPool)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
      call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceCur, 1)

      call mpas_pool_get_array(provisStatePool, 'normalVelocity', normalVelocityProvis, 1)
      call mpas_pool_get_array(provisStatePool, 'layerThickness', layerThicknessProvis, 1)
      call mpas_pool_get_array(provisStatePool, 'lowFreqDivergence', lowFreqDivergenceProvis, 1)

      call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
      call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)
      call mpas_pool_get_array(tendPool, 'lowFreqDivergence', lowFreqDivergenceTend)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

      call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)

      call mpas_threading_barrier()

      !$omp do schedule(runtime) private(k)
      do iEdge = 1, nEdges
         do k = 1, maxLevelEdgeTop(iEdge)
            normalVelocityProvis(k, iEdge) = normalVelocityCur(k, iEdge) + dt &
                                           * normalVelocityTend(k, iEdge)
         end do
      end do
      !$omp end do


      !$omp do schedule(runtime) private(k)
      do iCell = 1, nCells
         do k = 1, maxLevelCell(iCell)
            layerThicknessProvis(k, iCell) = layerThicknessCur(k, iCell) + dt &
                                           * layerThicknessTend(k, iCell)
         end do
      end do
      !$omp end do

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            configName = 'config_use_' // trim(groupItr % memberName)
            call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)

            if ( config_use_tracerGroup ) then
               call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupCur, 1)
               call mpas_pool_get_array(provisTracersPool, groupItr % memberName, tracersGroupProvis, 1)

               modifiedGroupName = trim(groupItr % memberName) // 'Tend'
               call mpas_pool_get_array(tracersTendPool, modifiedGroupName, tracersGroupTend)
               if ( associated(tracersGroupProvis) .and. associated(tracersGroupCur) .and. associated(tracersGroupTend) ) then
                  !$omp do schedule(runtime) private(k)
                  do iCell = 1, nCells
                     do k = 1, maxLevelCell(iCell)
                        tracersGroupProvis(:, k, iCell) = ( layerThicknessCur(k, iCell) * tracersGroupCur(:, k, iCell)  &
                                                 + dt * tracersGroupTend(:, k, iCell) &
                                                   ) / layerThicknessProvis(k, iCell)
                     end do

                  end do
                  !$omp end do
               end if
            end if
         end if
      end do

      if (associated(lowFreqDivergenceCur)) then
         !$omp do schedule(runtime)
         do iCell = 1, nCells
            lowFreqDivergenceProvis(:, iCell) = lowFreqDivergenceCur(:, iCell) + dt &
                                              * lowFreqDivergenceTend(:, iCell)
         end do
         !$omp end do
      end if

      if (config_prescribe_velocity) then
         !$omp do schedule(runtime)
         do iEdge = 1, nEdges
            normalVelocityProvis(:, iEdge) = normalVelocityCur(:, iEdge)
         end do
         !$omp end do
      end if

      if (config_prescribe_thickness) then
         !$omp do schedule(runtime)
         do iCell = 1, nCells
            layerThicknessProvis(:, iCell) = layerThicknessCur(:, iCell)
         end do
         !$omp end do
      end if
      call mpas_threading_barrier()

      call ocn_diagnostic_solve(dt, provisStatePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 1)
      call mpas_threading_barrier()

      ! ------------------------------------------------------------------
      ! Accumulating various parametrizations of the transport velocity
      ! ------------------------------------------------------------------
      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
         normalTransportVelocity(:, iEdge) = normalVelocityProvis(:, iEdge)
      end do
      !$omp end do
      call mpas_threading_barrier()

      ! Compute normalGMBolusVelocity, relativeSlope and RediDiffVertCoef if respective flags are turned on
      if (config_use_standardGM) then
         call ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, scratchPool)
      end if
      call mpas_threading_barrier()

      if (config_use_standardGM) then
         !$omp do schedule(runtime)
         do iEdge = 1, nEdges
            normalTransportVelocity(:, iEdge) = normalTransportVelocity(:, iEdge) + normalGMBolusVelocity(:,iEdge)
         end do
         !$omp end do
      end if
      call mpas_threading_barrier()
      ! ------------------------------------------------------------------
      ! End: Accumulating various parametrizations of the transport velocity
      ! ------------------------------------------------------------------

   end subroutine ocn_time_integrator_leap_diagnostic_update!}}}

!***********************************************************************************!   
!***********************************************************************************!   

   subroutine ocn_time_integrator_leap_accumulate_update(block, dt, err)!{{{
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err

      integer, pointer :: nCells, nEdges
      integer :: iCell, iEdge, k

      type (mpas_pool_type), pointer :: statePool, tendPool, meshPool
      type (mpas_pool_type), pointer :: tracersPool, tracersTendPool

      real (kind=RKIND), dimension(:, :), pointer :: normalVelocityNew, normalVelocityTend
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessNew, layerThicknessTend
      real (kind=RKIND), dimension(:, :), pointer :: highFreqThicknessNew, highFreqThicknessTend
      real (kind=RKIND), dimension(:, :), pointer :: lowFreqDivergenceNew, lowFreqDivergenceTend
!     real (kind=RKIND), dimension(:, :), pointer :: normalVelocityCur
!     real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCur
!     real (kind=RKIND), dimension(:, :), pointer :: highFreqThicknessCur
!     real (kind=RKIND), dimension(:, :), pointer :: lowFreqDivergenceCur

      real (kind=RKIND), dimension(:, :, :), pointer :: tracersGroupNew, tracersGroupTend

      integer, dimension(:), pointer :: maxLevelCell

      logical, pointer :: config_use_tracerGroup
      type (mpas_pool_iterator_type) :: groupItr
      character (len=StrKIND) :: modifiedGroupName
      character (len=StrKIND) :: configName

      err = 0

      call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

!     call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
!     call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
!     call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
!     call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceCur, 1)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)
      call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)
      call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceNew, 2)

      call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
      call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)

      call mpas_pool_get_array(tendPool, 'highFreqThickness', highFreqThicknessTend)
      call mpas_pool_get_array(tendPool, 'lowFreqDivergence', lowFreqDivergenceTend)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      !$omp do schedule(runtime) private(k)
      do iCell = 1, nCells
         do k = 1, maxLevelCell(iCell)
            layerThicknessNew(k, iCell) = layerThicknessNew(k, iCell) + dt * layerThicknessTend(k, iCell)
         end do
      end do
      !$omp end do

      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
         normalVelocityNew(:, iEdge) = normalVelocityNew(:, iEdge) + dt * normalVelocityTend(:, iEdge)
      end do
      !$omp end do

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            configName = 'config_use_' // trim(groupItr % memberName)
            call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)

            if ( config_use_tracerGroup ) then
               call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupNew, 2)

               modifiedGroupName = trim(groupItr % memberName) // 'Tend'
               call mpas_pool_get_array(tracersTendPool, modifiedGroupName, tracersGroupTend)
               if ( associated(tracersGroupNew) .and. associated(tracersGroupTend) ) then
                  !$omp do schedule(runtime) private(k)
                  do iCell = 1, nCells
                     do k = 1, maxLevelCell(iCell)
                        tracersGroupNew(:, k, iCell) = tracersGroupNew(:, k, iCell) + dt &
                                                * tracersGroupTend(:, k, iCell)
                     end do
                  end do
                  !$omp end do
               end if
            end if
         end if
      end do

      if (associated(highFreqThicknessNew)) then
         !$omp do schedule(runtime)
         do iCell = 1, nCells
            highFreqThicknessNew(:, iCell) = highFreqThicknessNew(:, iCell) + dt * highFreqThicknessTend(:, iCell)
         end do
         !$omp end do
      end if

      if (associated(lowFreqDivergenceNew)) then
         !$omp do schedule(runtime)
         do iCell = 1, nCells
            lowFreqDivergenceNew(:, iCell) = lowFreqDivergenceNew(:, iCell) + dt * lowFreqDivergenceTend(:, iCell)
         end do
         !$omp end do
      end if

   end subroutine ocn_time_integrator_leap_accumulate_update!}}}

!***********************************************************************************!   
!***********************************************************************************!   

   subroutine ocn_time_integrator_leap_timefilter(block, dt, err)!{{{
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err

      integer, pointer :: nCells, nEdges
      integer :: iCell, iEdge, k

      type (mpas_pool_type), pointer :: statePool, tendPool,meshPool,provisStatePool
      type (mpas_pool_type), pointer :: tracersPool, tracersTendPool

      real (kind=RKIND), dimension(:, :), pointer :: normalVelocityNew, normalVelocityCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessNew, layerThicknessCur
      real (kind=RKIND), dimension(:, :), pointer :: highFreqThicknessNew, highFreqThicknessCur
      real (kind=RKIND), dimension(:, :), pointer :: lowFreqDivergenceNew, lowFreqDivergenceCur
      real (kind=RKIND), dimension(:, :), pointer :: normalVelocityProvis
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessProvis
      real (kind=RKIND), dimension(:, :), pointer :: highFreqThicknessProvis
      real (kind=RKIND), dimension(:, :), pointer :: lowFreqDivergenceProvis

      real (kind=RKIND), dimension(:, :, :), pointer :: tracersGroupNew, tracersGroupCur

      integer, dimension(:), pointer :: maxLevelCell

      logical, pointer :: config_use_tracerGroup
      type (mpas_pool_iterator_type) :: groupItr
      character (len=StrKIND) :: modifiedGroupName
      character (len=StrKIND) :: configName

      err = 0

      call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
      call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
      call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceCur, 1)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)
      call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)
      call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceNew, 2)

      call mpas_pool_get_array(provisStatePool, 'normalVelocity', normalVelocityProvis, 1)
      call mpas_pool_get_array(provisStatePool, 'layerThickness', layerThicknessProvis, 1)
      call mpas_pool_get_array(provisStatePool, 'highFreqThickness', highFreqThicknessProvis, 1)
      call mpas_pool_get_array(provisStatePool, 'lowFreqDivergence', lowFreqDivergenceProvis, 1)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      !$omp do schedule(runtime) private(k)
      do iCell = 1, nCells
         do k = 1, maxLevelCell(iCell)
            layerThicknessNew(k, iCell) = layerThicknessCur(k, iCell)
         end do
      end do
      !$omp end do

      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
         normalVelocityNew(:, iEdge) = normalVelocityCur(:, iEdge) 
      end do
      !$omp end do

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            configName = 'config_use_' // trim(groupItr % memberName)
            call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)

            if ( config_use_tracerGroup ) then
               call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupNew, 2)
               call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupCur, 1)

               if ( associated(tracersGroupNew)) then
                  !$omp do schedule(runtime) private(k)
                  do iCell = 1, nCells
                     do k = 1, maxLevelCell(iCell)
                        tracersGroupNew(:, k, iCell) = tracersGroupCur(:, k, iCell) 
                     end do
                  end do
                  !$omp end do
               end if
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
      end if
   end subroutine ocn_time_integrator_leap_timefilter!}}}

!***********************************************************************************!   
!***********************************************************************************!   

   subroutine ocn_time_integrator_leap_cleanup(block, dt, err)!{{{
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err

      integer, pointer :: nCells, nEdges, indexTemperature, indexSalinity
      integer :: iCell, iEdge, k

      type (mpas_pool_type), pointer :: statePool, meshPool, forcingPool
      type (mpas_pool_type), pointer :: diagnosticsPool, scratchPool
      type (mpas_pool_type), pointer :: tracersPool

      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessNew, normalVelocityNew
      real (kind=RKIND), dimension(:, :), pointer :: normalTransportVelocity, normalGMBolusVelocity
      real (kind=RKIND), dimension(:, :, :), pointer :: tracersGroupNew

      integer, dimension(:), pointer :: maxLevelCell

      logical, pointer :: config_use_tracerGroup
      type (mpas_pool_iterator_type) :: groupItr
      character (len=StrKIND) :: modifiedGroupName
      character (len=StrKIND) :: configName

      logical, pointer :: config_use_standardGM

      err = 0

      call mpas_pool_get_config(block % configs, 'config_use_standardGM', config_use_standardGM)

      call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)

      call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)

      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupNew, 2)
            if ( associated(tracersGroupNew) ) then
               !$omp do schedule(runtime) private(k)
               do iCell = 1, nCells
                 do k = 1, maxLevelCell(iCell)
                   tracersGroupNew(:, k, iCell) = tracersGroupNew(:, k, iCell) / layerThicknessNew(k, iCell)
                 end do
               end do
               !$omp end do
            end if
         end if
      end do

      call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 2)
      call mpas_threading_barrier()

      call ocn_vmix_implicit(dt, meshPool, diagnosticsPool, statePool, forcingPool, scratchPool, err, 2)
      call mpas_threading_barrier()

   end subroutine ocn_time_integrator_leap_cleanup!}}}

!***********************************************************************************!   
!***********************************************************************************!   

end module ocn_time_integration_leap

! vim: foldmethod=marker
