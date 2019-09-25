module implicit_mod

  use, intrinsic :: iso_c_binding
  use parallel_mod, only : haltmp, abortmp, iam

contains

  subroutine advance_imp_nonstag(elem, edge1, edge2, edge3, red, deriv, cg,   &
       hybrid, dt, pmean, tl, nets, nete, xstate)

    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev
    use physical_constants, only : g
    use element_mod, only : element_t
    use parallel_mod, only : parallel_t, abortmp
    use edge_mod, only : edgevpack, edgevunpack
    use edgetype_mod, only : EdgeBuffer_t
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk
    use time_mod, only : timelevel_t, smooth
    use control_mod, only :  topology, test_case
    use shallow_water_mod, only : tc1_velocity
    use bndry_mod, only : bndry_exchangev
    use cg_mod, only : cg_t
    use precon_mod, only : pcg_presolver_nonstag
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use derived_type_mod ,only : derived_type, initialize
    use precon_type_mod ,only : precon_type, init_precon
    use ,intrinsic :: iso_c_binding
    use global_norms_mod
    use perf_mod, only : t_startf, t_stopf
#ifdef USE_FORTRILINOS
    use fortrilinos_mod
#endif

    implicit none

    type (element_t) ,intent(inout)   :: elem(:)
    type (EdgeBuffer_t)  , intent(in) :: edge1
    type (EdgeBuffer_t)  , intent(in) :: edge2
    type (EdgeBuffer_t)  , intent(in) :: edge3
    type (reductionbuffer_ordered_1d_t)  ,intent(in) :: red
    type (derivative_t)  , intent(in) :: deriv
    type (cg_t)                       :: cg
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=real_kind), intent(in) :: dt
    real (kind=real_kind), intent(in) :: pmean
    type (timeLevel_t)   , intent(in) :: tl
    integer              , intent(in) :: nets
    integer              , intent(in) :: nete

    ! =================
    ! Local
    ! =================

    ! Thread private working set ...

    type (parallel_t)          :: par
    real (c_double) ,allocatable ,dimension(:) :: xstate

    type (element_t)  :: pc_elem(size(elem))
    type (element_t)  :: jac_elem(size(elem))

!    type (element_t)  :: pc_elem(nets-nete+1)
!    type (element_t)  :: jac_elem(nets-nete+1)

    integer    :: i,j,k,n,ie,shiftv,shiftp
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx, lx
    integer    :: comm
    integer    :: state_size

! state_object is a derived data type passed from doloca thru cpp to calc_f
! in the form of a pointer
    type(derived_type) ,target         :: state_object
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object

    !type(c_funptr)                     :: res_calc
    !type(precon_type) ,target          :: pre_object
    !type(precon_type) ,pointer         :: pptr=>NULL()

    type(derived_type) ,target          :: pre_object
    type(derived_type) ,pointer         :: pptr=>NULL()
    type(c_ptr)                         :: c_ptr_to_pre

    type(derived_type) ,target          :: ajac_object
    type(derived_type) ,pointer         :: jptr=>NULL()
    type(c_ptr)                         :: c_ptr_to_jac

    integer(c_int) :: ierr = 0

    !    type(c_funptr)                     :: precon
#ifndef USE_FORTRILINOS
    interface
       subroutine noxsolve(VectorSize, StateVector, StateData, ierr) &
            bind(C,name='noxsolve')
         use, intrinsic :: iso_c_binding
         integer(c_int)               :: VectorSize  ! length of state vector
         real(c_double), dimension(*) :: StateVector ! 1d state vector
         type(c_ptr)                  :: StateData   ! derived type
         integer(c_int)               :: ierr        ! error flag
       end subroutine noxsolve
    end interface
#endif

    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep

    comm = par%comm
! n+1 guess is result at n - need this to match up with time_mod bootstrap style
    call t_startf('implicit header')
    call t_startf('implicit_pre_noxsolve')

    shiftv = np*np*nlev*(nete-nets+1)
    shiftp = 2*np*np*nlev*(nete-nets+1)
    lx = 0
    do ie=nets,nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                xstate(lx) = elem(ie)%state%v(i,j,1,k,n0)
                xstate(lx+shiftv) = elem(ie)%state%v(i,j,2,k,n0)
                xstate(lx+shiftp) = elem(ie)%state%p(i,j,k,n0)
             end do  !np
          end do  !np
       end do  !nlev
    end do !ie

    call t_stopf('implicit_pre_noxsolve')

    call t_startf('implicit_init')

    call initialize(state_object, lenx, elem, pmean, edge1, edge2, edge3, &
         hybrid, deriv, dt, tl, nets, nete)

    fptr => state_object
    c_ptr_to_object = c_loc(fptr)

    call t_stopf('implicit_init')

    ! ForTrilinos interface to use nox and loca, and returns xstate(n+1)

    call t_startf('noxsolve')
#ifndef USE_FORTRILINOS
    call noxsolve(size(xstate), xstate, c_ptr_to_object, ierr)
#else
    call noxsolve(size(xstate), xstate, fptr, ierr)
#endif
    if (ierr /= 0) call abortmp('Error in noxsolve: Newton failed to converge')
    call t_stopf('noxsolve')
    call t_startf('implicit_post_noxsolve')

!Moved conditionals out of do loops
      lx = 0
      do ie=nets,nete
        do k=1,nlev
          do j=1,np
            do i=1,np
              lx = lx+1
              elem(ie)%state%v(i,j,1,k,np1) = xstate(lx)
              elem(ie)%state%v(i,j,2,k,np1) = xstate(lx+shiftv)
              elem(ie)%state%p(i,j,k,np1)   = xstate(lx+shiftp)
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie

     call t_stopf('implicit_post_noxsolve')

       ! ==========================================
       ! If SW test case 1, velocity is prescribed.
       ! reset v back to initial values
       ! ==========================================

       ! TODO update with vortex and swirl possibly using set_prescribed_velocity
       do ie=nets,nete
       if (topology == "cube" .and. test_case=="swtc1") then
          do k=1,nlev
             elem(ie)%state%v(:,:,:,k,np1)=tc1_velocity(elem(ie)%spherep,elem(ie)%Dinv)
             elem(ie)%state%v(:,:,:,k,n0 )=elem(ie)%state%v(:,:,:,k,np1)
             elem(ie)%state%v(:,:,:,k,nm1)=elem(ie)%state%v(:,:,:,k,np1)
          end do
       end if
       end do !ie

    call t_stopf('implicit header')

    !$OMP BARRIER

  end subroutine advance_imp_nonstag

  subroutine sw_picard_schur_op(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_schur')

    use ,intrinsic :: iso_c_binding
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nelem
    use element_mod, only : element_t
    use edge_mod, only : edgevpack, edgevunpack
    use edgetype_mod, only : EdgeBuffer_t
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
         divergence_sphere, vorticity_sphere, divergence_sphere_wk
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case,tstep_type
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    real (c_double)       :: fxtemp(2*nelemd)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1

    real (kind=real_kind), dimension(np,np,2)  :: grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term

    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv

    real (kind=real_kind) ::  dt,dti, pmean
    real (kind=real_kind) ::  gam

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: lx

    call t_startf('Precon Schur')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0     = fptr%tl%n0
    np1    = fptr%tl%np1
    nm1    = fptr%tl%nm1
    dti    = 1.0d0/fptr%dt
    dt     = fptr%dt
    ns     = fptr%nets
    ne     = fptr%nete
    pmean  = fptr%pmean
    gam    = 0.5

    if (tstep_type==13) then ! CN
       dti=2.0d0*dti ! CN has a factor of 3/2 in the time step coefficient
    endif

    if (tstep_type==12) then ! BDF2
       dti=(3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2 in the time step coefficient
    endif

    lx = 0
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                fptr%base(ie)%state%p(i,j,k,np1) = xs(lx)!p =g*h in notes
             end do  !np
          end do  !np
       end do  !nlev
    end do !ie

    do ie=ns,ne

       do j=1,np
          do i=1,np
             mv(i,j)        = fptr%base(ie)%mp(i,j)
             spheremp(i,j)  = fptr%base(ie)%spheremp(i,j)
             rspheremp(i,j) = fptr%base(ie)%rspheremp(i,j)
             ps(i,j)        = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
          end do !np
       end do !np

       do k=1,nlev

          ! scalar -> latlon vector        (grad)del phi
          grade_n3 = gradient_sphere(fptr%base(ie)%state%p(:,:,k,np1),fptr%deriv,fptr%base(ie)%Dinv)

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

          do j=1,np
             do i=1,np

                vtens(i,j,1,k,ie)= spheremp(i,j)*( &
                     grade_n3(i,j,1)) ! + d(delp)/dx

                vtens(i,j,2,k,ie)= spheremp(i,j)*( &
                     grade_n3(i,j,2))!d(delp)/dy

             end do !np
          end do !np
       end do !nlev


       ! ===================================================
       ! Pack cube edges of tendencies, rotate velocities
       ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge2, vtens(1,1,1,1,ie),2*nlev,kptr,ie)
    end do !ie

    !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge2)
    !$OMP BARRIER

    do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge2, vtens(1,1,1,1,ie), 2*nlev, kptr, ie)

    end do !ie


    !!apply rspheremp
    !
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                vtens(i,j,1,k,ie) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                vtens(i,j,2,k,ie) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
             end do
          end do
       end do
    end do !ie
    ! end of Bt cal


    ! ===========================================================
    ! Compute velocity and pressure tendencies for all levels
    ! ===========================================================
    !F_diag_inv(B'p) (B'p is in vtens and matches along element interfaces so F_diag_inv can be applied locally )

    !  ApplyFDiagInv
    lx = 0

    if (topology == "cube" .and. test_case=="swtc1") then
       do ie=ns,ne
          do k=1,nlev
             do j=1,np
                do i=1,np
                   lx = lx+1
                   fptr%base(ie)%state%v(i,j,1,k,np1) = 0.0
                   fptr%base(ie)%state%v(i,j,2,k,np1) = 0.0
                end do
             end do
          end do
       end do
    else
       do ie=ns,ne
          do k=1,nlev
             do j=1,np
                do i=1,np
                   lx = lx+1
                   fptr%base(ie)%state%v(i,j,1,k,np1) = dt*vtens(i,j,1,k,ie)
                   fptr%base(ie)%state%v(i,j,2,k,np1) = dt*vtens(i,j,2,k,ie)
                end do
             end do
          end do
       end do
    end if

    ! we now need to apply B to vtens, as well as G to dh (a.k.a. fptr%base(ie)%state%p(:,:,k,np1), and sum these two quantities
    ! now compute Gdp -BF_diag_inv(B'p)
    !
    do ie=ns,ne
       do j=1,np
          do i=1,np
             spheremp(i,j)  = fptr%base(ie)%spheremp(i,j)
             rspheremp(i,j) = fptr%base(ie)%rspheremp(i,j)
          end do !np
       end do !np

       do k=1,nlev
          !
          !         ! ==============================================
          !         ! Compute kinetic energy term at each time level
          !         ! ==============================================
          !
          do j=1,np
             do i=1,np


                ! set u
                ulatlon(i,j,1) = fptr%base(ie)%state%v(i,j,1,k,n0)   !
                ulatlon(i,j,2) = fptr%base(ie)%state%v(i,j,2,k,n0)   !

                ! set du
                up(i,j,1) = fptr%base(ie)%state%v(i,j,1,k,np1)   !
                up(i,j,2) = fptr%base(ie)%state%v(i,j,2,k,np1)   !

                dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean +ps(i,j)
                dh(i,j)   = fptr%base(ie)%state%p(i,j,k,np1)

                pv1(i,j,1) = up(i,j,1)
                pv1(i,j,2) = 0.d00
                pv2(i,j,1) = 0.d00
                pv2(i,j,2) = up(i,j,2)

                dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )
                dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )

                dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )
                dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )

             end do !np
          end do !np

          !! create operators for full residual calc
          divpv1 = divergence_sphere(pv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar
          divpv2 = divergence_sphere(pv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar
          divdpv = divergence_sphere(dpv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar
          !
          !
          !          ! ==============================================
          !          ! Compute residual terms
          !          ! ==============================================
          !
          !
          !
          do j=1,np
             do i=1,np
                !
                ptens(i,j,k,ie) = spheremp(i,j)*(&
                                            !G
                     dh(i,j)*dti+&                         !delp/dt
                     divdpv(i,j) - &
                                !-BdiagFinvB'
                     (  ( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )*&
                     (divpv1(i,j)+divpv2(i,j))))

             end do !np
          end do !np
       end do !nlev

       ! ===================================================
       ! Pack cube edges of tendencies, rotate velocities
       ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge1, ptens(1,1,1,ie),nlev,kptr,ie)
    end do !ie

    !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge1)
    !$OMP BARRIER

    do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge1, ptens(1,1,1,ie), nlev, kptr, ie)
    end do !ie

    ! ===========================================================
    ! Compute velocity and pressure tendencies for all levels
    ! ===========================================================

    lx = 0
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie)
             end do
          end do
       end do
    end do !ie

    call t_stopf('Precon Schur')

  end subroutine sw_picard_schur_op



  subroutine sw_picard_DFinvBt_op(xs, nelemd,fx,nret, c_ptr_to_object) bind(C,name='sw_picard_DFinvBt')

    use ,intrinsic :: iso_c_binding
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nelem
    use element_mod, only : element_t
    use edge_mod, only : edgevpack, edgevunpack
    use edgetype_mod, only : EdgeBuffer_t
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
         divergence_sphere, vorticity_sphere, divergence_sphere_wk
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    integer(c_int) ,intent(in) ,value  :: nret
    real (c_double) ,intent(out)       :: fx(nret)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: spheremp,rspheremp,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens

    real (kind=real_kind), dimension(np,np,2)  :: grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np)    :: zeta_n ! relative vorticity

    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv

    real (kind=real_kind) ::  dt,dti, pmean
    real (kind=real_kind) ::  gam, tot_nv

    integer    :: i,j,k,n,ie, shiftv
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: lx

    call t_startf('Precon DFinvBt')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0    = fptr%tl%n0
    np1   = fptr%tl%np1
    nm1   = fptr%tl%nm1
    dti   = 1.0d0/fptr%dt
    dt    = fptr%dt
    ns    = fptr%nets
    ne    = fptr%nete
    pmean = fptr%pmean
    gam   = 0.5

    lx = 0
    shiftv = np*np*nlev*(ne-ns+1)
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                fptr%base(ie)%state%p(i,j,k,np1) = xs(lx)!p =g*h in notes
             end do  !np
          end do  !np
       end do  !nlev
    end do !ie

    do ie=ns,ne

       do j=1,np
          do i=1,np
             mv(i,j)          = fptr%base(ie)%mp(i,j)
             spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
             rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
          end do !np
       end do !np

       do k=1,nlev

          grade_n3 = gradient_sphere(fptr%base(ie)%state%p(:,:,k,np1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector        (grad)del phi

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

          do j=1,np
             do i=1,np

                vtens(i,j,1,k,ie)= spheremp(i,j)*( &
                     grade_n3(i,j,1)) ! + d(delp)/dx

                vtens(i,j,2,k,ie)= spheremp(i,j)*( &
                     grade_n3(i,j,2))!d(delp)/dy

             end do !np
          end do !np
       end do !nlev


       ! ===================================================
       ! Pack cube edges of tendencies, rotate velocities
       ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge2, vtens(1,1,1,1,ie),2*nlev,kptr, ie)
    end do !ie

    !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge2)
    !$OMP BARRIER

    do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge2, vtens(1,1,1,1,ie), 2*nlev, kptr, ie)

    end do !ie

    ! apply rspheremp
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                vtens(i,j,1,k,ie) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                vtens(i,j,2,k,ie) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
             end do
          end do
       end do
    end do !ie
    ! end of Bt cal


    ! ===========================================================
    ! Compute velocity and pressure tendencies for all levels
    ! ===========================================================

    lx = 0

    if (topology == "cube" .and. test_case=="swtc1") then
       do ie=ns,ne
          do k=1,nlev
             do j=1,np
                do i=1,np
                   lx = lx+1
                   fx(lx) = 0.0
                   fx(lx+shiftv) = 0.0
                end do
             end do
          end do
       end do !ie
    else
       do ie=ns,ne
          do k=1,nlev
             do j=1,np
                do i=1,np
                   lx = lx+1
                   fx(lx)        = dt*vtens(i,j,1,k,ie)
                   fx(lx+shiftv) = dt*vtens(i,j,2,k,ie)
                end do
             end do
          end do
       end do !ie
    end if

    call t_stopf('Precon DFinvBt')

  end subroutine sw_picard_DFinvBt_op



  subroutine sw_picard_block_11(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_block_11')

    use ,intrinsic :: iso_c_binding
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nelem
    use element_mod, only : element_t
    use edge_mod, only : edgevpack, edgevunpack
    use edgetype_mod, only : EdgeBuffer_t
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
         divergence_sphere, vorticity_sphere, divergence_sphere_wk
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case, tstep_type
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1

    real (kind=real_kind), dimension(np,np)    :: zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity

    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2

    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, tot_nv

    integer    :: i,j,k,n,ie,shiftv
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: lx

    call t_startf('precon 11')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0    = fptr%tl%n0
    np1   = fptr%tl%np1
    nm1   = fptr%tl%nm1
    dti   = 1.0d0/fptr%dt
    ns    = fptr%nets
    ne    = fptr%nete
    pmean = fptr%pmean
    gam   = 0.5

    if (tstep_type==13) then ! CN
       dti=2.0d0*dti ! CN has a factor of 3/2 in the time step coefficient
    endif

    if (tstep_type==12) then ! BDF2
       dti=(3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2 in the time step coefficient
    endif

    lx = 0
    shiftv = np*np*nlev*(ne-ns+1)

    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
                fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx+shiftv)
             end do  !np
          end do  !np
       end do  !nlev
    end do !ie

    do ie=ns,ne

       do j=1,np
          do i=1,np
             fcor(i,j)      = fptr%base(ie)%fcor(i,j)
             mv(i,j)        = fptr%base(ie)%mp(i,j)
             spheremp(i,j)  = fptr%base(ie)%spheremp(i,j)
             rspheremp(i,j) = fptr%base(ie)%rspheremp(i,j)
          end do !np
       end do !np

       do k=1,nlev

          ! ==============================================
          ! Compute kinetic energy term at each time level
          ! ==============================================

          do j=1,np
             do i=1,np

                ! set u
                ulatlon(i,j,1) = fptr%base(ie)%state%v(i,j,1,k,n0)   !
                ulatlon(i,j,2) = fptr%base(ie)%state%v(i,j,2,k,n0)   !

                ! set du
                up(i,j,1) = fptr%base(ie)%state%v(i,j,1,k,np1)   !
                up(i,j,2) = fptr%base(ie)%state%v(i,j,2,k,np1)   !

                delv1(i,j,1) = up(i,j,1)
                delv1(i,j,2) = 0.0d0
                delv2(i,j,1) = 0.0d0
                delv2(i,j,2) = up(i,j,2)

             end do !np
          end do !np

          ! create operators for full residual calc
          zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
          newton_vortdelv2 = vorticity_sphere(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
          newton_vortdelv1 = vorticity_sphere(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

          do j=1,np
             do i=1,np

                !picard linearization of 11 block
                !note this 11 block is representing a 2x2 system for up1 and up2

                vtens(i,j,1,k,ie) = spheremp(i,j)*( &
                     up(i,j,1)*dti - &                      ! delv1/dt
                     ulatlon(i,j,2)*newton_vortdelv1(i,j)-& ! u2*d(delv1)/dy & dx
                     up(i,j,2)*(zeta_n(i,j)+fcor(i,j))-&    ! delv2*(w+f)
                     ulatlon(i,j,2)*newton_vortdelv2(i,j))  ! v2*d(delv2)/dx &dy

                vtens(i,j,2,k,ie) = spheremp(i,j)*( &
                     up(i,j,1)*(zeta_n(i,j)+fcor(i,j))+&    ! delu1*(w+f)
                     ulatlon(i,j,1)*newton_vortdelv1(i,j)+& ! u1*d(delu1)/dx & dy
                     up(i,j,2)*dti +&                       ! delu2/dt
                     ulatlon(i,j,1)*newton_vortdelv2(i,j))  ! u1*d(delu2)/dy & dx

             end do !np
          end do !np
       end do !nlev

       ! ===================================================
       ! Pack cube edges of tendencies, rotate velocities
       ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,ie)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,ie)
    end do !ie

    !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge3)
    !$OMP BARRIER

    do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge3, ptens(1,1,1,ie), nlev, kptr, ie)
       kptr=nlev
       call edgeVunpack(fptr%edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, ie)

    end do !ie
    ! ===========================================================
    ! Compute velocity and pressure tendencies for all levels
    ! ===========================================================

    lx = 0

    if (topology == "cube" .and. test_case=="swtc1") then
       do ie=ns,ne
          do k=1,nlev
             do j=1,np
                do i=1,np
                   lx = lx+1
                   fx(lx) = 0.0
                   fx(lx+shiftv) = 0.0
                end do
             end do
          end do
       end do
    else
       do ie=ns,ne
          do k=1,nlev
             do j=1,np
                do i=1,np
                   lx = lx+1
                   fx(lx)        = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                   fx(lx+shiftv) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
                end do
             end do
          end do
       end do
    end if

    call t_stopf('precon 11')

  end subroutine sw_picard_block_11



  subroutine sw_picard_block_21(xs, nelemd,fx, nret, c_ptr_to_object) bind(C,name='sw_picard_block_21')

    use ,intrinsic :: iso_c_binding
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nelem
    use element_mod, only : element_t
    use edge_mod, only : edgevpack, edgevunpack
    use edgetype_mod, only : EdgeBuffer_t
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
         divergence_sphere, vorticity_sphere, divergence_sphere_wk
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    integer(c_int) ,intent(in) ,value  :: nret
    real (c_double) ,intent(out)       :: fx(nret)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: spheremp,rspheremp,mv
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1

    real (kind=real_kind), dimension(np,np)    :: zeta_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv

    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, tot_nv

    integer    :: i,j,k,n,ie, shiftv
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: lx

    call t_startf('picard 21')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0    = fptr%tl%n0
    np1   = fptr%tl%np1
    nm1   = fptr%tl%nm1
    dti   = 1.0d0/fptr%dt
    ns    = fptr%nets
    ne    = fptr%nete
    pmean = fptr%pmean
    gam   = 0.5

    lx = 0
    shiftv = np*np*nlev*(ne-ns+1)
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
                fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx+shiftv)
             end do  !np
          end do  !np
       end do  !nlev
    end do !ie

    do ie=ns,ne

       do j=1,np
          do i=1,np
             mv(i,j)        = fptr%base(ie)%mp(i,j)
             spheremp(i,j)  = fptr%base(ie)%spheremp(i,j)
             rspheremp(i,j) = fptr%base(ie)%rspheremp(i,j)
          end do !np
       end do !np

       do k=1,nlev

          ! ==============================================
          ! Compute kinetic energy term at each time level
          ! ==============================================

          do j=1,np
             do i=1,np

                ! set u
                ulatlon(i,j,1) = fptr%base(ie)%state%v(i,j,1,k,n0)   !
                ulatlon(i,j,2) = fptr%base(ie)%state%v(i,j,2,k,n0)   !

                ! set du
                up(i,j,1) = fptr%base(ie)%state%v(i,j,1,k,np1)   !
                up(i,j,2) = fptr%base(ie)%state%v(i,j,2,k,np1)   !

                ! set dh_n=p_notes=g(h+hmean)=phi
                !          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean +ps(i,j)
                !          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)

                pv1(i,j,1) = up(i,j,1)
                pv1(i,j,2) = 0.d00
                pv2(i,j,1) = 0.d00
                pv2(i,j,2) = up(i,j,2)

             end do !np
          end do !np

          ! create operators for full residual calc
          divpv1 = divergence_sphere(pv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar
          divpv2 = divergence_sphere(pv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

          do j=1,np
             do i=1,np

                ptens(i,j,k,ie) = spheremp(i,j)*(&
                     ( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )*&
                     (divpv1(i,j)+divpv2(i,j)))

             end do !np
          end do !np
       end do !nlev


       ! ===================================================
       ! Pack cube edges of tendencies, rotate velocities
       ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,ie)
    end do !ie

    !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge3)
    !$OMP BARRIER

    do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge3, ptens(1,1,1,ie), nlev, kptr, ie)

    end do !ie
    ! ===========================================================
    ! Compute velocity and pressure tendencies for all levels
    ! ===========================================================

    lx = 0
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie)
             end do
          end do
       end do
    end do !ie

    call t_stopf('picard 21')

  end subroutine sw_picard_block_21



  subroutine sw_jacobian_op(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_jacobian')

    use ,intrinsic :: iso_c_binding
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nelem
    use element_mod, only : element_t
    use edge_mod, only : edgevpack, edgevunpack
    use edgetype_mod, only : EdgeBuffer_t
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
         divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
         vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case, tstep_type
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    !    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1

    real (kind=real_kind), dimension(np,np,2)  :: grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    !    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, tot_nv

    integer    :: i,j,k,n,ie, shiftv, shiftp
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: lx

    !call flush(6)
    !write(6,*)'sw startf'
    !call flush(6)

    call t_startf('sw jacobian op')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0    = fptr%tl%n0
    np1   = fptr%tl%np1
    nm1   = fptr%tl%nm1
    dti   = 1.0d0/fptr%dt
    ns    = fptr%nets
    ne    = fptr%nete
    pmean = fptr%pmean
    gam   = 0.5d0

    if (tstep_type==13) then ! CN
       dti=2.0d0*dti ! CN has a factor of 3/2 in the time step coefficient
    endif

    if (tstep_type==12) then ! BDF2
       dti=(3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2 in the time step coefficient
    endif

    lx = 0
    shiftv = np*np*nlev*(ne-ns+1)
    shiftp = 2*np*np*nlev*(ne-ns+1)
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
                fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx+shiftv)
                fptr%base(ie)%state%p(i,j,k,np1)   = xs(lx+shiftp)
             end do !np
          end do !np
       end do !nlev
    end do !ie

    do ie=ns,ne

       do j=1,np
          do i=1,np
             fcor(i,j)      = fptr%base(ie)%fcor(i,j)
             mv(i,j)        = fptr%base(ie)%mp(i,j)
             spheremp(i,j)  = fptr%base(ie)%spheremp(i,j)
             rspheremp(i,j) = fptr%base(ie)%rspheremp(i,j)
             ps(i,j)        = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
          end do !np
       end do !np

       do k=1,nlev

          ! ==============================================
          ! Compute kinetic energy term at each time level
          ! ==============================================

          do j=1,np
             do i=1,np

                ! set u
                ulatlon(i,j,1) = fptr%base(ie)%state%v(i,j,1,k,n0)
                ulatlon(i,j,2) = fptr%base(ie)%state%v(i,j,2,k,n0)

                ! set du
                up(i,j,1) = fptr%base(ie)%state%v(i,j,1,k,np1)
                up(i,j,2) = fptr%base(ie)%state%v(i,j,2,k,np1)

                delv1(i,j,1) = up(i,j,1)
                delv1(i,j,2) = 0.0d0
                delv2(i,j,1) = 0.0d0
                delv2(i,j,2) = up(i,j,2)

                !test
                !         up(i,j,1)=v1p
                !          up(i,j,2)=v2p
                !end test

                ! where is ps???? shouldn't ps be included in dp? ps=g*h_s in comment above,
                ! ANSWER: NO, also no pmean, otherwise we would be adding these contributions at each nonlinear update
                !what is p
                !what is pmean
                !ps= surface height g*h_s in notes

                ! set dh_n=p_notes=g(h+hmean)=phi
                dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean +ps(i,j)
                !test unit vector
                !          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)

                ! set dp
                !dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)+ps(i,j)+pmean
                !dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)+pmean
                dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)

                pv1(i,j,1) = up(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
                pv1(i,j,2) = 0.d00

                pv2(i,j,1) = 0.d00
                pv2(i,j,2) = up(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )

                dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )
                dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )

                !dpv(i,j,1) = ulatlon(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,np1))
                !dpv(i,j,2) = ulatlon(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,np1))

             end do !np
          end do !np

          ! create operators for full residual calc
          grade_n1 = gradient_sphere(up(:,:,1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)delu_1
          grade_n2 = gradient_sphere(up(:,:,2),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)delu_2
          grade_n3 = gradient_sphere(dh,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector        (grad)del phi

          newton_grade_n1 = gradient_sphere(ulatlon(:,:,1)*up(:,:,1) ,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)u_1
          newton_grade_n2 = gradient_sphere(ulatlon(:,:,2)*up(:,:,2),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)u_2
          newton_grade_n3 = gradient_sphere(dh_n,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector        (grad)phi

          zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
          newton_vortdelv2 = vorticity_sphere(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
          newton_vortdelv1 = vorticity_sphere(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w

          divpv1 = divergence_sphere(pv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar
          divpv2 = divergence_sphere(pv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar
          divdpv = divergence_sphere(dpv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

          do j=1,np
             do i=1,np

                ! Jacobian
                vtens(i,j,1,k,ie) = spheremp(i,j)*( &
                     up(i,j,1)*dti - &                      ! delv1/dt
                     ulatlon(i,j,2)*newton_vortdelv1(i,j)+& ! u2*d(delv1)/dy & dx
                     newton_grade_n1(i,j,1)- &              ! delv1*d(v1)/dx
                     up(i,j,2)*(zeta_n(i,j)+fcor(i,j))-&    ! delv2*(w+f)
                     ulatlon(i,j,2)*newton_vortdelv2(i,j)+& ! v2*d(delv2)/dx &dy
                     newton_grade_n2(i,j,1)+ &              ! -delv2*d(v2)/dx
                     grade_n3(i,j,1))                       ! + d(delp)/dx

                vtens(i,j,2,k,ie) = spheremp(i,j)*( &
                     up(i,j,1)*(zeta_n(i,j)+fcor(i,j))+&    ! delu1*(w+f)
                     ulatlon(i,j,1)*newton_vortdelv1(i,j)+& ! u1*d(delu1)/dx & dy
                     newton_grade_n2(i,j,2)+&               ! delu1*d(u2)/dy
                     up(i,j,2)*dti +&                       ! delu2/dt
                     ulatlon(i,j,1)*newton_vortdelv2(i,j)+& ! u1*d(delu2)/dy & dx
                     newton_grade_n1(i,j,2)+ &              ! delu1*d(u1)/dy
                     grade_n3(i,j,2))                       ! d(delp)/dy

                ptens(i,j,k,ie) = spheremp(i,j)*(&
                     dh(i,j)*dti+&                          ! delp/dt
                     divpv1(i,j)+divpv2(i,j)+divdpv(i,j))

             end do !np
          end do !np
       end do !nlev

       ! ===================================================
       ! Pack cube edges of tendencies, rotate velocities
       ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,ie)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,ie)
    end do !ie

    !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge3)
    !$OMP BARRIER

    do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge3, ptens(1,1,1,ie), nlev, kptr, ie)
       kptr=nlev
       call edgeVunpack(fptr%edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, ie)

    end do !ie
    ! ===========================================================
    ! Compute velocity and pressure tendencies for all levels
    ! ===========================================================

    lx = 0
    if (topology == "cube" .and. test_case=="swtc1") then
       do ie=ns,ne
          do k=1,nlev
             do j=1,np
                do i=1,np
                   lx = lx+1
                   fx(lx) = 0.0
                   fx(lx+shiftv) = 0.0
                end do
             end do
          end do
       end do
    else
       do ie=ns,ne
          do k=1,nlev
             do j=1,np
                do i=1,np
                   lx = lx+1
                   fx(lx)        = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                   fx(lx+shiftv) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
                   fx(lx+shiftp) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie)
                end do
             end do
          end do
       end do
    end if

    call t_stopf('sw jacobian op')

  end subroutine sw_jacobian_op



  subroutine update_jac_state(xs, nelemd, c_ptr_to_object) bind(C,name='update_jac_state')

    use ,intrinsic :: iso_c_binding
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nelem
    use element_mod, only : element_t
    use derived_type_mod, only : derived_type
    use perf_mod, only : t_startf, t_stopf

    implicit none

    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer                            :: ns
    integer                            :: ne

    integer :: i, j, k, n, ie, shiftv, shiftp
    integer :: nm1, n0, np1
    integer :: lx

    ! write(6,*)'update jac state\n'
    ! call t_startf('update jac')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0  = fptr%tl%n0
    np1 = fptr%tl%np1
    nm1 = fptr%tl%nm1
    ns  = fptr%nets
    ne  = fptr%nete

    lx = 0
    shiftv = np*np*nlev*(ne-ns+1)
    shiftp = 2*np*np*nlev*(ne-ns+1)
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                fptr%base(ie)%state%v(i,j,1,k,n0) = xs(lx)
                fptr%base(ie)%state%v(i,j,2,k,n0) = xs(lx+shiftv)
                fptr%base(ie)%state%p(i,j,k,n0)   = xs(lx+shiftp)
             end do ! np
          end do ! np
       end do ! nlev
    end do !ie

    ! call t_stopf('update jac')

  end subroutine update_jac_state



  subroutine update_prec_state(xs, nelemd, c_ptr_to_object) bind(C,name='update_prec_state')

    use ,intrinsic :: iso_c_binding
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nelem
    use element_mod, only : element_t
    use derived_type_mod, only : derived_type
    use perf_mod, only : t_startf, t_stopf

    implicit none

    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    integer    :: i,j,k,n,ie, shiftv, shiftp
    integer    :: nm1,n0,np1
    integer    :: lx

    !write(6,*)'update prec state\n'
    !    call t_startf('residual lin')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0  = fptr%tl%n0
    np1 = fptr%tl%np1
    nm1 = fptr%tl%nm1
    ns  = fptr%nets
    ne  = fptr%nete

    lx = 0
    shiftv = np*np*nlev*(ne-ns+1)
    shiftp = 2*np*np*nlev*(ne-ns+1)
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                fptr%base(ie)%state%v(i,j,1,k,n0) = xs(lx)
                fptr%base(ie)%state%v(i,j,2,k,n0) = xs(lx+shiftv)
                fptr%base(ie)%state%p(i,j,k,n0)   = xs(lx+shiftp)
             end do !np
          end do !np
       end do !nlev
    end do !ie

  end subroutine update_prec_state



  subroutine get_jac_vector(xs, nelemd, c_ptr_to_object) bind(C,name='get_jac_vector')

    use ,intrinsic :: iso_c_binding
    use kinds, only : real_kind
    use dimensions_mod, only :  np, nlev, nelem
    use element_mod, only : element_t
    use derived_type_mod, only : derived_type
    use perf_mod, only : t_startf, t_stopf

    implicit none

    real (c_double)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    integer    :: i,j,k,n,ie,shiftv,shiftp
    integer    :: nm1,n0,np1
    integer    :: lx

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0 = fptr%tl%n0
    ns = fptr%nets
    ne = fptr%nete

    lx = 0
    shiftv = np*np*nlev*(ne-ns+1)
    shiftp = 2*np*np*nlev*(ne-ns+1)
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                xs(lx)        = fptr%base(ie)%state%v(i,j,1,k,n0)
                xs(lx+shiftv) = fptr%base(ie)%state%v(i,j,2,k,n0)
                xs(lx+shiftp) = fptr%base(ie)%state%p(i,j,k,n0)
             end do !np
          end do !np
       end do !nlev
    end do !ie

  end subroutine get_jac_vector



  subroutine test_id(xs, nelemd, fx, c_ptr_to_object) bind(C,name='test_id')
    use ,intrinsic :: iso_c_binding

    implicit none

    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(in)        :: xs(nelemd)
    real (c_double) ,intent(out)       :: fx(nelemd)
    type(c_ptr)                        :: c_ptr_to_object

    fx=xs

  end subroutine test_id



  subroutine get_discrete_params(cnets,cnete,cnp,cnlev, c_ptr_to_object) bind(C,name='get_discrete_params')

    use ,intrinsic :: iso_c_binding
    use kinds, only : real_kind
    use dimensions_mod, only :  np, nlev, nelem
    use element_mod, only : element_t
    use derived_type_mod, only : derived_type
    use perf_mod, only : t_startf, t_stopf

    implicit none

    integer (c_int)        :: cnets
    integer (c_int)        :: cnete
    integer (c_int)        :: cnp
    integer (c_int)        :: cnlev
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    cnets         = fptr%nets
    cnete         = fptr%nete

    cnp=np
    cnlev=nlev

  end subroutine get_discrete_params

  ! =============================================================================
  ! subroutines to stop run
  ! =============================================================================

  subroutine haltrun() bind(C,name='haltrun')
    ! ---------------------------------------------------------------------------
    ! wrapper to MPI finalize
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    implicit none

    call haltmp("halt called from c++")

  end subroutine haltrun


  subroutine abortrun() bind(C,name='abortrun')
    ! ---------------------------------------------------------------------------
    ! wrapper to MPI about
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    implicit none

    call abortmp("abort called from c++")

  end subroutine abortrun

end module implicit_mod
