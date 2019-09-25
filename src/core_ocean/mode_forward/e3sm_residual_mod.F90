module residual_mod

  use, intrinsic :: iso_c_binding

contains

  subroutine residual_c(xstate, fx, nelemd, c_ptr_to_object) bind(C,name='calc_f')
    use ,intrinsic :: iso_c_binding
    use derived_type_mod, only : derived_type
    implicit none
    real (c_double) ,intent(in)        :: xstate(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr
    call residual(xstate, fx, nelemd, fptr)
  end subroutine residual_c

  subroutine residual(xstate, fx, nelemd, fptr)
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nelem
    use element_mod, only : element_t
    use edge_mod, only : edgevpack, edgevunpack
    use edgetype_mod, only: EdgeBuffer_t
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case, tstep_type
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : derived_type
    use perf_mod, only : t_startf, t_stopf

    implicit none

    real (c_double) ,intent(in)        :: xstate(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    type(derived_type), pointer        :: fptr
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,ps,px
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,E_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um

    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam

    integer    :: i,j,k,n,ie,shiftv,shiftp
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lx

    call t_startf('FE implicit')
    call t_startf('FE_implicit_init')

    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    ns         = fptr%nets
    ne         = fptr%nete
    pmean      = fptr%pmean
    gam        = 0.5

!Moved conditionals out of do loops
    shiftv = np*np*nlev*(ne-ns+1)
    shiftp = 2*np*np*nlev*(ne-ns+1)

       call t_stopf('FE_implicit_init')
       call t_startf('FE_implicit_KE_resid_calc')
    lx = 0
    do ie=ns,ne
        do j=1,np
            do i=1,np
                fcor(i,j)        = fptr%base(ie)%fcor(i,j)
                spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
                ps(i,j)          = fptr%base(ie)%state%ps(i,j)
            end do
        end do
        do k=1,nlev

            ! ==============================================
            ! Compute kinetic energy term at each time level
            ! ==============================================

            do j=1,np
                do i=1,np

                    ! set u
                    ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   !
                    ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !

                    ! set du
                    lx = lx+1
                    up(i,j,1)      = xstate(lx)   !
                    up(i,j,2)      = xstate(lx+shiftv)   !
                    px(i,j)           = xstate(lx+shiftp)
                    ! set un-1
                    um(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,nm1)   !
                    um(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,nm1)   !


                    E(i,j)   = 0.5D0*(up(i,j,1)**2 + up(i,j,2)**2)  + &
                    px(i,j) + ps(i,j)
                    if (tstep_type==13) then  !compute extra terms needed for CN
                        E_n(i,j) = 0.5D0*(ulatlon(i,j,1)**2 + ulatlon(i,j,2)**2)  + &
                        fptr%base(ie)%state%p(i,j,k,n0) + ps(i,j)
                        pv_n(i,j,1) = ulatlon(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
                        pv_n(i,j,2) = ulatlon(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
                    endif

                    pv(i,j,1) = up(i,j,1)*( fptr%pmean + px(i,j) )
                    pv(i,j,2) = up(i,j,2)*( fptr%pmean + px(i,j) )

                end do
            end do

            !call flush(6)
            !stop
            grade = gradient_sphere(E,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector
            zeta = vorticity_sphere(up,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar
            div = divergence_sphere(pv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar
            ! residual time level n
            if (tstep_type==13) then
                !compute extra terms needed for CN
                grade_n = gradient_sphere(E_n,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector
                zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar
                div_n = divergence_sphere(pv_n,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar
                do j=1,np
                    do i=1,np
                       vtens_n(i,j,1,k,ie)=spheremp(i,j)* &
                       (ulatlon(i,j,2)*(fcor(i,j) + zeta_n(i,j)) - grade_n(i,j,1))
                       vtens_n(i,j,2,k,ie)=spheremp(i,j)* &
                       (-ulatlon(i,j,1)*(fcor(i,j) + zeta_n(i,j)) - grade_n(i,j,2))
                       ptens_n(i,j,k,ie) =  -spheremp(i,j)*div_n(i,j)
                       vtens(i,j,1,k,ie)=spheremp(i,j)* &
                       (up(i,j,2)*(fcor(i,j) + zeta(i,j)) - grade(i,j,1))
                       vtens(i,j,2,k,ie)=spheremp(i,j)* &
                       (-up(i,j,1)*(fcor(i,j) + zeta(i,j)) - grade(i,j,2))

                       ptens(i,j,k,ie) =  -spheremp(i,j)*div(i,j)

                       ptens(i,j,k,ie) = spheremp(i,j)*(px(i,j)- &
                       fptr%base(ie)%state%p(i,j,k,n0))*dti - 0.5*ptens(i,j,k,ie) &
                       - 0.5*ptens_n(i,j,k,ie)

                       vtens(i,j,1,k,ie) = spheremp(i,j)* &
                       (up(i,j,1)-ulatlon(i,j,1))*dti - 0.5*vtens(i,j,1,k,ie) &
                       - 0.5*vtens_n(i,j,1,k,ie)

                       vtens(i,j,2,k,ie) = spheremp(i,j)* &
                       (up(i,j,2)-ulatlon(i,j,2))*dti - 0.5*vtens(i,j,2,k,ie) &
                       - 0.5*vtens_n(i,j,2,k,ie)

                    end do
                end do
            else if (tstep_type==12) then !BDF2 2nd order
                do j=1,np
                    do i=1,np
                        vtens(i,j,1,k,ie)=spheremp(i,j)* &
                        (up(i,j,2)*(fcor(i,j) + zeta(i,j)) - grade(i,j,1))
                        vtens(i,j,2,k,ie)=spheremp(i,j)* &
                        (-up(i,j,1)*(fcor(i,j) + zeta(i,j)) - grade(i,j,2))

                        ptens(i,j,k,ie) =  -spheremp(i,j)*div(i,j)
                       if (nstep==0) then ! BE bootstrap
                           ptens(i,j,k,ie)   = spheremp(i,j)*(px(i,j)- &
                           fptr%base(ie)%state%p(i,j,k,n0))*dti - ptens(i,j,k,ie)

                           vtens(i,j,1,k,ie) = spheremp(i,j)* &
                           (up(i,j,1)-ulatlon(i,j,1))*dti - vtens(i,j,1,k,ie)

                           vtens(i,j,2,k,ie) = spheremp(i,j)* &
                           (up(i,j,2)-ulatlon(i,j,2))*dti - vtens(i,j,2,k,ie)

                           !      if (nstep==0) then ! CN bootstrap
                           !       ptens(i,j,k,ie) = spheremp(i,j)*(px(i,j)- &
                           !          fptr%base(ie)%state%p(i,j,k,n0))*dti - 0.5*ptens(i,j,k,ie) &
                           !          - 0.5*ptens_n(i,j,k,ie)
                           !
                           !       vtens(i,j,1,k,ie) = spheremp(i,j)* &
                           !         (up(i,j,1)-ulatlon(i,j,1))*dti - 0.5*vtens(i,j,1,k,ie) &
                           !          - 0.5*vtens_n(i,j,1,k,ie)
                           !
                           !       vtens(i,j,2,k,ie) = spheremp(i,j)* &
                           !         (up(i,j,2)-ulatlon(i,j,2))*dti - 0.5*vtens(i,j,2,k,ie) &
                           !          - 0.5*vtens_n(i,j,2,k,ie)

                       else

                           ptens(i,j,k,ie) = &
                           (1+gam)*spheremp(i,j)*(px(i,j) - &
                           fptr%base(ie)%state%p(i,j,k,n0))*dti - &
                           gam*spheremp(i,j)*(fptr%base(ie)%state%p(i,j,k,n0) - &
                           fptr%base(ie)%state%p(i,j,k,nm1))*dti - &
                           ptens(i,j,k,ie)

                           vtens(i,j,1,k,ie) = &
                           (1+gam)*spheremp(i,j)*(up(i,j,1)-ulatlon(i,j,1))*dti - &
                           gam*spheremp(i,j)*(ulatlon(i,j,1)-um(i,j,1))*dti - &
                           vtens(i,j,1,k,ie)

                           vtens(i,j,2,k,ie) = &
                           (1+gam)*spheremp(i,j)*(up(i,j,2)-ulatlon(i,j,2))*dti - &
                           gam*spheremp(i,j)*(ulatlon(i,j,2)-um(i,j,2))*dti - &
                           vtens(i,j,2,k,ie)
                       end if
                    end do
                end do

            else ! Backward Euler 1st order method
                do j=1,np
                    do i=1,np
                        vtens(i,j,1,k,ie)=spheremp(i,j)* &
                        (up(i,j,2)*(fcor(i,j) + zeta(i,j)) - grade(i,j,1))
                        vtens(i,j,2,k,ie)=spheremp(i,j)* &
                        (-up(i,j,1)*(fcor(i,j) + zeta(i,j)) - grade(i,j,2))

                        ptens(i,j,k,ie) =  -spheremp(i,j)*div(i,j)

                       ptens(i,j,k,ie)   = spheremp(i,j)*(px(i,j)- &
                       fptr%base(ie)%state%p(i,j,k,n0))*dti - ptens(i,j,k,ie)

                       vtens(i,j,1,k,ie) = spheremp(i,j)* &
                       (up(i,j,1)-ulatlon(i,j,1))*dti - vtens(i,j,1,k,ie)

                       vtens(i,j,2,k,ie) = spheremp(i,j)* &
                       (up(i,j,2)-ulatlon(i,j,2))*dti - vtens(i,j,2,k,ie)
                    end do
                end do
            endif
        end do

        ! ===================================================
        ! Pack cube edges of tendencies, rotate velocities
        ! ===================================================
        kptr=0
        call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,ie)
        kptr=nlev
        call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,ie)
    end do

   !pw++
   call t_stopf('FE_implicit_KE_resid_calc')
   !pw--

   !$OMP BARRIER
   !pw++
   call t_startf('FE_implicit_bndry_ex')
   !pw--
   !$OMP BARRIER
   call bndry_exchangeV(fptr%hybrid,fptr%edge3)
   !$OMP BARRIER
   !pw++
   call t_stopf('FE_implicit_bndry_ex')
   !pw--
   !$OMP BARRIER

   !pw++
   call t_startf('FE_implicit_bndry_unpack')
   !pw--


   do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge3, ptens(1,1,1,ie), nlev, kptr, ie)

       kptr=nlev
       call edgeVunpack(fptr%edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, ie)

   end do !ie

   call t_stopf('FE_implicit_bndry_unpack')

   ! ===========================================================
   ! Compute velocity and pressure tendencies for all levels
   ! ===========================================================

   call t_startf('FE_implicit_vel_pres')


   !Moved conditionals out of do loops
   lx = 0
   if (topology == "cube" .and. test_case=="swtc1") then
       do ie=ns,ne
           do k=1,nlev
               do j=1,np
                   do i=1,np
                       lx = lx+1
                       fx(lx) = 0.0
                       fx(lx+shiftv) = 0.0
                       fx(lx+shiftp) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie)
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
                       fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                       fx(lx+shiftv)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
                       fx(lx+shiftp) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie)
                   end do
               end do
           end do
       end do !ie
   end if

     call t_stopf('FE_implicit_vel_pres')
     call t_stopf('FE implicit')

  end subroutine residual

end module residual_mod
