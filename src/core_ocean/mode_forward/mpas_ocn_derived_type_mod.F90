module ocn_derived_type_mod
 
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_constants
  use mpas_dmpar
  
  implicit none
   
  type,public :: ocn_derived_type
!   type (domain_type) :: domain
!   type (dm_info) :: dminfo
    integer,pointer :: nCellsPtr,nEdgesPtr
    integer,pointer,dimension(:) :: nCellsArray,nEdgesArray
   
    real (kind=RKIND) :: dt   
    real (kind=RKIND),pointer,dimension(:) :: sshCur
  end type ocn_derived_type

!*****************************************************************
                              contains 
!*****************************************************************

  subroutine ocn_imp_initialize (object,dt)
    implicit none
    type (domain_type) :: domain
    type (block_type),pointer :: block
    type (mpas_pool_type),pointer :: statePool

    type (dm_info) :: dminfo
    type (ocn_derived_type) :: object
    real (kind=RKIND) :: dt
    
    ! Dimensions
    integer :: nCells, nEdges
    integer,pointer :: nCellsPtr    
    integer,pointer,dimension(:) :: nCellsArray
    real (kind=RKIND),pointer,dimension(:) :: sshCur
   
    object%dt = dt


    block => domain % blocklist
     
    do while (associated(block))

      call mpas_pool_get_dimension(block % dimensions, 'nCells', object%nCellsPtr)
      call mpas_pool_get_array(statePool, 'ssh', object%sshCur, 1)

!     object%nCellsPtr = nCellsPtr
!     object%sshCur => sshCur


      block => block % next
    end do

  end subroutine ocn_imp_initialize

end module
