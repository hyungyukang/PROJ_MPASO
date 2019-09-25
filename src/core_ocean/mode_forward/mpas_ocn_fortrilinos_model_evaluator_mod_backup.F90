module ocn_model_evaluator_mod
!---------------------------
  use forteuchos
  use fortpetra
  use fortrilinos
! use ocn_residual_mod
!---------------------------
  use, intrinsic :: iso_c_binding
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_constants
  use mpas_dmpar

! use mpas_vector_reconstruction
! use mpas_spline_interpolation
  use mpas_timer
  use mpas_threading
  use mpas_timekeeping
  use mpas_log

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

  real(scalar_type), parameter :: zero=0., one=1.

  type(TpetraMultiVector),private,allocatable :: solnvec
  type, extends(ForModelEvaluator) :: ocnModelEvaluator
    type(TeuchosComm), private :: comm
    type(TpetraMap), private :: x_map,y_map, f_map
    type(TpetraCrsGraph), private :: graph
    type(TpetraImport), private :: importer
    type(TpetraMultiVector), private :: x
!   type(TpetraMultiVector), private :: J_diagonal
    ! -------------------------------------------------------
    ! For MPASO 
    type(domain_type) :: domain
    integer :: nCells,nEdges
    integer, dimension(:), allocatable   :: nCellsArray, nEdgesArray
    integer, dimension(:), allocatable   :: nEdgesOnEdge, nEdgesOnCell,indexToCellID
    integer, dimension(:,:), allocatable :: cellsOnEdge,edgesOnCell, edgeSignOnCell
    real (kind=RKIND) :: dt
    real (kind=RKIND), dimension(:), allocatable :: dcEdge, bottomDepth, dvEdge, areaCell
    real (kind=RKIND), dimension(:), allocatable :: sshCur,sshSubcycleCur,normalBarotropicVelocityCur
    real (kind=RKIND), dimension(:), allocatable :: barotropicForcing,barotropicCoriolisTerm
    real (kind=RKIND), dimension(:), allocatable :: CGvec_r0
    integer :: s1_nCellsArray,s1_nEdgesArray
    integer :: s1_nEdgesOnCell
    integer :: s1_cellsOnEdge, s2_cellsOnEdge, s1_edgesOnCell, s2_edgesOncell
    integer :: s1_edgeSignOnCell, s2_edgeSignOnCell
    integer :: s1_dcEdge,s1_bottomDepth,s1_dvEdge,s1_areaCell
    integer :: s1_sshCur,s1_sshSubcycleCur,s1_normalBarotropicVelocityCur
    integer :: s1_barotropicForcing, s1_barotropicCoriolisTerm,s1_CGvec_r0
    ! -------------------------------------------------------
   
  contains
    procedure, private :: create_graph
!   procedure, private :: import_mpaso_variables
    procedure :: evaluate_residual       => ocnModelEvaluator_eval_resid
    procedure :: update_solution_vector  => ocnModelEvaluator_update_solution_vector
    procedure :: get_x_map               => ocnModelEvaluator_get_x_map 
    procedure :: get_f_map               => ocnModelEvaluator_get_f_map 
    procedure :: create_operator         => ocnModelEvaluator_create_operator
    procedure :: release                 => delete_ocnModelEvaluator
  end type

  interface ocnModelEvaluator
    procedure new_ocnModelEvaluator
  end interface 

!*********************************************************************
                                contains
!*********************************************************************

  function new_ocnModelEvaluator(comm,domain,xstate,dt,n_global_vec) result(self)
    ! -------------------------------------------------------------- !
    use,intrinsic :: iso_c_binding
    type(domain_type) :: domain
    real(kind=RKIND) :: dt
    type(ocnMOdelEvaluator) :: self
    type(TeuchosComm), intent(in) :: comm 
    integer(global_size_type), intent(in) :: n_global_vec
    integer :: i,j
    integer(size_type) :: num_vecs=1
    integer(size_type) :: col
    integer :: lclrow
    real(scalar_type) :: val
    ! -------------------------------------------------------------- !
    type (block_type), pointer :: block
    type (mpas_pool_type), pointer :: statePool
    type (mpas_pool_type), pointer :: tracersPool
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: diagnosticsPool
    type (mpas_pool_type), pointer :: tendPool
    type (mpas_pool_type), pointer :: tracersTendPool

    integer :: nCells, nEdges
    integer :: iCell,iEdge,cell1,cell2
    integer, pointer :: nCellsPtr,nEdgesPtr
    integer, dimension(:), pointer :: nCellsArray, nEdgesArray
    integer, dimension(:), pointer :: nEdgesOnEdge, nEdgesOnCell, indexToCellID
    integer, dimension(:,:), pointer :: cellsOnEdge,edgesOnCell, edgeSignOnCell
    real (kind=RKIND), dimension(:), pointer :: dcEdge, bottomDepth, dvEdge, areaCell
    real (kind=RKIND), dimension(:), pointer :: sshCur,sshSubcycleCur,normalBarotropicVelocityCur
    real (kind=RKIND), dimension(:), pointer :: barotropicForcing,barotropicCoriolisTerm
    real (kind=RKIND), dimension(:), pointer :: CGvec_r0

    real (kind=RKIND), dimension(n_global_vec) :: xstate

    integer,dimension(2) :: shape2d
    ! -------------------------------------------------------------- !
   
    print*, 'PRINT in new_evaluator' 

    self = ForModelEvaluator()
    self%comm = comm

    ! owned_space
    self%x_map = TpetraMap(n_global_vec, comm)
!   self%y_map = self%x_map

!   self%importer = TpetraImport(self%x_map,self%y_map)
    self%importer = TpetraImport(self%x_map,self%x_map)
   
    ! residual space
    self%f_map = self%x_map


    self%domain = domain
    self%dt = dt

    ! Initialize the graph for W CrsMatrix object
!   self%graph = self%create_graph(self%x_map,self%y_map)
!   self%graph = self%create_graph(self%x_map,domain)
    self%graph = self%create_graph(self%x_map)

!   print*, 'After graph'
    ! Allocate space for domain and solution vectors
    self%x = TpetraMultiVector(self%x_map,num_vecs)

!   call self%xvec%doImport(self%node_coords,self%importer,TpetraINSERT)

    do i = 1,n_global_vec
      lclrow = i
      col = 1
      val = xstate(i)
      call self%x%replaceLocalValue(lclrow,col,val)
    end do

    ! -------------------------------------------------------------- !  ! For MPASO variables input
 
  end function new_ocnModelEvaluator

  !===================================================================

! function create_graph(self, x_map,y_map) result(graph)
! function create_graph(self, x_map, domain) result(graph)
  function create_graph(self, x_map) result(graph)
    ! -------------------------------------------------------------- !
    class(ocnModelEvaluator) :: self
    type(TpetraCrsGraph) :: graph
!   type(TpetraMap) :: x_map, y_map
    type(TpetraMap) :: x_map !, y_map
    integer(global_ordinal_type) :: gblrow, gblinds(1)
    integer(size_type) :: num_ent_per_row,lrow
    ! --- !
!   type(domain_type) :: domain
    type (block_type), pointer :: block
    type (mpas_pool_type), pointer :: meshPool
    integer, dimension(:), pointer :: nCellsArray, nEdgesArray
    integer, dimension(:), pointer :: nEdgesOnEdge, nEdgesOnCell, indexToCellID
    integer, dimension(:,:), pointer :: cellsOnEdge,edgesOnCell, edgeSignOnCell
    real (kind=RKIND), dimension(:), pointer :: dcEdge, bottomDepth, dvEdge, areaCell
    integer :: nvtxs,nedges,nvtxs1 
    integer,parameter :: MAX_NUM_OF_COLS=10
    integer,parameter :: MAX_LINE_LENGTH=100
    integer,dimension(:),allocatable :: rownum,rowptr,part,options,vwgt
    character(len=MAX_LINE_LENGTH) line
    integer i, io, reading_unit,j,ncol,k,objval,nparts,itmp1,itmp2,isum1
    integer :: iCell,nCells,cell1,cell2,iEdge
    double precision, dimension(MAX_NUM_OF_COLS) :: test_array
    ! -------------------------------------------------------------- !

    print*, 'PRINT in create_graph' 

    num_ent_per_row = 9
    graph = TpetraCrsGraph(x_map,x_map,num_ent_per_row,TpetraStaticProfile)

    block => self%domain % blocklist
    do while (associated(block))

      call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
      call mpas_pool_get_subpool(block % structs, 'mesh'  , meshPool)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell'  , nEdgesOnCell  )
      call mpas_pool_get_array(meshPool, 'edgesOnCell'   , edgesOnCell   )
      call mpas_pool_get_array(meshPool, 'cellsOnEdge'   , cellsOnEdge   )
 
      nCells = nCellsArray(1)

      do iCell = 1, nCells
        gblrow = x_map%getGlobalElement(iCell)
 
        do i = 1, nEdgesOnCell(iCell)

          iEdge = edgesOnCell(i,iCell)
          if ( cellsOnEdge(1,iEdge) < nCells+1 .and. cellsOnEdge(2,iEdge) < nCells+1 ) then
          gblinds(1) = cellsOnEdge(1,iEdge)
          call graph%insertGlobalIndices(gblrow,gblinds)
          gblinds(1) = cellsOnEdge(2,iEdge)
          call graph%insertGlobalIndices(gblrow,gblinds)
          endif

        end do ! i

      end do ! iCell

      block => block % next
    end do  ! block

    call graph%fillComplete()

  end function create_graph

  !===================================================================

  subroutine ocnModelEvaluator_eval_resid(self,x,f)
    
    !-----------------------------------------------------------------
    class(ocnModelEvaluator),intent(in) :: self
    class(TpetraMultiVector),intent(in) :: x
    class(TpetraMultiVector),intent(in) :: f
    integer(size_type), parameter :: ione = 1
    real(scalar_type), dimension(:), pointer :: xdata
    real(scalar_type), dimension(:), pointer :: udata
    real(scalar_type), dimension(:), pointer :: fdata
    integer(c_int)                        :: nvec
    integer(size_type) :: col
    integer :: lclrow
    real(scalar_type) :: val
    !-------
    type(domain_type) :: domain
    type(MPAS_TimeInterval_type) :: timeStep
    real(kind=RKIND) :: dt
    integer :: ierr,err_tmp
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

    integer,pointer :: nCellsPtr, nEdgesPtr
    integer, dimension(:),pointer :: nCellsArray
    integer, dimension(:),pointer :: nEdgesArray
    integer, dimension(:),pointer :: nEdgesOnCell
    integer, dimension(:,:),pointer :: cellsOnEdge
    integer, dimension(:,:),pointer :: edgesOnCell
    integer, dimension(:,:),pointer :: edgeSignOnCell
    real (kind=RKIND) :: sshTendb1,sshTendb2,sshTendAx,sshEdgeCur,thicknessSumCur,sshDiffCur
    real (kind=RKIND) :: sshEdgeMid,sshEdgeLag,thicknessSumMid,thicknessSumLag,sshDiffLag
    real (kind=RKIND) :: fluxb1,fluxb2,fluxAx,sshCurArea,sshLagArea
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

    character (len=*), parameter :: subcycleGroupName = 'subcycleFields'
    !-----------------------------------------------------------------

!   print*, 'STOP before resid' ; STOP
!   print*, 'PRINT in resid *********' 

!   udata => solnvec%getData(ione)
!   print*, udata(1:4)

    call self%update_solution_vector(x)
    call    f%putScalar(zero)

    xdata => self%x%getData(ione)
    udata => solnvec%getData(ione)

    print*, xdata(1:4)
    print*, udata(1:4)
    stop

    !**************************************************************!

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

!     sshSubcycleCur(1:nCells) = xdata(1:nCells)

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
             sshEdgeCur = 0.5_RKIND * (        sshCur(cell1) +         sshCur(cell2))
!            sshEdgeLag = 0.5_RKIND * (sshSubcycleCur(cell1) + sshSubcycleCur(cell2))
             sshEdgeLag = 0.5_RKIND * (udata(cell1) + udata(cell2))
             sshEdgeMid = 0.5_RKIND * (           sshEdgeLag + sshEdgeCur           )

             ! method 1, matches method 0 without pbcs, works with pbcs.
             thicknessSumCur = sshEdgeCur + min(bottomDepth(cell1), bottomDepth(cell2))
             thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1), bottomDepth(cell2))
             thicknessSumMid = sshEdgeMid + min(bottomDepth(cell1), bottomDepth(cell2))

             ! nabla (ssh^0)
             sshDiffCur = (        sshCur(cell2)-        sshCur(cell1)) / dcEdge(iEdge)
!            sshDiffLag = (sshSubcycleCur(cell2)-sshSubcycleCur(cell1)) / dcEdge(iEdge)
             sshDiffLag = (udata(cell2)-udata(cell1)) / dcEdge(iEdge)

             !--------------------------------------------------------------!
                fluxb1 = thicknessSumMid * normalBarotropicVelocityCur(iEdge)
                fluxb2 = thicknessSumLag * (0.5_RKIND*gravity*sshDiffCur + (-barotropicCoriolisTerm(iEdge)-barotropicForcing(iEdge)) )
                fluxAx = thicknessSumLag * sshDiffLag
   
             sshTendb1 = sshTendb1 + edgeSignOnCell(i, iCell) * fluxb1 * dvEdge(iEdge)
             sshTendb2 = sshTendb2 + edgeSignOnCell(i, iCell) * fluxb2 * dvEdge(iEdge)
             sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) * fluxAx * dvEdge(iEdge)
             !--------------------------------------------------------------!

           end do ! i
             
             sshTendb1 = (4.0_RKIND/(gravity*self%dt)) * sshTendb1
             sshTendb2 = (2.0_RKIND/(gravity   )) * sshTendb2

             sshCurArea = (4.0_RKIND/(gravity*self%dt**2.0)) *         sshCur(iCell) * areaCell(iCell)
!            sshLagArea = (4.0_RKIND/(gravity*self%dt**2.0)) * sshSubcycleCur(iCell) * areaCell(iCell)
             sshLagArea = (4.0_RKIND/(gravity*self%dt**2.0)) * udata(iCell) * areaCell(iCell)

             CGvec_r0(iCell) =-(-sshCurArea - sshTendb1 + sshTendb2)   &
                              +(-sshLagArea - sshTendAx) 

          lclrow = iCell
          col = 1
          val = CGvec_r0(iCell)
          call f%sumIntoLocalValue(lclrow,col,val)

      end do ! iCell
    

      ! -------------------------------------------------------------------------------

      block => block % next
    end do  ! block
        
!   call mpas_threading_barrier()
!   call mpas_timer_start("si halo iter r0")
!   call mpas_dmpar_exch_group_create(self%domain, subcycleGroupName)
!   call mpas_dmpar_exch_group_add_field(self%domain, subcycleGroupName, 'CGvec_r0', 1)
!   call mpas_threading_barrier()
!   call mpas_dmpar_exch_group_full_halo_exch(self%domain, subcycleGroupName)
!   call mpas_dmpar_exch_group_destroy(self%domain, subcycleGroupName)
!   call mpas_timer_stop("si halo iter r0")
!
!   print*, 'After resid calling************'
!   stop
    !**************************************************************!
    
    ! Replace value ------------------------------------------------
   
    nullify(xdata)
    nullify(udata)

  end subroutine ocnModelEvaluator_eval_resid

  !===================================================================


  subroutine ocnModelEvaluator_update_solution_vector(self,xp)
!   !----------------------------------------------------------------
    class(ocnModelEvaluator),intent(in) :: self
    class(TpetraMultiVector),intent(in) :: xp
    integer(size_type) :: num_vecs=1
!   !----------------------------------------------------------------

!   print*, 'PRINT in solution_vector'

    if (.not.allocated(solnvec)) then
      allocate(solnvec, source=TpetraMultiVector(self%x_map,num_vecs))
    endif
    call solnvec%doImport(xp, self%importer, TpetraREPLACE)
  end subroutine ocnModelEvaluator_update_solution_vector

  !==================================================================

  function ocnModelEvaluator_get_x_map(self) result(map)
    class(ocnModelEvaluator),intent(in) :: self
    type(TpetraMap) :: map

    print*, 'PRINT in get_x_map'

    map = self%x_map
  end function ocnModelEvaluator_get_x_map

  !==================================================================

  function ocnModelEvaluator_get_f_map(self) result(map)
    class(ocnModelEvaluator),intent(in) :: self
    type(TpetraMap) :: map

    print*, 'PRINT in get_f_map'

    map = self%x_map
  end function ocnModelEvaluator_get_f_map

  !==================================================================

  function ocnModelEvaluator_create_operator(self) result(op)
    class(ocnModelEvaluator),intent(in) :: self
    type(TpetraCrsMatrix) :: matrix
    type(TpetraOperator) :: op
  
    print*, 'PRINT in create_operator'

    matrix = TpetraCrsMatrix(self%graph)
    
    op = matrix_to_operator(matrix)
    
  end function ocnModelEvaluator_create_operator

  !==================================================================

  subroutine delete_ocnModelEvaluator(self)
    class(ocnModelEvaluator),intent(inout) :: self
 
    print*, 'PRINT in release'
!   print*, 'STOP in release' ; STOP

!   call self%ForModelEvaluator%release()
    call self%comm%release()
    call self%x_map%release()
    call self%f_map%release()
    call self%graph%release()
    call self%importer%release()
    call self%x%release()
    call self%importer%release()
   
    call self%ForModelEvaluator%release()
  end subroutine delete_ocnModelEvaluator

  !==================================================================

end module
