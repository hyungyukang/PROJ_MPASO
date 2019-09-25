module model_evaluator_mod
  use forteuchos
  use fortpetra
  use fortrilinos
  use derived_type_mod, only : derived_type
  use residual_mod, only : residual

  use, intrinsic :: ISO_C_BINDING
  implicit none

  type, extends(ForModelEvaluator) :: E3SMModelEvaluator
    type(TpetraMap),    private :: x_map
    type(derived_type), pointer, private :: state_data

  contains
    procedure :: evaluate_residual      => E3SMModelEvaluator_eval_resid
    procedure :: get_x_map              => E3SMModelEvaluator_get_x_map
    procedure :: get_f_map              => E3SMModelEvaluator_get_f_map
    procedure :: release                => delete_E3SMModelEvaluator

    procedure :: update_state_data      => E3SMModelEvaluator_update_state_data
    ! procedure :: update_solution_vector => E3SMModelEvaluator_update_solution_vector
  end type

  interface E3SMModelEvaluator
    procedure new_E3SMModelEvaluator
  end interface

contains

  function new_E3SMModelEvaluator(map) &
      result(self)
    use ,intrinsic :: iso_c_binding
    type(E3SMModelEvaluator)    :: self
    type(TpetraMap), intent(in) :: map

    self = ForModelEvaluator()

    ! owned space
    self%x_map = map

  end function new_E3SMModelEvaluator

  subroutine E3SMModelEvaluator_eval_resid(self, x, f)
    class(E3SMModelEvaluator), intent(in)   :: self
    class(TpetraMultiVector), intent(in)    :: x
    class(TpetraMultiVector), intent(inout) :: f

    integer(size_type), parameter :: ione = 1

    real(c_double), dimension(:), pointer :: xdata
    real(c_double), dimension(:), pointer :: fdata
    integer(c_int)                        :: nelemd

    nelemd = x%getLocalLength()

    xdata => x%getData(ione)
    fdata => f%getData(ione)

    call residual(xdata, fdata, nelemd, self%state_data)
  end subroutine E3SMModelEvaluator_eval_resid

  subroutine E3SMModelEvaluator_update_state_data(self, data)
    class(E3SMModelEvaluator) :: self
    type(derived_type), target, intent(in) :: data

    self%state_data => data
  end subroutine E3SMModelEvaluator_update_state_data

  function E3SMModelEvaluator_get_x_map(self) result(map)
    class(E3SMModelEvaluator), intent(in) :: self
    type(TpetraMap) :: map

    map = self%x_map
  end function E3SMModelEvaluator_get_x_map

  function E3SMModelEvaluator_get_f_map(self) result(map)
    class(E3SMModelEvaluator), intent(in) :: self
    type(TpetraMap) :: map

    map = self%x_map
  end function E3SMModelEvaluator_get_f_map

  subroutine delete_E3SMModelEvaluator(self)
    class(E3SMModelEvaluator), intent(inout) :: self

    ! TODO
    ! Release all member data

    ! Call base class release()
    call self%ForModelEvaluator%release()
  end subroutine


end module
