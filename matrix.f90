module matrix_class
  use data_desc_class
  implicit none

  type matrix
     integer                  :: m         ! Number of rows in the matrix
     integer                  :: n         ! Number of cols in the matrix
     real(8), allocatable     :: data(:,:) ! Matrix entries
     type(data_desc), pointer :: desc      ! Pointer to the FE mesh descriptor
  end type matrix

contains
  
  subroutine matrix_create ( desc, m ) 
    implicit none
    ! Parameters
    type(data_desc), intent(in), target :: desc
    type (matrix), intent(out)            :: m

    ! Locals
    integer :: k, j, i 

    ! Point to the descriptor of the distributed FE mesh
    m%desc => desc 

    m%m = m%desc%node_num
    m%n = m%desc%node_num
    allocate ( m%data(1:m%m,1:m%n) )
    m%data = 0.0
    
  end subroutine matrix_create

  subroutine matrix_destroy ( m ) 
    implicit none
    ! Parameters
    type (matrix), intent(inout)     :: m
   
    m%m = -1
    m%n = -1
    deallocate ( m%data )

  end subroutine matrix_destroy

  subroutine matrix_print ( m ) 
    implicit none
    ! Parameters
    type (matrix), intent(in) :: m
    
    ! Locals 
    integer :: i 

    do i=1,m%n
       write(*,'(20(f10.2,1x))') m%data(i,:)
    end do
  end subroutine matrix_print

end module matrix_class
