module data_desc_class
  use mpi
  implicit none

  ! Derived data type which describes the
  ! data distribution of a uniform 2D FE mesh
  ! of P1 elements (triangles)
  type data_desc
      integer  :: nx_g, ny_g 
      integer  :: nx_l, ny_l
      integer  :: node_num
      integer  :: elem_num

      integer  :: me, np 

      integer  :: up_pid
      integer  :: up_off   ! Node identifier of the first node shared by this process
                           ! and the one inmediately on top of it

      integer  :: down_pid     
      integer  :: down_off ! Node identifier of the first node shared by this process
                           ! and the one inmediately below it
      
      integer  :: mpi_comm
  end type data_desc

contains

  subroutine data_desc_create (nx,ny,comm,desc)
    implicit none
    ! Parameters
    integer, intent(in)          :: nx, ny, comm
    type(data_desc), intent(out) :: desc

    ! Locals
    integer :: info
    integer :: i
    
    desc%nx_g = nx
    desc%ny_g = ny

    ! Duplicate communicator in order to avoid conflict
    ! among multiple instances of type(data_desc) within
    ! the same parallel program. This will also avoid 
    ! conflicts among this code and MPI-parallel third party 
    ! libraries in case of the mpi_comm_world communicator.
    call mpi_comm_dup(comm,desc%mpi_comm,info)

    ! How many processes are in desc%mpi_comm ?
    call mpi_comm_size(desc%mpi_comm,desc%np,info) 

    ! Which process identifier am I ?
    call mpi_comm_rank(desc%mpi_comm,desc%me,info) 

    ! Compute up and down process identifiers
    desc%down_pid = mod ( desc%me - 1, desc%np )
    if ( desc%down_pid == -1 ) then
       desc%down_pid = desc%np -1
    end if
    desc%up_pid   = mod ( desc%me + 1, desc%np ) 

    ! Compute number of local nodes in each direction
    ! The topmost process (i.e., desc%me == np -1) gets 
    ! (potentially) additional rows of nodes of the finite
    ! element mesh. In particular, those remaining from the
    ! integer division of ny_g and the number of processes 
    desc%nx_l = desc%nx_g
    desc%ny_l = (desc%ny_g-1)/desc%np + 1
    if ( desc%me == desc%np -1 ) then
       desc%ny_l = desc%ny_l + mod ( desc%ny_g-1, desc%np )
    end if

    desc%node_num = desc%nx_l * desc%ny_l
    desc%elem_num =  2 * ( desc%nx_l - 1 ) * ( desc%ny_l  - 1 )

    ! Compute local node identifiers which are shared among
    ! me and desc%down_pid and me and desc%up_pid

    ! Here is the disposition of local node identifiers
    ! for a local 3x10 portion:
    !
    !   21---22---23---24----------------------------30
    !    |\ 8 |\10 |\12 |
    !    | \  | \  | \  |
    !    |  \ |  \ |  \ |  \ |
    !    |  7\|  9\| 11\|   \|
    !   11---12---13---14---15---16---17---18---19---20
    !    |\ 2 |\ 4 |\ 6 |\  8|                   |\ 18|
    !    | \  | \  | \  | \  |                   | \  |
    !    |  \ |  \ |  \ |  \ |      ...          |  \ |
    !    |  1\|  3\|  5\| 7 \|                   |17 \|
    !    1----2----3----4----5----6----7----8----9---10
    
    desc%down_off = 1 
    desc%up_off   = desc%nx_l*(desc%ny_l-1)+1

  end subroutine data_desc_create
  
  subroutine data_desc_destroy (desc)
    implicit none
    ! Parameters
    type(data_desc), intent(inout) :: desc

    desc%nx_g = -1 
    desc%ny_g = -1
    desc%nx_l = -1 
    desc%ny_l = -1 

    desc%node_num = -1
    desc%elem_num = -1
    
    desc%me = -1
    desc%np = -1 
    
    desc%up_pid = -1 
    desc%up_off = -1 
    
    desc%down_pid = -1
    desc%down_off   = -1 
    
    desc%mpi_comm = -1

  end subroutine data_desc_destroy

end module data_desc_class
