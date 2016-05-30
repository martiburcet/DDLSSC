program fem2d_poisson_driver 
  use mpi
  use mesh_class
  use data_desc_class
  use fem_assembly
  use solver 
  use timer_class
  implicit none

  type (data_desc) :: desc
  type (mesh)      :: msh
  type (matrix)    :: A 
  type (precond)   :: M 
  type (vector)    :: x, b 
  type (timer)     :: t
  integer          :: nx, ny
  real(8)          :: xl,xr,yb,yt 

  integer          :: me, np, info

  real(8)          :: err
  integer          :: it
  logical          :: conv
  integer          :: i 
  integer          :: prec_type 

  ! Initialize MPI library
  call mpi_init(info) 
  if (info /= mpi_success) then
     write(0,*) 'Error in initalizing MPI, bailing out',info 
     call mpi_abort(mpi_comm_world, -1, info)
  end if

  call read_parameters(nx,ny,xl,xr,yb,yt,prec_type)

  ! Create distributed grid descriptor
  call data_desc_create (nx, ny, mpi_comm_world, desc) 

  call mesh_create ( xl, xr, yb, yt, desc, msh )

  call matrix_create  ( desc, A )
  call vector_create  ( desc, b ) 
  call vector_create  ( desc, x ) 

  call assembly_linear_system ( msh, A, b )

  do i=1,10
    call timer_create (t, 'Precond. set-up time' )
    call timer_start (t)
    call precond_create ( A, prec_type, M ) 
    call timer_stop(t)
    call timer_report(t)

    call timer_create (t, 'Parallel PCG time' )
    call timer_start (t)
    ! Set initial solution x
    x%data = 0.0
    call pcg ( A, M, b, x, trace=1, itmax=5000, rtol=1.0e-6, atol=0.0, err=err, it=it, conv=conv)
    call timer_stop(t)
    call timer_report(t,.false.)
  
    call precond_destroy(M)
  end do

  call check_solution ( msh, x )

  ! call matrix_print ( A )

  ! Destroy (free) derived data types
  call vector_destroy (x)
  call vector_destroy (b)
  call matrix_destroy (A)
  call mesh_destroy ( msh ) 
  call data_desc_destroy (desc)

  ! Finalize MPI library
  call mpi_finalize(info)

contains

  subroutine read_parameters (nx,ny,xl,xr,yb,yt,prec_type)
    implicit none
    ! Output parameters
    integer, intent(out)  :: nx, ny
    real(8), intent(out)  :: xl, xr, yb, yt 
    integer, intent(out)  :: prec_type

    character(len=256)           :: program_name
    character(len=256)           :: argument
    integer                      :: numargs, iargc, me, ierr

    call mpi_comm_rank(mpi_comm_world, me, ierr)

    numargs = iargc()
    call getarg(0, program_name)
    if (.not. (numargs == 7) ) then
       call mpi_barrier(mpi_comm_world, ierr)
       if ( me == 0 ) then
          write (6,*) 'Usage: ', trim(program_name), ' [nx ny xl xr yb yt identity|diagonal|neumann]'
       end if
       call mpi_abort(mpi_comm_world, -1, ierr)
    end if

    call getarg(1,argument)
    read (argument,*) nx

    call getarg(2,argument)
    read (argument,*) ny

    call getarg(3,argument)
    read (argument,*) xl

    call getarg(4,argument)
    read (argument,*) xr

    call getarg(5,argument)
    read (argument,*) yb 

    call getarg(6,argument)
    read (argument,*) yt

    call getarg(7,argument)
    if ( argument == 'identity' ) then
       prec_type = identity
    else if ( argument == 'diagonal' ) then
       prec_type = diagonal
    else if ( argument == 'neumann' ) then
       prec_type = neumann 
    else
       call mpi_barrier(mpi_comm_world, ierr)
       if ( me == 0 ) then
          write (6,*) 'Usage: ', trim(program_name), ' [nx ny xl xr yb yt identity|diagonal|neumann]'
       end if
       call mpi_abort(mpi_comm_world, -1, ierr)
    end if

  end subroutine read_parameters
   
end program fem2d_poisson_driver 
