module timer_class
  use mpi
  implicit none

  integer, parameter :: max_message = 256

  type timer
      integer                :: mpi_comm     ! MPI Communicator
      integer                :: root     ! PID responsible of gathering/reporting timings
      character(max_message) :: message  ! Concept being measured (e.g., assembly)
      real(8)                :: start    ! last call to start
      real(8)                :: stop     ! last call to stop
      real(8)                :: accum    ! sum of all stop-start 
  end type timer

contains
  
    subroutine timer_create ( t, message, mpi_comm, root )
#ifdef MPI_MOD
      use mpi
#endif
      implicit none 
#ifdef MPI_H
      include 'mpif.h'
#endif
      ! Parameters
      type(timer)    , intent(out)          :: t
      character*(*)  , intent(in)           :: message
      integer        , intent(in), optional :: mpi_comm
      integer        , intent(in), optional :: root
      
      ! Locals
      integer                               :: mpi_comm_, root_

      mpi_comm_ = mpi_comm_world
      root_  = 0 
      if ( present(mpi_comm) ) mpi_comm_ = mpi_comm
      if ( present(root)  ) root_  = root 
      
      t%mpi_comm   = mpi_comm_
      t%root    = root_ 
      t%message = trim(message)
      t%start   = 0.0 
      t%stop    = 0.0
      t%accum   = 0.0
      ! t%accum   = 1.79769E+308 ! Largest double precision number
    end subroutine timer_create

    subroutine timer_init ( t )
      implicit none 
      ! Parameters
      type(timer), intent(inout)          :: t
      t%start  = 0.0 
      t%stop   = 0.0
      t%accum  = 0.0
      t%accum  = 1.79769E+308 ! Largest double precision number
    end subroutine timer_init

    subroutine timer_start ( t )
      implicit none 
      ! Parameters
      type(timer), intent(inout)          :: t
      ! Locals
      integer :: ierr
      call mpi_barrier (t%mpi_comm, ierr)
      t%start  = mpi_wtime()
    end subroutine timer_start

    subroutine timer_stop ( t )
      implicit none 
      ! Parameters
      type(timer), intent(inout)          :: t
      t%stop = mpi_wtime()
      
      if ( t%stop - t%start >= 0.0) then
         t%accum = t%accum + (t%stop - t%start)
      end if
      t%start  = 0.0 
      t%stop   = 0.0

    end subroutine timer_stop
   
    subroutine timer_report ( t, header )
      implicit none 
      ! Parameters
      type(timer), intent(inout) :: t
      logical, intent(in), optional  :: header 
      
      ! Locals
      character(len=*), parameter    :: fmt_header = '(a25,1x,3(2x,a15),3(2x,a15))'
      character(len=*), parameter    :: fmt_data   = '(a25,1x,3(2x,es15.9),3(2x,es15.9))'
      real(8)                        :: accum_max, accum_min, accum_sum
      integer                        :: my_id, num_procs, ierr
      logical                        :: header_

      call mpi_comm_size ( t%mpi_comm, num_procs, ierr )
      call mpi_comm_rank ( t%mpi_comm, my_id    , ierr )

      accum_max = t%accum
      accum_min = t%accum
      accum_sum = t%accum
      header_ = .true.
      if (present(header)) header_ = header 

      if ( header_ ) then 
        if ( my_id == t%root ) write(*,fmt_header) '', 'Min (secs.)', 'Max (secs.)', 'Avg (secs.)'
      end if         
  
      call mpi_allreduce ( accum_min, accum_min, 1, mpi_real8, mpi_max, t%mpi_comm, ierr  )
      call mpi_allreduce ( accum_max, accum_max, 1, mpi_real8, mpi_min, t%mpi_comm, ierr  )
      call mpi_allreduce ( accum_sum, accum_sum, 1, mpi_real8, mpi_sum, t%mpi_comm, ierr  )

      if ( my_id == t%root ) write(*,fmt_data) adjustl(t%message), accum_min, accum_max, accum_sum/num_procs

    end subroutine timer_report

end module timer_class
