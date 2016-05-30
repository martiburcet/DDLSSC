program hello_world
  use mesh_class
  use data_desc_class

#ifdef ENABLE_OPENMP
  use omp_lib 
#endif

  implicit none

  ! Locals
  integer :: info, me, np

#ifdef ENABLE_OPENMP
  integer :: tid
#endif

  call mpi_init(info) 
  if (info /= mpi_success) then
     write(0,*) 'Error in initalizing MPI, bailing out',info 
     call mpi_abort(mpi_comm_world, -1, info)
  end if

  call mpi_comm_size(mpi_comm_world,np,info)
  call mpi_comm_rank(mpi_comm_world,me,info)
  
#ifdef ENABLE_OPENMP
!$OMP PARALLEL PRIVATE(tid) SHARED (me,np) 
  tid = omp_get_thread_num()
  write(*,'(a,i2,a,i2,a,i2,a)') 'Hello world from thread', tid, ' process id', me, ' of', np, ' processors'
!$OMP END PARALLEL
#else
  write(*,'(a,i2,a,i2,a)') 'Hello world from process id', me, ' of', np, ' processors'
#endif

  call mpi_finalize(info)

end program hello_world
