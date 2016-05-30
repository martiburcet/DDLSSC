module matrix_vector
  use matrix_class
  use vector_class
  use mpi
  implicit none

contains

  subroutine matvec ( A, x, y )
    implicit none
    ! Parameters 
    type(matrix), intent(in)    :: A 
    type(vector), intent(in)    :: x
    type(vector), intent(inout) :: y

    ! Locals
    integer              :: ierr
    integer 		     :: i
    real(8), allocatable :: psum(:)

    y%data(:)=0.0d0
    
     ! La idea, es multiplicar por columnas ya que en fortran se obtiene mayor rendimiento debido al storage del propio programa
     
#ifdef ENABLE_OPENMP

    !$omp parallel default(none) shared(x,A,y) private(i,psum)  
	allocate (psum(A%m))
    psum = 0.0d0
     ! Number of rows in the matrix (m) , Number of cols in the matrix (n)
     ! Debido a que la operación 'reduction' no se puede aplicar ya que sólo es aplicable a constantes y, en nuestro caso es un vector, vamos 	   a utilizar lo orden 'critical' pag 22. http://www.openmp.org/presentations/miguel/F95_OpenMPv1_v2.pdf

     !$omp do schedule(static)
 	 do i=1, A%n 
 		psum = psum + A%data(:,i) * x%data(i)	
  	 end do
     !$omp end do
     
     !$omp critical      
     y%data=y%data+psum
     !$omp end critical
     
     deallocate (psum)
     !$omp end parallel 
#else 
    do i=1,A%n
        y%data = y%data+ A%data(:,i) * x%data(i) 
    end do
#endif
     
  end subroutine matvec

end module matrix_vector
