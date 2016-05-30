module vector_class
  use mpi
  use data_desc_class
  implicit none

  type vector
     integer                  :: n         ! Number of vector entries
     real(8), allocatable     :: data(:)   ! Vector entries
     type(data_desc), pointer :: desc      ! Pointer to the FE mesh descriptor
  end type vector

  integer, parameter :: MAX_TAG      = 10000
  integer            :: current_tag  = 0 

contains
  
  subroutine vector_create ( desc, v )
    implicit none
    ! Parameters
    type(data_desc), intent(in), target :: desc
    type (vector)  , intent(out)        :: v

    ! Locals
    integer :: k, j, i 

    ! Point to the descriptor of the distributed FE mesh
    v%desc => desc 
    
    v%n = v%desc%node_num
    allocate ( v%data(1:v%n) )
    
    v%data = 0.0

  end subroutine vector_create

  subroutine vector_axpby ( alpha, beta, x, y )
!*****************************************************************************
!
!  vector_axpby performs the vector update y = beta*y + alpha*x,
!  where beta and alpha are scalars; x and y are vectors.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) alpha, beta, scalars involved in the
!    vector operator
!   
!    Input, type(vector) x. Input vector x.
!
!    Input/Output, type(vector) y. Updated vector y.
!
    implicit none
    ! Parameters
    real(8)     , intent(in)    :: alpha, beta
    type(vector), intent(in)    :: x
    type(vector), intent(inout) :: y

    y%data = beta * y%data + alpha * x%data
 
  end subroutine vector_axpby

  subroutine vector_copy ( x, y )
!*****************************************************************************
!
!  vector_copy copies x into y 
!
!  Parameters:
!
!    Input, type(vector) x. Input vector x.
!
!    Input/Output, type(vector) y. Vector y resulting from the copy.
!
    implicit none
    ! Parameters
    type(vector), intent(in)    :: x
    type(vector), intent(inout) :: y

    y%data = x%data
 
  end subroutine vector_copy

  subroutine vector_dot ( x, y, alpha )
    implicit none
    ! Parameters 
    type(vector), intent(in)  :: x
    type(vector), intent(in)  :: y
    real(8)     , intent(out) :: alpha
    
    ! Locals
    integer :: ierr
    integer :: i

    ! Usando la función DOT_PRODUCT implementaremos y(tras)*x

    alpha = DOT_PRODUCT(y%data, x%data)
    
    ! Hacemos un MPI_allreduce para sumar todas las contribuciones en el dot product de los subdominios.
   call MPI_ALLREDUCE(MPI_IN_PLACE,         & ! Send data
                  alpha,                & ! Recv data
                  1,                    & ! dimension
                  MPI_DOUBLE_PRECISION, & ! Type
                  MPI_SUM,              & ! Operation
                  x%desc%mpi_comm,      & ! Communicator
                  ierr)

  end subroutine vector_dot


  subroutine vector_nrm2 ( x, nrm2 )
    implicit none
    ! Parameters 
    type(vector), intent(in)  :: x
    real(8)     , intent(out) :: nrm2

    ! Locals
    integer              :: ierr
    integer				 :: i
    type (vector)		 :: fullx

    nrm2=0.d0

    ! La subroutine vector_dot multiplica dos vectores, por lo que, dándole como parámetros de entrada fullx & fullx 
    !(ya que x es partial sumed) tendremos su multiplicación que haciendo la raíz cuadrada, finalmente, obtendremos la norma.
    call vector_comm (x, fullx)
    call vector_dot (fullx, fullx, nrm2)
    nrm2= sqrt(nrm2)

  end subroutine vector_nrm2

  subroutine vector_comm ( x, y)
  !*****************************************************************************
!
!  vector_comm sums-up all subdomain contributions (partially 
!  summed entries) to vector entries of x shared among several 
!  subdomains, i.e., vector entries corresponding to FE mesh nodes
!  laying on the interface, and stores the sum into y. In the case 
!  we are facing (1D FE mesh data distribution), a given node on the 
!  interface is shared among at most 2 different subdomains.
!
!  Parameters:
!
!    Input, type(vector) x. Input (partially summed) vector x. 
!
!    Input/Output, type(vector) y. Output (fully summed) vector y.  

    implicit none
    ! Parameters 
    type(vector), intent(in)     :: x
    type(vector), intent(inout)  :: y
       
    ! Locals
    integer:: ierr
    integer istat(MPI_STATUS_SIZE)
    real(8), allocatable::uprecv(:)
    real(8), allocatable::dwnrecv(:)
    real(8), allocatable::upsend(:)
    real(8), allocatable::dwnsend(:)

    ! Alocatamos todas las variables que van a compartir la información entre los subdominios.

    allocate(uprecv(1:x%desc%nx_l+1))
    allocate(dwnrecv(1:x%desc%nx_l+1))
    allocate(upsend(1:x%desc%nx_l))
    allocate(dwnsend(1:x%desc%nx_l))

    ! Usando la subroutine create_vector, creamos los vectores upsend and dwnsend que contendrán la información de los nodos.
    ! Es importante comentar que para declarar el vector upsend, es necesario identificar cual será el primer nodo local dentro del 
    ! dominio que mandará la información y cual será el último. 
    ! El problema encontrado en upsend, se evita en el caso de dwnsend, ya que siempre empezará en 1 y terminará en nx locales.

    call vector_create(x%desc, y)
    upsend = x%data(x%desc%up_off:x%desc%node_num)
    dwnsend = x%data(1:x%desc%nx_l)


    ! En el caso en que tengamos más de un procesador, es decir, nuestro dominio está subdivido en varios subdominios, se deberán 
    ! comunicar entre ellos.

    ! Nuestro vector y de salida, tendrá el mismo tamaño que nuestro vector x de entrada por la aplicación de los restrictors operators.
    call vector_copy(x,y)

    ! Hacemos las comunicaciones de nuestros nodos en caso de que haya más de un procesador a través de un Sendrecv.
    if (x%desc%np.gt.1) then
		call MPI_Sendrecv(upsend,  			& ! send buff
			x%desc%nx_l, 			& ! dimension
			MPI_DOUBLE_PRECISION,		& ! type
			x%desc%up_pid,			& ! send to ...
			0,				& ! stag
			dwnrecv, 			& ! recv buff
			x%desc%nx_l, 			& ! dimension
			MPI_DOUBLE_PRECISION, 		& ! type		
			x%desc%down_pid,  		& ! recv from ...
			0,  				& ! rtag
			x%desc%mpi_comm,		& ! communicator
			istat, ierr)

		call MPI_Sendrecv(dwnsend,  			& ! send buff
			x%desc%nx_l, 			& ! dimension
			MPI_DOUBLE_PRECISION,		& ! type
			x%desc%down_pid,		& ! send to ...
			0,				& ! stag
			uprecv, 			& ! recv buff
			x%desc%nx_l, 			& ! dimension
			MPI_DOUBLE_PRECISION, 		& ! type		
			x%desc%up_pid,  		& ! recv from ...
			0,  				& ! rtag
			x%desc%mpi_comm,		& ! communicator
			istat, ierr)

	! Parte superior del subdominio 0.
		if (x%desc%me == 0) then
            y%data(x%desc%up_off:x%desc%node_num) = y%data(x%desc%up_off:x%desc%node_num) + uprecv

	! Parte inferior del mayor subdominio.
        else if (x%desc%me == x%desc%np-1) then
            y%data(1:x%desc%nx_l) = y%data(1:x%desc%nx_l) + dwnrecv
	
	! Dominios intermedios tanto inferior como superiormente.        
		else
            y%data(1:x%desc%nx_l) = y%data(1:x%desc%nx_l) + dwnrecv
            y%data(x%desc%up_off:x%desc%node_num) = y%data(x%desc%up_off:x%desc%node_num) + uprecv
        end if

    end if

    ! Desalocatamos todas las variables que hemos creado para liberar memoria :-)
    deallocate(uprecv)
    deallocate(dwnrecv)
    deallocate(upsend)
    deallocate(dwnsend)
   

   
  end subroutine vector_comm


  subroutine vector_weight ( x )
!*****************************************************************************
!
!  vector_weight multiplies each component of x shared among
!  several subdomains (i.e., vector entries corresponding to FE 
!  mesh nodes laying on the interface) by the inverse of the 
!  number of subdomains (i.e., by 1/2) 
!
!  Parameters:
!
!    Input/Output, type(vector) x. Output vector x.  
!
    implicit none
    ! Parameters 
    type(vector), intent(inout)  :: x 

    ! Locals
    integer              :: ierr

    ! En el caso de que haya mas de un procesador, es decir tenemos varios subdominios.
     if (x%desc%np.gt.1) then
    ! Para el subodominio 0, solo tenemos interface superior, por lo que, se dividirán los nodos superiores.
		if (x%desc%me == 0) then
    	    x%data(x%desc%up_off:x%desc%node_num) = x%data(x%desc%up_off:x%desc%node_num)/2.0d0
    ! Para el subodominio (np-1) superior, solo tenemos interface inferior, por lo que, se dividirán los nodos inferiores.
    	else if (x%desc%me == x%desc%np-1) then
    	    x%data(1:x%desc%nx_l) = x%data(1:x%desc%nx_l)/2.0d0
    ! Para subdominios intermedios deberemos dividir tanto en nodos superiores como inferiores.
    	else
    	    x%data(1:x%desc%nx_l) = x%data(1:x%desc%nx_l)/2.0d0
    	    x%data(x%desc%up_off:x%desc%node_num) = x%data(x%desc%up_off:x%desc%node_num)/2.0d0
    	end if
     end if

  end subroutine vector_weight


  subroutine vector_destroy ( v ) 
    implicit none
    ! Parameters
    type (vector), intent(inout)     :: v
   
    v%n = -1
    deallocate ( v%data )

  end subroutine vector_destroy

end module vector_class
