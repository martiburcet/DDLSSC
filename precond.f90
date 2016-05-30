module precond_class
  use data_desc_class
  use matrix_class
  use vector_class
  implicit none

  type precond
     type(data_desc), pointer :: desc        ! Pointer to the FE mesh descriptor
     integer                  :: type
     type(vector)             :: inv_diag    ! Stores the inverse of the diagonal of A
     real(8), allocatable     :: chol(:,:)   ! Stores the cholesky factor of A in its 
                                             ! upper triangle
  end type precond

  integer, parameter  :: identity = 0  ! M = I
  integer, parameter  :: diagonal = 1  ! M = D, where D=diag(A)
  integer, parameter  :: neumann  = 2  ! M = sum_i=1^P = I_i A_i^-1 I_i^T

  interface

!!$      SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
!!$
!!$  -- LAPACK routine (version 3.1) --
!!$    Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!!$     November 2006
!!$
!!$     .. Scalar Arguments ..
!!$      CHARACTER          UPLO
!!$      INTEGER            INFO, LDA, N
!!$     ..
!!$     .. Array Arguments ..
!!$      DOUBLE PRECISION   A( LDA, * )
!!$     ..
!!$
!!$  Purpose
!!$  =======
!!$
!!$  DPOTRF computes the Cholesky factorization of a real symmetric
!!$  positive definite matrix A.
!!$
!!$  The factorization has the form
!!$     A = U**T * U,  if UPLO = 'U', or
!!$     A = L  * L**T,  if UPLO = 'L',
!!$  where U is an upper triangular matrix and L is lower triangular.
!!$
!!$  This is the block version of the algorithm, calling Level 3 BLAS.
!!$
!!$  Arguments
!!$  =========
!!$
!!$  UPLO    (input) CHARACTER*1
!!$          = 'U':  Upper triangle of A is stored;
!!$          = 'L':  Lower triangle of A is stored.
!!$
!!$  N       (input) INTEGER
!!$          The order of the matrix A.  N >= 0.
!!$
!!$  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!!$          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!!$          N-by-N upper triangular part of A contains the upper
!!$          triangular part of the matrix A, and the strictly lower
!!$          triangular part of A is not referenced.  If UPLO = 'L', the
!!$          leading N-by-N lower triangular part of A contains the lower
!!$          triangular part of the matrix A, and the strictly upper
!!$          triangular part of A is not referenced.
!!$
!!$          On exit, if INFO = 0, the factor U or L from the Cholesky
!!$          factorization A = U**T*U or A = L*L**T.
!!$
!!$  LDA     (input) INTEGER
!!$          The leading dimension of the array A.  LDA >= max(1,N).
!!$
!!$  INFO    (output) INTEGER
!!$          = 0:  successful exit
!!$          < 0:  if INFO = -i, the i-th argument had an illegal value
!!$          > 0:  if INFO = i, the leading minor of order i is not
!!$                positive definite, and the factorization could not be
!!$                completed.
!!$
!!$  =====================================================================
!!$
     SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
       implicit none
       ! .. Scalar Arguments ..
       CHARACTER, INTENT(IN)     :: UPLO
       INTEGER  , INTENT(IN)     :: INFO, LDA, N
       ! *     ..
       !       .. Array Arguments ..
       REAL(8) , INTENT(INOUT)  :: A( LDA, * )
     END SUBROUTINE DPOTRF

!!$      SUBROUTINE DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
!!$
!!$  -- LAPACK routine (version 3.1) --
!!$     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!!$     November 2006
!!$
!!$     .. Scalar Arguments ..
!!$      CHARACTER          UPLO
!!$      INTEGER            INFO, LDA, LDB, N, NRHS
!!$     ..
!!$     .. Array Arguments ..
!!$      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!!$     ..
!!$
!!$  Purpose
!!$  =======
!!$
!!$  DPOTRS solves a system of linear equations A*X = B with a symmetric
!!$  positive definite matrix A using the Cholesky factorization
!!$  A = U**T*U or A = L*L**T computed by DPOTRF.
!!$
!!$  Arguments
!!$  =========
!!$
!!$  UPLO    (input) CHARACTER*1
!!$          = 'U':  Upper triangle of A is stored;
!!$          = 'L':  Lower triangle of A is stored.
!!$
!!$  N       (input) INTEGER
!!$          The order of the matrix A.  N >= 0.
!!$
!!$  NRHS    (input) INTEGER
!!$          The number of right hand sides, i.e., the number of columns
!!$          of the matrix B.  NRHS >= 0.
!!$
!!$  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!!$          The triangular factor U or L from the Cholesky factorization
!!$          A = U**T*U or A = L*L**T, as computed by DPOTRF.
!!$
!!$  LDA     (input) INTEGER
!!$          The leading dimension of the array A.  LDA >= max(1,N).
!!$
!!$  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!!$          On entry, the right hand side matrix B.
!!$          On exit, the solution matrix X.
!!$
!!$  LDB     (input) INTEGER
!!$          The leading dimension of the array B.  LDB >= max(1,N).
!!$
!!$  INFO    (output) INTEGER
!!$          = 0:  successful exit
!!$          < 0:  if INFO = -i, the i-th argument had an illegal value
!!$
!!$  =====================================================================
!!$
     SUBROUTINE DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
       implicit none
       ! .. Scalar Arguments ..
       CHARACTER, INTENT(IN)     :: UPLO
       INTEGER  , INTENT(IN)     :: N, NRHS, LDA, LDB, INFO
       ! *     ..
       !       .. Array Arguments ..
       REAL(8) , INTENT(IN)    :: A( LDA, * )
       REAL(8) , INTENT(INOUT) :: B( LDB, * )
     END SUBROUTINE DPOTRS

  end interface 


contains
  
  subroutine precond_create ( A, type, prec ) 
    implicit none
    ! Parameters
    type (matrix)  , intent(in), target :: A
    integer        , intent(in)         :: type
    type (precond) , intent(out)        :: prec

    ! Locals
    integer              :: ierr
    integer		 		 :: i
    integer 			 :: j
    type (vector)        :: diag
    type (vector)        :: diagGlob

    ! Primero creamos el vector para la matriz diagonal y para la inversa de la diagonal, pasando la información de la malla a través de la 	descripción de A.
    call vector_create (A%desc, diag)
    call vector_create (A%desc, prec%inv_diag)

    ! Alocatamos cholesky dandole las dimensiones de la matriz A y lo inicializamos.
    allocate ( prec%chol(1:A%m,1:A%n) )
    prec%chol = 0.0d0

    ! La descripción de la malla del precondicionador apunta a la descripción de A (siendo A un puntero a la información de la malla).
    ! Entendemos que con => pasamos una direccion de memoria, no un valor.
    prec%desc => A%desc
    ! Type indica el tipo de precondicionador que queremos usar, identidad, diagonal o neuman
    prec%type = type
    
    ! El precondicionador identidad no esta dentro del tipo precond, por lo que hemos supuesto que no se implementa, ya que si se utiliza
    ! la identidad como precondicionador, al multiplicar por A, el resultado es A. Así que si selecionamos type=0=identity, simplemente
    ! calculamos sin precondicionador. Lo hemos deducido porque no está en el tipo precond.
    
    ! Calculamos la matriz diagonal
    if (type == diagonal) then        
		do i=1,A%m
            diag%data(i) = A%data(i,i)
        end do
    ! Llamamos al comunicador del vector para obtener el vector global de elementos de la diagonal, con la aportacion de todas las 
    ! contribuciones de las interfaces.
        call vector_comm(diag, diagGlob)
    ! Invertimos diagGlob y obtenemos el precondicionador, que es la inversa de la diagonal.
        do i=1,A%m            
            prec%inv_diag%data(i) = 1.0d0/diagGlob%data(i)
        end do
    ! Para el caso de neumann se realiza una cholesky factorization para poder invertir la matriz mas facilmente. Para ello se utilizan
    ! las librerias LAPACK que se enuncian en la parte superior del código.
    else if (type == neumann) then
	! Asignamos los valores de la matriz al precondicionador de cholesky
        prec%chol = A%data
        call DPOTRF('U', A%n, prec%chol, A%m, ierr)
    end if

  end subroutine precond_create

  subroutine precond_apply (M,r,z)
    implicit none
    ! Parameters
    type (precond), intent(in)        :: M
    type (vector) , intent(in)        :: r
    type (vector) , intent(inout)     :: z

    ! Locals
    integer              :: ierr,i
    type (vector)        :: zsub ! z subassembled. Antes de aplicar los pesos.
    
    ! Aplicamos el precondicionador. Primero utilizamos vector_comm para sumar todos los residuales locales.
    ! z=Mr, en caso de la identidad directamente z=r
    if (M%type == identity) then
        call vector_comm(r, z)
    ! En el caso de diagonal, z=Mr=diag(A)r. Primero calculamos z=r(suma de todas las contribuciones) con el vector_comm y 
    !luego multiplicamos por el precondicionador:
    ! z=diag(A)*z. 
    else if (M%type == diagonal) then
        call vector_comm(r, z)
        do i=1,z%desc%node_num
            z%data(i) = M%inv_diag%data(i) * z%data(i)
        end do
    ! Neumann. Calculamos la inversa del precondicionador usando la libreria LAPACK. Como hemos factorizado M con Lower triangular, pues
    ! utilizaremos lower cholesky factorization.
    else if (M%type == neumann) then
    ! Primero hacemos la suma de la contribución de todos los residuos, como en los apartados anteriores. En este caso al ser neuman-neuman 
    ! deberemos aplicar la matriz de pesos: E=WR'. y de momento solo hemos multiplicado por los restriction operators R. A esta suma de 
    ! contribución de todos los subdominios la llamaremos zsub (z subassembled)         
		call vector_comm(r, zsub)
		call vector_weight (zsub)
    ! Utilizando LAPACK calculamos cholesky factorization para invertir.
        call DPOTRS('U', r%desc%node_num, 1, M%chol, r%desc%node_num, zsub%data, zsub%desc%node_num, ierr)
    ! Una vez calculado la nueva zsub, hay que ensamblarla y aplicar otra vez el weight operator	
        call vector_comm(zsub, z)
        call vector_weight(z)
    end if
  end subroutine precond_apply

  subroutine precond_destroy ( prec ) 
    implicit none
    ! Parameters
    type (precond), intent(inout)     :: prec
  
    if ( prec%type == diagonal ) then
      call vector_destroy ( prec%inv_diag ) 
    else if ( prec%type == neumann ) then
      deallocate(prec%chol)
    end if
 
  end subroutine precond_destroy

end module precond_class
