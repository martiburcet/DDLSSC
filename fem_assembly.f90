module fem_assembly
  use mesh_class
  use matrix_class
  use vector_class
  implicit none

  ! Constants
  real (8), parameter :: pi = 3.141592653589793D+00
  
contains

  subroutine assembly_linear_system (msh, A, b) 
    implicit none
    ! Paramters
    type(mesh)  , intent(in)    :: msh
    type(matrix), intent(inout) :: A
    type(vector), intent(inout) :: b

    ! Locals
    integer :: e, i1, i2, i3, ti1, ti2, ti3, nti1, nti2, nti3
    integer :: j1, j2, j3, tj1, tj2, tj3, ntj1, ntj2, ntj3
    integer :: q1, q2, q3, nq1, nq2, nq3
    integer :: k, i, j
    real(8) :: area, xq, yq, wq, qi, dqidx, dqidy, rhs, qj, dqjdx, dqjdy
    real(8) :: u
    
    !
    !  ASSEMBLE THE SYSTEM
    !
    !  Assemble the coefficient matrix A and the right-hand side B of the
    !  finite element equations, ignoring boundary conditions.
    !

    b%data(1:msh%desc%node_num) = 0.0D+00
    A%data(1:msh%desc%node_num,1:msh%desc%node_num) = 0.0D+00

    do e = 1, msh%desc%elem_num

       i1 = msh%elem_node(1,e)
       i2 = msh%elem_node(2,e)
       i3 = msh%elem_node(3,e)

       area = 0.5D+00 * &
            ( msh%xy(1,i1) * ( msh%xy(2,i2) - msh%xy(2,i3) ) &
            + msh%xy(1,i2) * ( msh%xy(2,i3) - msh%xy(2,i1) ) &
            + msh%xy(1,i3) * ( msh%xy(2,i1) - msh%xy(2,i2) ) )
       !
       !  Consider each quadrature point.
       !  Here, we use the midside nodes as quadrature points.
       !
       do q1 = 1, 3

          q2 = mod ( q1, 3 ) + 1

          nq1 = msh%elem_node(q1,e)
          nq2 = msh%elem_node(q2,e)

          xq = 0.5D+00 * ( msh%xy(1,nq1) + msh%xy(1,nq2) )
          yq = 0.5D+00 * ( msh%xy(2,nq1) + msh%xy(2,nq2) )
          wq = 1.0D+00 / 3.0D+00
          !
          !  Consider each test function in the element.
          !
          do ti1 = 1, 3

             ti2 = mod ( ti1,     3 ) + 1
             ti3 = mod ( ti1 + 1, 3 ) + 1

             nti1 = msh%elem_node(ti1,e)
             nti2 = msh%elem_node(ti2,e)
             nti3 = msh%elem_node(ti3,e)

             qi = 0.5D+00 * ( &
                  ( msh%xy(1,nti3) - msh%xy(1,nti2) ) * ( yq - msh%xy(2,nti2) ) &
                  - ( msh%xy(2,nti3) - msh%xy(2,nti2) ) * ( xq - msh%xy(1,nti2) ) ) / area
             dqidx = - 0.5D+00 * ( msh%xy(2,nti3) - msh%xy(2,nti2) ) / area
             dqidy =   0.5D+00 * ( msh%xy(1,nti3) - msh%xy(1,nti2) ) / area

             rhs = 2.0D+00 * pi * pi * sin ( pi * xq ) * sin ( pi * yq )

             b%data(nti1) = b%data(nti1) + area * wq * rhs * qi
             !
             !  Consider each basis function in the element.
             !
             do tj1 = 1, 3

                tj2 = mod ( tj1,     3 ) + 1
                tj3 = mod ( tj1 + 1, 3 ) + 1

                ntj1 = msh%elem_node(tj1,e)
                ntj2 = msh%elem_node(tj2,e)
                ntj3 = msh%elem_node(tj3,e)

                qj = 0.5D+00 * ( &
                     ( msh%xy(1,ntj3) - msh%xy(1,ntj2) ) * ( yq - msh%xy(2,ntj2) ) &
                     - ( msh%xy(2,ntj3) - msh%xy(2,ntj2) ) * ( xq - msh%xy(1,ntj2) ) ) / area
                dqjdx = - 0.5D+00 * ( msh%xy(2,ntj3) - msh%xy(2,ntj2) ) / area
                dqjdy =   0.5D+00 * ( msh%xy(1,ntj3) - msh%xy(1,ntj2) ) / area

                A%data(nti1,ntj1) = A%data(nti1,ntj1) &
                     + area * wq * ( dqidx * dqjdx + dqidy * dqjdy )

             end do

          end do

       end do

    end do

    !
    !  BOUNDARY CONDITIONS
    !
    !  If the K-th variable is at a boundary node, replace the K-th finite
    !  element equation by a boundary condition that sets the variable to U(K).
    !
    k = 0

    do j = 1, msh%desc%ny_l

       do i = 1, msh%desc%nx_l

          k = k + 1

          if ( i == 1 .or. &
               i == msh%desc%nx_l .or. &
               (j == 1 .and.  msh%desc%me == 0) .or. &
               (j == msh%desc%ny_l .and. msh%desc%me == msh%desc%np-1 ) ) then

             call exact ( msh%xy(1,k), msh%xy(2,k), u )

             A%data(k,1:A%n) = 0.0D+00
             A%data(k,k)      = 1.0D+00
             b%data(k)        = u

             b%data(1:k-1) =  b%data(1:k-1) - u * A%data(1:k-1,k)
             b%data(k+1:A%m) =  b%data(k+1:A%m) - u * A%data(k+1:A%m,k)

             A%data(1:k-1,k) = 0.0
             A%data(k+1:A%m,k) = 0.0 

          end if
       end do
    end do


  end subroutine assembly_linear_system

  subroutine exact ( x, y, u )
!*****************************************************************************
!
!  EXACT calculates the exact solution of the Poisson problem at hand
!
!  Discussion:
!
!    The function specified here depends on the problem being
!    solved.  The user must be sure to change both EXACT and RHS
!    or the program will have inconsistent data.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of a point
!    in the region, at which the exact solution is to be evaluated.
!
!    Output, real ( kind = 8 ) U the value of
!    the exact solution U at the point (X,Y).
!
    implicit none
    ! Parameters
    real (8), intent(in)  :: x
    real (8), intent(in)  :: y
    real (8), intent(out) :: u

    u = sin ( pi * x ) * sin ( pi * y ) + x
    
    return
  end subroutine exact


  subroutine check_solution ( msh, x )
!*****************************************************************************
!
!  check_solution compares the computed approximate solution
!  (using linear FEs+PCG) and the exact solution of the 2D Poisson 
!  equation 
!
!
!  Parameters:
!
!    Input, type(mesh)   msh, the FE mesh
!    Input, type(vector) x  , the computed approximate solution 
!

    implicit none
    ! Parameters
    type (mesh)  , intent(in) :: msh
    type (vector), intent(in) :: x 
    
    ! Locals
    integer :: i, j, k
    real(8) :: u 

    !
    !  COMPARE computed and exact solutions.
    !
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     K     I     J          X           Y            U               U             Error'
    write ( *, '(a)' ) '                                                   exact          computed              '
    write ( *, '(a)' ) ' '
    
    k = 0

    do j = 1, msh%desc%ny_l
       do i = 1, msh%desc%nx_l
          
          k = k + 1
          
          call exact ( msh%xy(1,k), msh%xy(2,k), u )
          write ( *, &
               '(2x,i4,2x,i4,2x,i4,2x,f10.2,2x,f10.2,2x,f14.6,2x,f14.6,2x,f14.6)' ) &
               k, i, j, msh%xy(1,k), msh%xy(2,k), u, x%data(k), abs ( u - x%data(k) ) 
          
       end do
       write ( *, '(a)' ) ' '
    end do

  end subroutine check_solution
    
end module fem_assembly
