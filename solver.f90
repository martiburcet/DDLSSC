module solver
  use mesh_class
  use matrix_class
  use precond_class
  use vector_class
  use matrix_vector
  implicit none

  integer, parameter :: i_rtol_ = 1, i_atol_    = 2, & 
                      & i_rn2_  = 3, i_err1_ = 4, i_ub1_ = 5 
  
  integer, parameter :: i_trace_ = 1, i_itmax_ = 2 
  
  integer, parameter :: siz_val_ = 5, siz_ctrl_ = 2

  type itconv_data
    integer :: controls(siz_ctrl_)
    real(8) :: values  (siz_val_)
  end type itconv_data
contains

subroutine itconv_data_init ( me, methdname,trace,itmax,& 
     &                        r, rtol,atol,p_itc_dat)
  implicit none 
  integer, intent(in)            :: me
  character(len=*), intent(in)   :: methdname
  integer, intent(in)            :: trace, itmax
  type(vector), intent(in)       :: r
  real, intent(in)               :: rtol, atol
  type(itconv_data), intent(out) :: p_itc_dat    

  p_itc_dat%controls(:) = 0
  p_itc_dat%values(:)   = 0.0

  p_itc_dat%controls( i_trace_ ) = trace
  p_itc_dat%controls( i_itmax_ ) = itmax


  call vector_nrm2 ( r, p_itc_dat%values(i_rn2_) )
  p_itc_dat%values(i_ub1_) = rtol*p_itc_dat%values(i_rn2_) + atol

  p_itc_dat%values( i_atol_   )  = atol
  p_itc_dat%values( i_rtol_   )  = rtol
  if ((p_itc_dat%controls( i_trace_ ) > 0).and. (me == 0)) &
       &  call log_header(methdname) 

  return

end subroutine itconv_data_init

subroutine log_header( methdname )
  implicit none 

  ! Parameters
  character(len=*), intent(in)   :: methdname

  ! Local variables
  character(len=*), parameter    :: fmt1='(a18,1x,a4,3(2x,a15))'
  character(len=*), parameter    :: fmt2='(a18,1x,a4,3(2x,a15),3(2x,a15))'
  integer, parameter             :: outlen=18 
  character(len=len(methdname))  :: mname
  character(len=outlen)          :: outname
  
  mname = adjustl(trim(methdname))
  write(outname,'(a)') mname(1:min(len_trim(mname),outlen-1))//':'
  write(*,fmt1) adjustl(outname),'Iteration','Error Estimate','Tolerance'

end subroutine log_header

subroutine log_conv ( me, methdname, it, trace, err1, ub1 )
  implicit none 

  ! Parameters
  integer, intent(in)           :: me
  character(len=*), intent(in)  :: methdname
  integer, intent(in)           :: it
  integer, intent(in)           :: trace
  real, intent(in)              :: err1, ub1


  ! Local variables
  character(len=*), parameter   :: fmt1='(a18,1x,i4,3(2x,es16.9))'
  character(len=*), parameter   :: fmt2='(a18,1x,i4,3(2x,es16.9),3(2x,es16.9))'
  integer, parameter            :: outlen=18 
  character(len=len(methdname)) :: mname
  character(len=outlen)         :: outname

  if ( (trace/=0) .and.(mod(it,trace) == 0).and.(me == 0)) then 
     mname = adjustl(trim(methdname))
     write(outname,'(a)') mname(1:min(len_trim(mname),outlen-1))//':'
     write(*,fmt1) adjustl(outname), it, err1, ub1
  endif

end subroutine log_conv


subroutine itconv_data_end ( me, methdname, it, p_itc_dat, err1, iter, conv )
  !use psb_base_mod
  implicit none
  integer, intent(in)               :: me
  character(len=*), intent(in)      :: methdname
  integer, intent(in)           :: it
  type(itconv_data), intent(in)     :: p_itc_dat

  real, optional, intent(out)   :: err1
  integer,  optional, intent(out)   :: iter
  logical,  optional, intent(out)   :: conv

  character(len=*), parameter  :: fmt11='(a,2x,es16.9,1x,a,1x,i4,1x,a)'
  character(len=*), parameter  :: fmt12='(a,3(2x,es16.9))'

  character(len=*), parameter  :: fmt21='(a,2x,es16.9,1x,es16.9,1x,a,1x,i4,1x,a)'
  character(len=*), parameter  :: fmt22='(a,3(2x,es16.9),3(2x,es16.9))'

  if ( present(conv) ) conv = .true.
     if ( p_itc_dat%values(i_err1_) > p_itc_dat%values(i_ub1_) ) then
        if ( me == 0 )  then 
           write(*,fmt11) trim(methdname)//' failed to converge to ', p_itc_dat%values(i_ub1_),&
                & ' in ',it,' iterations. '
           write(*,fmt12) 'Last iteration error estimate: ',&
                & p_itc_dat%values(i_err1_)
        end if
        if (present(conv)) conv = .false.
     else
        if ( me == 0 .and. p_itc_dat%controls(i_trace_) /= 0  )  then 
           write(*,fmt11) trim(methdname)//' converged to ', p_itc_dat%values(i_ub1_),&
                & ' in ',it,' iterations. '
           write(*,fmt12) 'Last iteration error estimate: ',&
                & p_itc_dat%values(i_err1_)
        end if
     end if
     if ( present(err1) ) err1 = p_itc_dat%values(i_err1_)
     if ( present(iter) ) iter = it

end subroutine itconv_data_end


function itconv_data_check ( me, methdname, it, r, p_itc_dat )
  implicit none 
  ! Parameters
  logical                              :: itconv_data_check
  integer, intent(in)                  :: me
  character(len=*)     , intent(in)    :: methdname
  integer              , intent(in)    :: it
  type(vector)         , intent(in)    :: r
  type(itconv_data)    , intent(inout) :: p_itc_dat

  itconv_data_check = .false.

  ! Compute || r(i) ||p_itc_dat%values(i_err1_)
  call vector_nrm2 ( r, p_itc_dat%values(i_err1_) )
  
  itconv_data_check = ( p_itc_dat%values(i_err1_) <= p_itc_dat%values(i_ub1_) )

  itconv_data_check = (itconv_data_check .or. (p_itc_dat%controls(i_itmax_) <= it))

  if ( (p_itc_dat%controls(i_trace_) > 0).and.&
       & ((mod(it,p_itc_dat%controls(i_trace_)) == 0).or.itconv_data_check)) then 
     call log_conv(me,methdname,it,p_itc_dat%controls(i_trace_), &
          p_itc_dat%values(i_err1_), p_itc_dat%values(i_ub1_) )  
  end if

  return
end function itconv_data_check


!=============================================================================
! Preconditioned Conjugate Gradient
!=============================================================================
subroutine pcg( A,M,b,x,trace,itmax,rtol,atol,it,err,conv)
  !-----------------------------------------------------------------------------
  ! This routine performs pcg iterations on Ax=b with preconditioner M. 
  !-----------------------------------------------------------------------------
  implicit none

  ! Mandatory parameters
  type(matrix) , intent(inout) :: A        ! Matrix
  type(precond), intent(inout) :: M        ! Preconditioner
  type(vector) , intent(in)    :: b        ! RHS
  type(vector) , intent(inout) :: x        ! Approximate solution

  ! Optional parameters
  integer, intent(in)   :: trace       ! Message every trace iterations
  integer, intent(in)   :: itmax       ! Max. # of iterations
  real   , intent(in)   :: rtol        ! Relative tolerance
  real   , intent(in)   :: atol        ! Absolute tolerance
  integer, intent(out)  :: it          ! # of iterations to converge
  real   , intent(out)  :: err         ! Error estimates
  logical, intent(out)  :: conv        ! Converged?


  ! Locals
  type (itconv_data) :: stopdat      ! Object to store conv. info
  integer        :: trace_       ! Message every trace iterations
  integer        :: itmax_       ! Max. # of iterations
  real           :: rtol_        ! Relative tolerance
  real           :: atol_        ! Absolute tolerance
  real           :: r_z, Ap_p, alpha, beta
  integer        :: iter
  integer            :: me, np
  type(vector)       :: r,p,Ap,z     ! Working vectors


  ! 0) Initialize optional arguments
  trace_ = trace
  itmax_ = itmax
  rtol_  = rtol
  atol_  = atol

  call vector_create ( A%desc, r  )
  call vector_create ( A%desc, Ap )
  call vector_create ( A%desc, z  )
  call vector_create ( A%desc, p  )

  me = A%desc%me 
  np = A%desc%np

  ! 1) Compute initial residual
  ! 1.a) r=Ax
  call matvec ( A, x, r )

  ! 1.b) r=b-r
  call vector_axpby( 1.0, -1.0, b, r )

  ! 2) z=inv(M)r
  call precond_apply ( M , r , z )

  ! 3) <r,z>
  call vector_dot( r, z, r_z )

  ! 4) Initializations:
  ! p=z
  call vector_copy(z,p) 

  ! Register convergence info
  call itconv_data_init (   me  , &
                           'PCG', &
                           trace_, &
                           itmax_, & 
                           r, &
                           rtol_, &
                           atol_, &
                           stopdat  )
  ! 5) Iteration
  iter = 0
  loop_pcg: do
     iter = iter + 1

 
     ! Ap = A*p
     call matvec(A,p,Ap)

     ! <Ap,p>
     call vector_dot (Ap, p, Ap_p)

     if (Ap_p /= 0.0) then
        alpha = r_z / Ap_p
     else
        alpha = 0.0
     end if

     ! x = x + alpha*p
     call vector_axpby ( alpha, 1.0, p, x )

     ! r = r - alpha*Ap
     call vector_axpby (-alpha, 1.0, Ap, r )

     if (itconv_data_check ( me  ,  &
                             'PCG',  &  
                             iter,  & 
                             r,  &
                             stopdat  ) ) then
        exit loop_pcg
     end if

     ! z = inv(M) r
     call precond_apply(M,r,z)

     beta = 1.0/r_z
     call vector_dot(r,z,r_z)
     beta = beta*r_z

     ! p = z + beta*p
     call vector_axpby(1.0,beta,z,p)

  end do loop_pcg

  call vector_destroy ( r  )
  call vector_destroy ( Ap )
  call vector_destroy ( z  )
  call vector_destroy ( p  )

  ! Return convergence info
  call itconv_data_end ( me, 'PCG', iter, stopdat, &  
                         err, it, conv )

end subroutine pcg
    
end module solver
