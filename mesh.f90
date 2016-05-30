module mesh_class
  use data_desc_class
  implicit none

  type mesh
      real(8) :: xl_g, xr_g, yb_g, yt_g
      real(8) :: xl_l, xr_l, yb_l, yt_l

      integer, allocatable :: elem_node(:,:)
      real(8), allocatable :: xy(:,:)
      
      type(data_desc), pointer :: desc

  end type mesh

contains
  
  subroutine mesh_create ( xl, xr, yb, yt, desc, m ) 
    implicit none
    ! Parameters
    real(8), intent(in)                 :: xl, xr, yb, yt
    type(data_desc), intent(in), target :: desc
    type (mesh), intent(out)            :: m

    ! Locals
    integer :: k, j, i 

    ! Point to the descriptor of the distributed FE mesh
    m%desc => desc 

    m%xl_g = xl
    m%xr_g = xr
    m%yb_g = yb
    m%yt_g = yt
   
    m%xl_l = xl 
    m%xr_l = xr
    m%yb_l = yb + ((yt-yb)/(desc%ny_g-1))*((desc%ny_g-1)/desc%np)*desc%me
    m%yt_l = m%yb_l + ((yt-yb)/(desc%ny_g-1))*(desc%ny_l-1) 
    
    ! write (*,'(f8.3,f8.3,f8.3,f8.3)') m%xl_g, m%xr_g, m%yb_g, m%yt_g
    ! write (*,'(f8.3,f8.3,f8.3,f8.3)') m%xl_l, m%xr_l, m%yb_l, m%yt_l 

    allocate ( m%xy(1:2,1:m%desc%node_num) )
    
    k = 0
    do j = 1, m%desc%ny_l
       do i = 1, m%desc%nx_l
          
          k = k + 1
          
          m%xy(1,k) = ( real ( m%desc%nx_l - i,     kind = 8 ) * m%xl_l   &
                      + real (      i - 1, kind = 8 ) * m%xr_l ) &
                      / real ( m%desc%nx_l     - 1, kind = 8 )
          
          m%xy(2,k) = ( real ( m%desc%ny_l - j,     kind = 8 ) * m%yb_l   &
                      + real (      j - 1, kind = 8 ) * m%yt_l ) &
                      / real ( m%desc%ny_l     - 1, kind = 8 )

          ! write (*,'(i,i,f8.3,f8.3)') i, j, m%xy(1,k), m%xy(2,k)
       end do
    end do

    allocate ( m%elem_node(1:3,1:m%desc%elem_num) )
    k = 0

    do j = 1, m%desc%ny_l - 1
       do i = 1, m%desc%nx_l - 1
          
          k = k + 1
          m%elem_node(1,k) = i     + ( j - 1 ) * m%desc%nx_l
          m%elem_node(2,k) = i + 1 + ( j - 1 ) * m%desc%nx_l
          m%elem_node(3,k) = i     +   j       * m%desc%nx_l
          
          ! write (*,'(i,i,i,i)') k,  m%elem_node(1,k), m%elem_node(2,k), m%elem_node(3,k)
    
          k = k + 1
          m%elem_node(1,k) = i + 1 +   j       * m%desc%nx_l
          m%elem_node(2,k) = i     +   j       * m%desc%nx_l
          m%elem_node(3,k) = i + 1 + ( j - 1 ) * m%desc%nx_l

          ! write (*,'(i,i,i,i)') k,  m%elem_node(1,k), m%elem_node(2,k), m%elem_node(3,k)
       end do
    end do
    
  end subroutine mesh_create

  subroutine mesh_destroy ( m ) 
    implicit none
    ! Parameters
    type (mesh), intent(inout)     :: m
   
    ! nullify pointer
    nullify(m%desc)

    m%xl_g = 0.0
    m%xr_g = 0.0
    m%yb_g = 0.0
    m%yt_g = 0.0
   
    m%xl_l = 0.0
    m%xr_l = 0.0
    m%yb_l = 0.0
    
    ! write (*,'(f8.3,f8.3,f8.3,f8.3)') m%xl_g, m%xr_g, m%yb_g, m%yt_g
    ! write (*,'(f8.3,f8.3,f8.3,f8.3)') m%xl_l, m%xr_l, m%yb_l, m%yt_l 

    deallocate ( m%xy )
    deallocate ( m%elem_node )

  end subroutine mesh_destroy

end module mesh_class
