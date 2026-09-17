module set_precision
  use iso_fortran_env, only : real64, int8, int16, int32, int64
  implicit none
  private
  public :: dp, i1, i2, i4, i8
  integer, parameter :: dp  = real64
  integer, parameter :: i1  = int8
  integer, parameter :: i2  = int16
  integer, parameter :: i4  = int32
  integer, parameter :: i8  = int64
end module set_precision

module set_constants
  use set_precision, only : dp
  implicit none
  private
  public :: zero, one
  public :: max_text_line_length
  real(dp), parameter :: zero      = 0.0_dp
  real(dp), parameter :: one       = 1.0_dp
  integer,  parameter :: max_text_line_length = 1024
end module set_constants

module timer_derived_type

  use set_precision, only : dp
  use set_constants, only : zero
  implicit none
  private
  public :: basic_timer_t

  type :: basic_timer_t
    private
    real(dp)         :: time_start   = zero
    real(dp), public :: time_elapsed = zero
  contains
    private
    procedure, public, pass :: tic => timer_tick
    procedure, public, pass :: toc => timer_tock
  end type basic_timer_t

contains

  function get_time()
    integer(kind=8) :: ticks, ticks_per_sec, max_ticks
    real(dp) :: get_time
    call system_clock( count      = ticks,                                     &
                      count_rate = ticks_per_sec,                              &
                      count_max  = max_ticks )
    if ( ticks_per_sec == 0 ) then
      get_time = zero
    else
      get_time = real(ticks,dp) / real(ticks_per_sec,dp)
    end if
  end function get_time

  subroutine timer_tick( this )
    class(basic_timer_t), intent(inout) :: this
    this%time_elapsed = zero
    this%time_start   = get_time()
  end subroutine timer_tick

  function timer_tock( this )
    class(basic_timer_t), intent(in) :: this
    real(dp)                         :: timer_tock
    timer_tock = get_time() - this%time_start
  end function timer_tock

end module timer_derived_type

module string_stuff
  implicit none
  private
  public :: generate_newline_string
  public :: write_integer_tuple
contains
  subroutine write_integer_tuple(integer_list,out_string,plus,delim)
    integer, dimension(:),  intent(in)  :: integer_list
    character(*),           intent(out) :: out_string
    logical,      optional, intent(in)  :: plus
    character(*), optional, intent(in)  :: delim
    integer :: j, sz
    character(1) :: p
    p = ' '
    if ( present(plus) ) then
      if ( plus ) p = '+'
    end if
    sz = size(integer_list)
    out_string=''
    if (sz<1) return
    if ( integer_list(1)>0 ) then
      write(out_string,'(A,I0)') trim(out_string)//p,integer_list(1)
    elseif ( integer_list(1)==0 ) then
      write(out_string,'(A,I0)') trim(out_string)//' ',integer_list(1)
    else
      write(out_string,'(A,I0)') trim(out_string),integer_list(1)
    end if
    if ( present(delim) ) then
      do j = 2,sz
        if ( integer_list(j)>0 ) then
          write(out_string,'(A,I0)') trim(out_string)//delim//p,integer_list(j)
        elseif ( integer_list(j)==0 ) then
          write(out_string,'(A,I0)') trim(out_string)//delim//' ',integer_list(j)
        else
          write(out_string,'(A,I0)') trim(out_string)//delim,integer_list(j)
        end if
      end do
    else
      do j = 2,sz
        if ( integer_list(j)>0 ) then
          write(out_string,'(A,I0)') trim(out_string)//', '//p,integer_list(j)
        elseif ( integer_list(j)==0 ) then
          write(out_string,'(A,I0)') trim(out_string)//',  ',integer_list(j)
        else
          write(out_string,'(A,I0)') trim(out_string)//', ',integer_list(j)
        end if
      end do
    end if
  end subroutine write_integer_tuple
  
  subroutine generate_newline_string(strings,out_fmt)
    character(*), dimension(:), intent(in)  :: strings
    character(*)              , intent(out) :: out_fmt
    integer :: j, sz
    sz = size(strings)
    out_fmt = '('//trim(strings(1))
    do j=2,sz
        out_fmt=trim(out_fmt)//',/,'//trim(strings(j))
    end do
    out_fmt=trim(out_fmt)//')'
  end subroutine generate_newline_string
end module string_stuff

module index_conversion
  implicit none
  private
  public :: global2local  
contains
  pure function global2local(iG,nSub) result(iSub)
    integer,               intent(in) :: iG
    integer, dimension(:), intent(in) :: nSub
    integer, dimension(size(nSub)) :: iSub
    integer :: i, nDims, p, iGtmp, iTmp
    nDims = size(nSub)
    if (nDims==1) then
      iSub(1) = iG
      return
    end if
    p = product(nSub)
    iGtmp = iG
    do i = nDims,1,-1
      p = p/nSub(i)
      iTmp = mod(iGtmp-1,p) + 1
      iSub(i) = (iGtmp-iTmp)/p + 1
      iGtmp = iTmp
    end do
  end function global2local
end module index_conversion

module combinatorics
  implicit none
  private
  public :: rand_int_in_range
  public :: nchoosek
  public :: get_exponents
  public :: mono_rank_grlex
  public :: mono_rank_grlex2
  public :: mono_rank_grlex3
  public :: mono_unrank_grlex
  public :: mono_unrank_grlex2
  public :: mono_unrank_grlex3
  public :: get_diff_idx
  public :: get_deriv_coef_v1, get_deriv_coef_v2
  ! public :: fact

  integer, parameter :: n_max = 10
  integer :: i_, j_
  integer, parameter, dimension(0:n_max)   :: int_   = [ 1, ( i_, i_ = 1, n_max ) ]
  integer, parameter, dimension(0:n_max)   :: fact_  = [ (product(int_(0:i_)), i_ = 0, n_max ) ]
  integer, parameter, dimension(0:2*n_max) :: pfact_ = [ ( 1, i_ = 0, n_max-1 ), fact_ ]
  integer, parameter, dimension(0:n_max, 0:n_max):: fall_ = reshape( &
          [( ( fact_(i_) / merge( pfact_(n_max + i_ - j_), 1, j_ <= i_ ), &
           i_ = 0, n_max ), j_ = 0, n_max ) ], [n_max+1, n_max+1] )
contains

  impure elemental function rand_int_in_range(lo,hi) result(num)
    use set_precision, only : dp
    integer, intent(in) :: lo, hi
    integer             :: num
    real(dp) :: harvest
    call random_number(harvest)
    num = nint( harvest*real(hi-lo,dp) + real(lo,dp) )
  end function rand_int_in_range

  pure function get_deriv_coef_v1(n_dim,exponents,diff_order) result(dcoef)
    integer,               intent(in) :: n_dim
    integer, dimension(:), intent(in) :: exponents, diff_order
    integer                           :: dcoef
    integer :: d, i
    dcoef = 1
    do d = 1,n_dim
      if ( exponents(d) > n_max .or. exponents(d)<diff_order(d) ) then
        dcoef = 0
        return
      end if
      do i = exponents(d),exponents(d)-diff_order(d)+1,-1
        dcoef = dcoef * i
      end do
    end do
  end function get_deriv_coef_v1

  pure function get_deriv_coef_v2(n_dim,exponents,diff_order) result(dcoef)
    integer,               intent(in) :: n_dim
    integer, dimension(:), intent(in) :: exponents, diff_order
    integer                           :: dcoef
    integer :: d
    ! dcoef = 1
    ! do d = 1,n_dim
    !   if ( exponents(d) > n_max .or. exponents(d)<diff_order(d) ) then
    !     dcoef = 0
    !     return
    !   end if
    !   dcoef = dcoef * fall_(exponents(d),diff_order(d))
    ! end do
    dcoef = 0
    if ( any( (exponents > n_max).or.(exponents<diff_order) ) ) return
    ! dcoef = 1
    ! do d = 1,n_dim
    !   dcoef = dcoef * fall_(exponents(d),diff_order(d))
    ! end do
    dcoef = product( fall_(exponents,diff_order) )
  end function get_deriv_coef_v2

  pure function nchoosek( n, k ) result( c )
    integer, intent(in) :: n, k
    integer             :: c
    integer :: i
    c = 0
    if (k>n) return

    c = 1
    do i = 1, min(n-k,k)
      c = c * ( n - (i-1) )
      c = c / i
    end do
  end function nchoosek

  pure subroutine get_exponents(n_dim,degree,n_terms,exponents,idx,diff_idx)
    use index_conversion, only : global2local
    integer, intent(in) :: n_dim, degree, n_terms
    integer, dimension(n_dim,n_terms), intent(out) :: exponents
    integer, dimension(0:degree),      intent(out) :: idx
    integer, dimension(n_dim,n_terms), optional, intent(out) :: diff_idx
    integer :: curr_total_degree, i, j, cnt, N_full_terms
    integer, dimension(n_dim) :: tmp_exp, nsub
    cnt = 0
    do curr_total_degree = 0,degree
      ! idx(curr_total_degree+1) = cnt + 1
      nSub = curr_total_degree + 1
      N_full_terms = (curr_total_degree+1) ** n_dim
      do j = 0,N_full_terms-1
        tmp_exp = global2local(j+1,nsub)-1
        if ( sum(tmp_exp) == curr_total_degree ) then
          cnt = cnt + 1
          exponents(:,cnt) = tmp_exp
        end if
      end do
      idx(curr_total_degree) = cnt
      tmp_exp = 0
    end do

    ! determine corresponding gradient terms of a given term
    if (present(diff_idx)) then
      diff_idx = -1 ! last terms (idx(degree)+1:idx(degree+1)) are not defined
      if ( degree==0) return
      do j = 1,idx(degree-1)
        tmp_exp = exponents(:,j)
        curr_total_degree = sum(tmp_exp)
        cnt = 0
        do i = idx(curr_total_degree)+1,idx(curr_total_degree+1)
          if ( sum( abs( exponents(:,i) - tmp_exp ) )==1 ) then
            cnt = cnt + 1
            diff_idx(cnt,j) = i
          end if
        end do
      end do
    end if
  end subroutine get_exponents

  pure function mono_rank_grlex( m, x ) result(r)
    integer, intent(in)               :: m
    integer, intent(in), dimension(:) :: x
    integer                           :: r
    integer, dimension(m-1) :: xs
    integer :: i, j
    integer :: ks, n, nm, ns, tim1
    r = -1
    if ( m < 1 ) return
    if ( any(x(1:m) < 0) ) return
    nm = sum(x(1:m))
    ns = nm + m - 1
    ks = m - 1
    xs(1) = x(1) + 1
    do i = 2, ks
      xs(i) = xs(i-1) + x(i) + 1
    end do
    ! xs = [(sum(x(1:i)+1), i = 1,m-1)]
    r = 1
    do i = 1, ks
      if ( i == 1 ) then
        tim1 = 0
      else
        tim1 = xs(i-1)
      end if
      if ( tim1 + 1 <= xs(i) - 1 ) then
        do j = tim1 + 1, xs(i) - 1
          r = r + nchoosek( ns - j, ks - i )
        end do
      end if
    end do
    do n = 0, nm - 1
      r = r + nchoosek( n + m - 1, n )
    end do
  end function mono_rank_grlex

  pure function mono_rank_grlex2( n_dim, e ) result(r)
    integer, intent(in)               :: n_dim
    integer, intent(in), dimension(:) :: e
    integer                           :: r
    integer, dimension(n_dim-1) :: es
    integer :: i, j
    integer :: ns, ks, total_degree, esim1

    r = -1
    if ( n_dim < 1 ) return
    if ( any( e(1:n_dim) < 0 ) ) return
    r = 1
    total_degree = sum( e(1:n_dim) )
    ! set of elements 
    ns = n_dim + total_degree - 1
    ks = n_dim - 1
    es = [(sum(e(1:i)+1), i = 1,n_dim-1)]
    r = 1
    do j = 1, es(1) - 1
      r = r + nchoosek( ns - j, ks - 1 )
    end do
    do i = 2, ks
      if ( es(i-1) + 1 <= es(i) - 1 ) then
        do j = es(i-1) + 1, es(i) - 1
          r = r + nchoosek( ns - j, ks - i )
        end do
      end if
    end do
    do i = 0, total_degree - 1
      r = r + nchoosek( i + n_dim - 1, i )
    end do
  end function mono_rank_grlex2

  pure function mono_rank_grlex3( n_dim, e ) result(r)
    integer, intent(in)               :: n_dim
    integer, intent(in), dimension(:) :: e
    integer                           :: r
    integer, dimension(n_dim-1) :: es
    integer :: i, j
    integer :: ns, ks, total_degree, esim1

    r = -1
    if ( n_dim < 1 ) return
    if ( any( e(1:n_dim) < 0 ) ) return
    r = 1
    total_degree = sum( e(1:n_dim) )
    ! set of elements 
    ns = n_dim + total_degree - 1
    ks = n_dim - 1
    ! es = [(sum(e(1:i)+1), i = 1,n_dim-1)]
    es = [(sum(e(n_dim:n_dim+1-i:-1)+1), i = 1,n_dim-1)]
    r = 1
    do j = 1, es(1) - 1
      r = r + nchoosek( ns - j, ks - 1 )
    end do
    do i = 2, ks
      if ( es(i-1) + 1 <= es(i) - 1 ) then
        do j = es(i-1) + 1, es(i) - 1
          r = r + nchoosek( ns - j, ks - i )
        end do
      end if
    end do
    do i = 0, total_degree - 1
      r = r + nchoosek( i + n_dim - 1, i )
    end do
  end function mono_rank_grlex3

  pure function mono_unrank_grlex( m, r ) result(x)
    integer,               intent(in)  :: m
    integer,               intent(in)  :: r
    integer, dimension(m)              :: x

  integer :: i, j
  integer :: ks, nksub, nm, ns, rtmp, r1, r2
  integer, dimension(m-1) :: xs

  x = -1
  if ( m < 1 ) return !  Ensure that m >= 1
  if ( r < 1 ) return !  Ensure that r >= 1

!  special case for m == 1
  if ( m == 1 ) then
    x(1) = r - 1
    return
  end if
!  Determine the appropriate value of NM.
!  Do this by adding up the number of compositions of sum 0, 1, 2, 
!  ..., without exceeding RANK.  Moreover, RANK - this sum essentially
!  gives you the rank of the composition within the set of compositions
!  of sum NM.  And that's the number you need in order to do the
!  unranking.
!
  r1 = 1
  nm = -1
  do
    nm = nm + 1
    rtmp = nchoosek( nm + m - 1, nm )
    if ( r < r1 + rtmp ) exit
    r1 = r1 + rtmp
  end do

  r2 = r - r1

  ks = m - 1
  ns = nm + m - 1

  nksub = nchoosek( ns, ks )

  j = 1

  do i = 1, ks
    rtmp = nchoosek( ns - j, ks - i )
    do while ( rtmp <= r2 .and. rtmp > 0 )
      r2 = r2 - rtmp
      j = j + 1
      rtmp = nchoosek( ns - j, ks - i )
    end do
    xs(i) = j
    j = j + 1
  end do

  x(1) = xs(1) - 1
  do i = 2, m - 1
    x(i) = xs(i) - xs(i-1) - 1
  end do
  x(m) = ns - xs(ks)
end function mono_unrank_grlex

pure function mono_unrank_grlex2( n_dim, r ) result(e)
    integer,               intent(in)  :: n_dim
    integer,               intent(in)  :: r
    integer, dimension(n_dim)              :: e

  integer :: i, j
  integer :: ks, nksub, nm, ns, rtmp, r1, r2
  integer, dimension(n_dim-1) :: es

  e = -1
  if ( n_dim < 1 ) return !  Ensure that n_dim >= 1
  if ( r < 1 ) return     !  Ensure that r >= 1

!  special case for n_dim == 1
  e = r - 1
  if ( n_dim == 1 ) return
  r1 =  1
  nm = -1
  do
    nm = nm + 1
    rtmp = nchoosek( nm + n_dim - 1, nm )
    if ( r < r1 + rtmp ) exit
    r1 = r1 + rtmp
  end do
  r2 = r - r1
  ks = n_dim - 1
  ns = nm + n_dim - 1
  nksub = nchoosek( ns, ks )
  j = 1
  do i = 1, ks
    rtmp = nchoosek( ns - j, ks - i )
    do while ( rtmp <= r2 .and. rtmp > 0 )
      r2 = r2 - rtmp
      j = j + 1
      rtmp = nchoosek( ns - j, ks - i )
    end do
    es(i) = j
    j = j + 1
  end do
  e(1) = es(1) - 1
  do i = 2, n_dim - 1
    e(i) = es(i) - es(i-1) - 1
  end do
  e(n_dim) = ns - es(ks)
end function mono_unrank_grlex2

pure function mono_unrank_grlex3( n_dim, r ) result(e)
    integer,               intent(in)  :: n_dim
    integer,               intent(in)  :: r
    integer, dimension(n_dim)              :: e
    integer :: i, j
    integer :: ks, nksub, nm, ns, rtmp, r1, r2
    integer, dimension(n_dim-1) :: es

  e = -1
  if ( n_dim < 1 ) return !  Ensure that n_dim >= 1
  if ( r < 1 ) return     !  Ensure that r >= 1

!  special case for n_dim == 1
  e = r - 1
  if ( n_dim == 1 ) return
  r1 =  1
  nm = -1
  do
    nm = nm + 1
    rtmp = nchoosek( nm + n_dim - 1, nm )
    if ( r < r1 + rtmp ) exit
    r1 = r1 + rtmp
  end do
  r2 = r - r1
  ks = n_dim - 1
  ns = nm + n_dim - 1
  nksub = nchoosek( ns, ks )
  j = 1
  do i = 1, ks
    rtmp = nchoosek( ns - j, ks - i )
    do while ( rtmp <= r2 .and. rtmp > 0 )
      r2 = r2 - rtmp
      j = j + 1
      rtmp = nchoosek( ns - j, ks - i )
    end do
    es(i) = j
    j = j + 1
  end do
  e(n_dim) = es(1) - 1
  do i = 2, n_dim - 1
    e(n_dim+1-i) = es(i) - es(i-1) - 1
  end do
  e(1) = ns - es(ks)
end function mono_unrank_grlex3

pure function get_diff_idx(n_dim,e,order) result(diff_idx)
  integer,               intent(in) :: n_dim
  integer, dimension(:), intent(in) :: e, order
  integer                           :: diff_idx
  diff_idx = mono_rank_grlex3(n_dim,e(1:n_dim)+order(1:n_dim))
end function get_diff_idx

end module combinatorics

! program main
!   use set_precision, only : dp
!   use combinatorics, only : nchoosek, get_exponents
!   use combinatorics, only : mono_rank_grlex, mono_rank_grlex2, mono_rank_grlex3, mono_unrank_grlex2, mono_unrank_grlex3
!   use combinatorics, only : get_diff_idx
!   use timer_derived_type, only : basic_timer_t
!   implicit none
!   type(basic_timer_t) :: timer
!   integer :: n_terms, n_dim, degree, i, d
!   integer, dimension(:,:), allocatable :: exponents, diff_idx, diff_idx2
!   integer, dimension(:),   allocatable :: idx, exp, order

!   n_dim  = 3
!   degree = 4
!   n_terms = nchoosek( n_dim + degree, degree )
!   allocate( exponents(n_dim,n_terms), idx(0:degree), diff_idx(n_dim,n_terms), diff_idx2(n_dim,n_terms), order(n_dim) )

!   call get_exponents(n_dim,degree,n_terms,exponents,idx,diff_idx=diff_idx)

!   ! do i = 1,n_terms
!   !   write(*,*) i, ':', exponents(:,i), mono_rank_grlex(n_dim,exponents(:,i)), &
!   !                                      mono_rank_grlex2(n_dim,exponents(:,i)), &
!   !                                      mono_rank_grlex3(n_dim,exponents(:,i)), &
!   !                                      '|', mono_unrank_grlex2(n_dim,mono_rank_grlex2(n_dim,exponents(:,i)))
!   !                                     !  '|', mono_unrank_grlex2(n_dim,mono_rank_grlex(n_dim,exponents(:,i)))
!   ! end do

!   ! do i = 1,n_terms
!   !   write(*,*) i, ':', exponents(:,i) - mono_unrank_grlex2(n_dim,mono_rank_grlex2(n_dim,exponents(:,i)))
!   !   write(*,*) i, ':', exponents(:,i) - mono_unrank_grlex3(n_dim,i)
!   ! end do

!   do i = 1,n_terms
!     do d = 1,n_dim
!       order = 0
!       order(d) = 1
!       diff_idx2(d,i) = get_diff_idx(n_dim,exponents(:,i),order)
!     end do
!   end do

!   do i = 1,n_terms
!     write(*,*) i, ':', diff_idx(:,i), '|', diff_idx2(:,i)
!   end do

!   ! do i = 1,n_terms
!   !   write(*,*) i, ':', exponents(:,i) - mono_unrank_grlex2(n_dim,mono_rank_grlex2(n_dim,exponents(:,i)))
!   !   write(*,*) i, ':', exponents(:,i) - mono_unrank_grlex3(n_dim,i)
!   ! end do


!   deallocate( exponents, idx, diff_idx )
! end program main

! program main
!   use set_precision, only : dp
!   use combinatorics, only : nchoosek, get_exponents
!   use combinatorics, only : mono_rank_grlex, mono_rank_grlex2, mono_rank_grlex3, mono_unrank_grlex2, mono_unrank_grlex3
!   use combinatorics, only : get_diff_idx
!   use combinatorics, only : get_deriv_coef_v1, get_deriv_coef_v2
!   use timer_derived_type, only : basic_timer_t
!   implicit none
!   type(basic_timer_t) :: timer
!   integer :: n_terms, n_dim, degree, i, j, dcoef1, dcoef2, max_coef
!   integer, dimension(:,:), allocatable :: exponents, diff_idx, diff_idx2
!   integer, dimension(:),   allocatable :: idx, exp, order

!   n_dim  = 7
!   degree = 8
!   n_terms = nchoosek( n_dim + degree, degree )
!   allocate( exponents(n_dim,n_terms), idx(0:degree), diff_idx(n_dim,n_terms), diff_idx2(n_dim,n_terms), order(n_dim) )

!   call get_exponents(n_dim,degree,n_terms,exponents,idx,diff_idx=diff_idx)

!   ! do i = 1,n_terms
!   do i = idx(degree-1)+1,idx(degree)
!     exp = exponents(:,i)
!     do j = idx(0)+1,idx(degree-1)
!       order = exponents(:,j)
!       dcoef1 = get_deriv_coef_v1(n_dim,exp,order)
!       dcoef2 = get_deriv_coef_v2(n_dim,exp,order)
!       max_coef = max(max_coef,dcoef1)
!       ! write(*,'(I2," (",2(I2,","),I2,"), (",2(I2,","),I2,") : ",I0)') i, exp, order, dcoef1
!       if (dcoef1 /= dcoef2 ) then
!         write(*,*) 'Error!'
!         stop
!       end if
!     end do
!   end do

!   write(*,*) max_coef
!   deallocate( exponents, idx, diff_idx )
! end program main


program main
  use set_precision, only : dp
  use set_constants, only : zero, one
  use combinatorics, only : rand_int_in_range
  use combinatorics, only : nchoosek, get_exponents
  use combinatorics, only : mono_rank_grlex, mono_rank_grlex2, mono_rank_grlex3, mono_unrank_grlex2, mono_unrank_grlex3
  use combinatorics, only : get_diff_idx
  use combinatorics, only : get_deriv_coef_v1, get_deriv_coef_v2
  use timer_derived_type, only : basic_timer_t
  implicit none
  type(basic_timer_t) :: timer
  integer :: n_terms, n_dim, degree, dcoef, junk
  integer :: i, j, k, m, n_iter, n_samples
  integer, dimension(:,:), allocatable :: exponents, diff_idx, rand_order
  integer, dimension(:),   allocatable :: idx, exp, order, rand_exp, rand_idx
  real(dp) :: t_avg, xnsamp
  n_samples = 10000
  n_iter = 1000000
  n_dim  = 3
  degree = 8
  n_terms = nchoosek( n_dim + degree, degree )
  allocate( exponents(n_dim,n_terms), idx(0:degree), diff_idx(n_dim,n_terms), order(n_dim) )
  allocate( rand_order(n_dim,n_iter), rand_idx(n_iter), rand_exp(n_iter) )
  call get_exponents(n_dim,degree,n_terms,exponents,idx,diff_idx=diff_idx)

  xnsamp = one / real(n_samples,dp)
  ! t_avg = zero
  ! do m = 1,n_samples
  !   call timer%tic()
  !   do k = 1,n_iter
  !     junk = 0
  !     do i = 1,n_terms
  !       exp = exponents(:,i)
  !       do j = idx(0)+1,idx(degree-1)
  !         order = exponents(:,j)
  !         dcoef = get_deriv_coef_v1(n_dim,exp,order)
  !         junk = junk + dcoef
  !       end do
  !     end do
  !   end do
  !   t_avg = t_avg + timer%toc()
  ! end do
  ! write(*,*) t_avg*xnsamp

  rand_idx   = rand_int_in_range(spread(1,1,n_iter),spread(n_terms,1,n_iter))
  rand_order = rand_int_in_range(spread(spread(0,1,n_dim),2,n_iter),exponents(:,rand_idx))
  t_avg = zero
  do m = 1,n_samples
    junk = 0
    call timer%tic()
    do k = 1,n_iter
      dcoef = get_deriv_coef_v2(n_dim,exponents(:,rand_idx(k)),rand_order(:,k))
      junk = junk + dcoef
    end do
    t_avg = t_avg + timer%toc()
  end do
  write(*,*) t_avg*xnsamp

  deallocate( rand_order, rand_exp )

  deallocate( exponents, idx, diff_idx )
end program main