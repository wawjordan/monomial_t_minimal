module set_precision
  use iso_fortran_env, only : real64, int32, int64
  implicit none
  private
  public :: dp, i4, i8
  integer, parameter :: dp  = real64
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
  public :: nchoosek
  public :: get_exponents
contains

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

end module combinatorics

module monomial_basis_derived_type
  implicit none
  private
  public :: monomial_basis_t
  type :: monomial_basis_t
    private
    integer, public :: total_degree
    integer, public :: n_dim
    integer, public :: n_terms
    integer, public, allocatable, dimension(:)   :: idx
    integer, public, allocatable, dimension(:,:) :: exponents
    integer, public, allocatable, dimension(:,:) :: diff_idx
  contains
    private
    procedure, public, pass :: eval_m  => evaluate_monomial
    procedure, public, pass :: deval_m => evaluate_monomial_derivative
    procedure, public, pass :: destroy => destroy_monomial_basis_t
    procedure, public, pass :: check_gradient_indexing
  end type monomial_basis_t

  interface monomial_basis_t
    procedure constructor
  end interface monomial_basis_t

contains

  pure function constructor( total_degree, n_dim ) result(this)
    use combinatorics, only : nchoosek, get_exponents
    integer, intent(in) :: total_degree, n_dim
    type(monomial_basis_t) :: this

    call this%destroy()

    this%total_degree  = total_degree
    this%n_dim   = n_dim
    this%n_terms = nchoosek( n_dim + total_degree, total_degree )
    allocate( this%exponents( this%n_dim, this%n_terms ) )
    allocate( this%idx(0:this%total_degree) )
    allocate( this%diff_idx( this%n_dim, this%n_terms ) )
    call get_exponents( this%n_dim, this%total_degree, this%n_terms,           &
                        this%exponents, this%idx, diff_idx=this%diff_idx )
  end function constructor
  
  pure subroutine destroy_monomial_basis_t(this)
    class(monomial_basis_t), intent(inout) :: this
    if ( allocated(this%exponents) ) deallocate( this%exponents )
    if ( allocated(this%idx) )       deallocate( this%idx )
    if ( allocated(this%diff_idx) )  deallocate( this%diff_idx )
  end subroutine destroy_monomial_basis_t

  pure subroutine evaluate_monomial(this,term,x,val,coef)
    use set_precision, only : dp
    use set_constants, only : one
    class(monomial_basis_t), intent(in)  :: this
    integer,                 intent(in)  :: term
    real(dp), dimension(:),  intent(in)  :: x
    real(dp),                        intent(out) :: val
    integer,                         intent(out) :: coef
    integer :: d, i
    val  = one ! << x^[a] >>
    coef = 1   ! << [a]! >>
    do d = 1,this%n_dim
      do i = this%exponents(d,term),1,-1
        val  = val * x(d)
        coef = coef * i
      end do
    end do
  end subroutine evaluate_monomial

  pure subroutine evaluate_monomial_derivative( this, term, x, order,          &
                                                dval, dcoef, coef )
    use set_precision, only : dp
    use set_constants, only : zero, one
    class(monomial_basis_t),         intent(in)  :: this
    integer,                         intent(in)  :: term
    real(dp), dimension(:),          intent(in)  :: x
    integer,  dimension(:),          intent(in)  :: order
    real(dp),                        intent(out) :: dval
    integer,                         intent(out) :: dcoef, coef
    integer :: d, i
    
    dcoef = 1 ! D^[b](x^[a]) = ([a]!)/([a]-[b])! x^[a-b] => << ([a]!)/([a]-[b])! >>
    coef  = 1 ! << ([a]-[b])! >>
    dval  = zero ! << D^[b](x^[a]) >>
    if ( any( this%exponents(:,term)-order(1:this%n_dim) < 0 ) ) return

    dval  = one
    do d = 1,this%n_dim
      do i = this%exponents(d,term),this%exponents(d,term)-order(d)+1,-1
        dcoef = dcoef * i
      end do
      do i = this%exponents(d,term)-order(d),1,-1
        dval  = dval * x(d)
        coef = coef * i
      end do
    end do
  end subroutine evaluate_monomial_derivative

  subroutine check_gradient_indexing( this )
    use set_constants, only : max_text_line_length
    use string_stuff, only : write_integer_tuple
    class(monomial_basis_t), intent(in)  :: this
    integer :: d, term
    integer, dimension(this%n_dim) :: grad_idx
    character(max_text_line_length) :: tmp_string, out_string
    write(*,*) this%idx
    do d = this%total_degree-1,1,-1
      ! for each term with this total degree:
      do term = this%idx(d-1)+1,this%idx(d)
        grad_idx = this%diff_idx(:,term) ! indices to extract gradient information
        out_string=''
        call write_integer_tuple([d,term],tmp_string)
        write(out_string,'(A,A)') trim(out_string),trim(tmp_string)//' : '
        call write_integer_tuple(this%exponents(:,term),tmp_string)
        write(out_string,'(A,A)') trim(out_string),'['//trim(tmp_string)//'] : '
        call write_integer_tuple(grad_idx,tmp_string)
        write(out_string,'(A,A)') trim(out_string),'('//trim(tmp_string)//')'
        write(*,'(A)') trim(out_string)
      end do
    end do
    term = 1
    grad_idx = this%diff_idx(:,term) ! indices to extract gradient information
    out_string=''
    call write_integer_tuple([d,term],tmp_string)
    write(out_string,'(A,A)') trim(out_string),trim(tmp_string)//' : '
    call write_integer_tuple(this%exponents(:,term),tmp_string)
    write(out_string,'(A,A)') trim(out_string),'['//trim(tmp_string)//'] : '
    call write_integer_tuple(grad_idx,tmp_string)
    write(out_string,'(A,A)') trim(out_string),'('//trim(tmp_string)//')'
    write(*,'(A)') trim(out_string)

    ! constant term
  end subroutine check_gradient_indexing

end module monomial_basis_derived_type

program main
  use set_precision, only : dp
  use monomial_basis_derived_type, only : monomial_basis_t
  implicit none
  integer :: n_dim, rec_degree
  type(monomial_basis_t) :: p
  n_dim      = 2
  rec_degree = 4
  p = monomial_basis_t( rec_degree, n_dim )
  call p%check_gradient_indexing()
  call p%destroy()
end program main