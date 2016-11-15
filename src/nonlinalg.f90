module nonlinalg
  use iso_fortran_env, only: wp => real64
  use linalg, only: linsolve_quick
  implicit none

  ! real(wp), parameter :: pi = 4._wp * datan(1._wp)

  interface
    function fun_interf(n, x) result(y)
      import wp
      integer,  intent(in)                  :: n
      real(wp), intent(in), dimension(n)    :: x
      real(wp),             dimension(n)    :: y
    end function
  end interface

contains

  function nonlinsolve(myfun, n, x0) result(x)
    integer,  intent(in)                  :: n
    real(wp), intent(in), dimension(n)    :: x0
    procedure(fun_interf)                 :: myfun

    integer                   :: ii
    real(wp), dimension(n)    :: x, z, b
    real(wp), dimension(n,n)  :: jac
    real(wp), parameter       :: eps = epsilon(1e0)

    x = x0

    do ii = 1, 100

      b = myfun(n,x)
      if ( ii == 1 .or. modulo(ii, 1) == 0 ) jac = jacobian(myfun, n, x)

      call linsolve_quick(n, jac, 1, -b, z)

      ! print '(8(f13.6))', x, norm2(z), norm2(myfun(n, x))

      x = x + z

      ! if ( norm2(myfun(n, x)) < eps ) exit
      if ( norm2(z) < eps ) exit
    end do

  end function nonlinsolve

  function jacobian(fun, n, x) result(dydx)
    integer,  intent(in)                  :: n
    real(wp), intent(in), dimension(n)    :: x
    real(wp),             dimension(n,n)  :: dydx
    procedure(fun_interf)                 :: fun

    integer                   :: jj
    real(wp), dimension(n)    :: xx, fun_plus, fun_minus
    real(wp), parameter       :: eps = epsilon(1e0)

    ! print*,
    ! print*, "Now Calculating jacobian"
    ! print*,

    dydx = 0._wp

    ! do jj = 1,n
    !   print*, dydx(jj,:)
    ! end do
    ! stop

    !$omp parallel do &
    !$omp& private (jj, xx, fun_plus, fun_minus) shared(dydx, x, n)
    do jj = 1,n
      xx = x
      xx(jj) = xx(jj) + eps
      fun_plus = fun(n, xx)

      xx = x
      xx(jj) = xx(jj) - eps
      fun_minus = fun(n, xx)

      dydx(:,jj) = (fun_plus - fun_minus ) / (2._wp * eps)
    end do
    !$omp end parallel do

    ! do jj = 1,n
    !   print*, dydx(jj,:)
    ! end do
    ! stop

    return
  end function jacobian

end module nonlinalg
