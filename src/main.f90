program main
  use iso_fortran_env, only: wp => real64
  use misc, only: myfun
  use runge_kutta, only: midpoint, heun, ralston
  implicit none

  integer :: ii

  ! real(wp) :: t, y, dy
  ! t = 0_wp
  ! y = 0_wp
  ! dy = myfun(t, y)
  ! print *, t, y, dy

  integer, parameter :: N = 2
  real(wp), dimension(N) :: c, b
  real(wp), dimension(N,N) :: a

  call ralston(N, a, b, c)

  print *, c
  print *, b

  do ii = 1, N
    print *, a(ii, :)
  end do

end program main
