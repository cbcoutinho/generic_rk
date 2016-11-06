program main
  use iso_fortran_env, only: wp => real64
  use misc, only: mysub
  use runge_kutta, only: rk
  use lib_array, only: linspace
  use lib_constants, only: pi => pi_dp
  implicit none

  integer                         :: ii, ios
  integer, parameter              :: N = 50, num_eq = 2
  real(wp)                        :: dt
  real(wp), dimension(N)          :: t
  real(wp), dimension(num_eq)     :: y0
  real(wp), dimension(N,num_eq)   :: y

  call linspace(0._wp, 4._wp*pi, t)
  dt = t(2) - t(1)
  y(1,:) = [0._wp, 1._wp]

  do ii = 1, N-1
    y0 = y(ii, :)
    call rk(num_eq, t(ii), dt, y0, mysub)
    y(ii+1, :) = y0
  end do

  open(unit=21, file='data.out', iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file 21"
  do ii = 1, N
    write(21,*) t(ii), y(ii, :)
  end do
  close(unit=21, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 21"



end program main
