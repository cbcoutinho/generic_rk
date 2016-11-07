program main
  use iso_fortran_env, only: wp => real64
  use misc, only: mysub
  use runge_kutta, only: rk_wrapper
  use lib_array, only: linspace
  use lib_constants, only: pi => pi_dp
  implicit none

  integer                 :: ii, ios

  integer, parameter      :: n = 3
  real(wp)                :: t0, tend, dt
  ! real(wp), dimension(num_t)    :: t
  real(wp), dimension(n)  :: y0
  ! real(wp), dimension(num_t,n)  :: y



  open(unit=21, file='data.out', iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file 21"

  t0 = 0.0_wp
  dt = 0.1_wp
  ! tend = 10_wp*pi
  tend = 100._wp
  y0 = [0.1_wp, 0.0_wp, 0.0_wp]

  ! call linspace(t0, tend, t)
  ! dt = t(2) - t(1)
  ! do ii = 1, num_t-1
  !   y0 = y(ii, :)
  !   call rk_explicit(n, t(ii), dt, y0, mysub)
  !   y(ii+1, :) = y0
  ! end do
  ! do ii = 1, num_t
  !   write(21,*) t(ii), y(ii, :)
  ! end do

  write(21,*) t0, y0

  call rk_wrapper(n, t0, tend, dt, y0, mysub)






  close(unit=21, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 21"



end program main
