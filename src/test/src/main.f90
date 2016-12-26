program main
  use iso_fortran_env, only: wp => real64
  use misc, only: mysub, pi
  use runge_kutta, only: rk_wrapper
  use lib_array, only: linspace
  implicit none

  integer                       :: ii, ios
  integer, parameter            :: n = 3, num_t = 10
  real(wp)                      :: t0, tend
  real(wp), dimension(num_t)    :: t
  real(wp), dimension(n)        :: y0
  real(wp), dimension(num_t,n)  :: y


  ! Set time vector `t`
  t0   = 0.0_wp
  tend = 100._wp
  ! tend = 4_wp*pi
  call linspace(t0, tend, t)

  ! Set y vector `y` using y0
  ! y0 = [0._wp, 0.1_wp]
  ! y0 = [2._wp, 0._wp]
  ! y0 = [ 0._wp ]
  ! y0 = [0.1_wp, 10.0_wp, 10.0_wp]
  ! y0      = [1._wp, 0._wp, -1._wp]
  y0      = [1._wp, 0._wp, 0._wp]
  y(1,:)  = y0



  open(unit=22, file='raw.out', iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file 22"

  write(22,*) t0, 0._wp, y0
  call rk_wrapper(n, num_t, t, y, mysub)

  close(unit=22, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 22"







  open(unit=21, file='data.out', iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file 21"

  do ii = 1, num_t
      write(21,*) t(ii), y(ii, :)
  end do

  close(unit=21, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 21"


end program main
