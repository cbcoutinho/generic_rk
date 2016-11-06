program main
  use iso_fortran_env, only: wp => real64
  use misc, only: myfun
  implicit none

  integer :: ii
  real(wp) :: t, y, dy

  t = 0_wp
  y = 0_wp

  dy = myfun(t, y)

  print *, t, y, dy

end program main
