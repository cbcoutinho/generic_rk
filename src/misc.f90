module misc
  use iso_fortran_env, only: wp => real64
  implicit none

contains

  elemental function myfun(t, y) result(dy)
    real(wp), intent(in)  :: t, y
    real(wp)              :: dy

    dy = cos(t)

    return
  end function myfun

end module misc
