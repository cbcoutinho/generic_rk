module misc
  use iso_fortran_env, only: wp => real64
  implicit none

contains

  subroutine mysub(n, t, y, dy)
    integer, intent(in) :: n
    real(wp), intent(in)  :: t
    real(wp), intent(in), dimension(n) :: y
    real(wp), intent(out), dimension(n) :: dy

    dy(1) = y(2)
    dy(2) = -y(1)

    return
  end subroutine mysub

end module misc
