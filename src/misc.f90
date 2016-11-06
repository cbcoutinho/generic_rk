module misc
  use iso_fortran_env, only: wp => real64
  implicit none

contains

  subroutine mysub(N, t, y, dy)
    integer, intent(in) :: N
    real(wp), intent(in)  :: t
    real(wp), intent(in), dimension(N) :: y
    real(wp), intent(out), dimension(N) :: dy

    ! dy(1) = dcos(t)
    ! dy(2) = dsin(t)

    dy(1) = y(2)
    dy(2) = -y(1)

    return
  end subroutine mysub

end module misc
