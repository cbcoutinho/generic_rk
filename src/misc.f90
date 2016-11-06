module misc
  use iso_fortran_env, only: wp => real64
  implicit none

contains

  subroutine mysub(t, y, dy)
    real(wp), intent(in)  :: t
    real(wp), intent(in), dimension(2) :: y
    real(wp), intent(out), dimension(2) :: dy

    ! dy(1) = dcos(t)
    ! dy(2) = dsin(t)

    dy(1) = y(2)
    dy(2) = -y(1)

    return
  end subroutine mysub

end module misc
