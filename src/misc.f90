module misc
  use iso_fortran_env, only: wp => real64
  use lib_constants, only: pi => pi_dp
  implicit none

contains

  subroutine mysub(n, t, y, dy)
    integer, intent(in)                 :: n
    real(wp), intent(in)                :: t
    real(wp), intent(in), dimension(n)  :: y
    real(wp), intent(out), dimension(n) :: dy

    ! call simple_trig(n, t, y, dy)
    ! call simple_ode(n, t, y, dy)
    ! call vanderpol(n, t, y, dy)
    call lorenz(n, t, y, dy)


    return
  end subroutine mysub

  subroutine simple_trig(n, t, y, dy)
    integer, intent(in)                 :: n
    real(wp), intent(in)                :: t
    real(wp), intent(in), dimension(n)  :: y
    real(wp), intent(out), dimension(n) :: dy

    dy(1) = dcos(t)
    dy(2) = -dsin(t)

    return
  end subroutine simple_trig

  subroutine simple_ode(n, t, y, dy)
    integer, intent(in)                 :: n
    real(wp), intent(in)                :: t
    real(wp), intent(in), dimension(n)  :: y
    real(wp), intent(out), dimension(n) :: dy

    dy(1) = y(2)
    dy(2) = -y(1)

    return
  end subroutine simple_ode

  subroutine vanderpol(n, t, y, dy)
    integer, intent(in)                 :: n
    real(wp), intent(in)                :: t
    real(wp), intent(in), dimension(n)  :: y
    real(wp), intent(out), dimension(n) :: dy

    real(wp), parameter                 :: mu = 50.53_wp
    real(wp), parameter                 :: A = 3.2_wp
    real(wp), parameter                 :: omega = 2._wp*pi/10_wp

    dy(1) = y(2)
    dy(2) = mu * (1-y(1)*y(1)) * y(2) - y(1) + A*dsin(t*omega)

    return
  end subroutine vanderpol

  subroutine lorenz(n, t, y, dy)
    integer, intent(in)                 :: n
    real(wp), intent(in)                :: t
    real(wp), intent(in), dimension(n)  :: y
    real(wp), intent(out), dimension(n) :: dy

    real(wp), parameter                 :: sigma = 10._wp
    real(wp), parameter                 :: rho = 28._wp
    real(wp), parameter                 :: beta = 8._wp/3._wp

    dy(1) = sigma * (y(2)-y(1))
    dy(2) = y(1) * (rho-y(3)) - y(2)
    dy(3) = y(1)*y(2) - beta*y(3)

    return
  end subroutine lorenz

end module misc
