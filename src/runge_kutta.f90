module runge_kutta
  use iso_fortran_env, only: wp => real64
  implicit none

contains

  subroutine midpoint(N, a, b, c)
    integer, intent(in) :: N
    real(wp) :: alpha
    real(wp), intent(out), dimension(N) :: c, b
    real(wp), intent(out), dimension(N,N) :: a

    alpha = 0.5_wp
    call two_stage(N, a, b, c, alpha)

  end subroutine midpoint

  subroutine heun(N, a, b, c)
    integer, intent(in) :: N
    real(wp) :: alpha
    real(wp), intent(out), dimension(N) :: c, b
    real(wp), intent(out), dimension(N,N) :: a

    alpha = 1._wp
    call two_stage(N, a, b, c, alpha)

  end subroutine heun

  subroutine ralston(N, a, b, c)
    integer, intent(in) :: N
    real(wp) :: alpha
    real(wp), intent(out), dimension(N) :: c, b
    real(wp), intent(out), dimension(N,N) :: a

    alpha = 2._wp/3._wp
    call two_stage(N, a, b, c, alpha)

  end subroutine ralston

  subroutine two_stage(N, a, b, c, alpha)
    integer, intent(in) :: N
    real(wp), intent(in) :: alpha
    real(wp), intent(out), dimension(N) :: c, b
    real(wp), intent(out), dimension(N,N) :: a

    c = [0._wp, alpha]
    b = [(1._wp - (1._wp/(2._wp*alpha))), (1._wp/(2._wp*alpha))]

    a = 0._wp
    a(2,1) = alpha

  end subroutine two_stage

end module runge_kutta
