module runge_kutta
  use iso_fortran_env, only: wp => real64
  implicit none

contains

  subroutine rk(t, dt, y, dysub)
    ! Dummy arguments
    real(wp), intent(in)                    :: t, dt
    real(wp), intent(inout),  dimension(2)  :: y
    real(wp),                 dimension(2)  :: dy

    interface
      subroutine dysub(t, y, dy)
        import wp
        real(wp), intent(in)                  :: t
        real(wp), intent(in),   dimension(2)  :: y
        real(wp), intent(out),  dimension(2)  :: dy
      end subroutine
    end interface

    ! Local arguments
    integer                               :: ii, jj, n, m
    real(wp)                              :: dummy_t
    real(wp), dimension(size(y))          :: dummy_y
    real(wp), dimension(:),   allocatable :: b, c
    real(wp), dimension(:,:), allocatable :: a, k


    call classic_rk4(a, b, c)

    n = size(b)
    m = size(y)
    ! print *, c
    ! print *, b
    !
    ! do ii = 1, n
    !   print *, a(ii, :)
    ! end do

    allocate(k(n, m))
    k = 0._wp

    call dysub(t, y, dy)
    k(1, :) = dy

    ! print *, t
    ! print *, y
    ! print *, dy
    ! print *,

    dummy_y = y

    do ii = 2, N

      dummy_t = t + dt*c(ii)
      do jj = 1, ii-1
        dummy_y = dummy_y + dt*a(ii,jj)*k(jj,:)
      end do

      call dysub(dummy_t, dummy_y, dy)
      k(ii, :) = dy

    end do

    y = y + [( dt*dot_product(b,k(:,jj)), jj = 1, m )]

    deallocate(a, b, c, k)

    return
  end subroutine rk

  subroutine classic_rk4(a, b, c)
    real(wp), intent(out), dimension(:),    allocatable :: c, b
    real(wp), intent(out), dimension(:,:),  allocatable :: a

    allocate(a(4,4), b(4), c(4))

    c = [0._wp, 0.5_wp, 0.5_wp, 1._wp]
    b = [1._wp/6._wp, 1._wp/3._wp, 1._wp/3._wp, 1._wp/6._wp ]

    a = 0._wp
    a(2,1) = 0.5_wp
    a(3,2) = 0.5_wp
    a(4,3) = 1._wp

    return
  end subroutine classic_rk4

  subroutine three_eighths_rk4(a, b, c)
    real(wp), intent(out), dimension(:),    allocatable :: c, b
    real(wp), intent(out), dimension(:,:),  allocatable :: a

    allocate(a(4,4), b(4), c(4))

    c = [0._wp, 1._wp/3._wp, 2._wp/6._wp, 1._wp]
    b = [1._wp/8._wp, 3._wp/8._wp, 3._wp/8._wp, 1._wp/8._wp ]

    a = 0._wp
    a(2,1) = 0.5_wp
    a(3,1:2) = [-1._wp/3._wp, 1._wp]
    a(4,1:3) = [1._wp, -1._wp, 1._wp]

    return
  end subroutine three_eighths_rk4

  subroutine midpoint(a, b, c)
    real(wp)                                            :: alpha
    real(wp), intent(out), dimension(:),    allocatable :: c, b
    real(wp), intent(out), dimension(:,:),  allocatable :: a

    allocate(a(2,2), b(2), c(2))

    alpha = 0.5_wp
    call two_stage(a, b, c, alpha)

    return
  end subroutine midpoint

  subroutine heun(a, b, c)
    real(wp)                                            :: alpha
    real(wp), intent(out), dimension(:),    allocatable :: c, b
    real(wp), intent(out), dimension(:,:),  allocatable :: a

    allocate(a(2,2), b(2), c(2))

    alpha = 1._wp
    call two_stage(a, b, c, alpha)

    return
  end subroutine heun

  subroutine ralston(a, b, c)
    real(wp)                                            :: alpha
    real(wp), intent(out), dimension(:),    allocatable :: c, b
    real(wp), intent(out), dimension(:,:),  allocatable :: a

    allocate(a(2,2), b(2), c(2))

    alpha = 2._wp/3._wp
    call two_stage(a, b, c, alpha)

    return
  end subroutine ralston

  subroutine two_stage(a, b, c, alpha)
    real(wp), intent(in) :: alpha
    real(wp), intent(out), dimension(2) :: c, b
    real(wp), intent(out), dimension(2,2) :: a

    c = [0._wp, alpha]
    b = [(1._wp - (1._wp/(2._wp*alpha))), (1._wp/(2._wp*alpha))]

    a = 0._wp
    a(2,1) = alpha

    return
  end subroutine two_stage

end module runge_kutta
