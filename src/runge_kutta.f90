module runge_kutta
  use iso_fortran_env, only: wp => real64
  implicit none

contains

  subroutine rk_wrapper(n, t0, tend, dt, y0, dysub)
    ! Dummy arguments
    integer,  intent(in)                    :: n
    real(wp), intent(in)                    :: t0, tend, dt
    real(wp), intent(inout),  dimension(n)  :: y0
    real(wp),                 dimension(n)  :: dy

    interface
      subroutine dysub(n, t, y, dy)
        import wp
        integer,  intent(in)                  :: n
        real(wp), intent(in)                  :: t
        real(wp), intent(in),   dimension(n)  :: y
        real(wp), intent(out),  dimension(n)  :: dy
      end subroutine
    end interface

    ! Local arguments
    integer                               :: ii, jj, kk, m
    real(wp), parameter                   :: eps1 = 1d-8 ! sqrt(epsilon(1._wp))
    real(wp), parameter                   :: min_dt = 1d-10, max_dt = 1d-1
    real(wp)                              :: t
    real(wp), dimension(n)                :: y

    y = y0
    t = t0

    do while ( t <= tend )
      call rk_adaptive(n, t, dt, y, dysub)
      ! print*, t, y
      write(21, *) t, y
      ! stop
    end do

    return
  end subroutine rk_wrapper

  subroutine rk_adaptive(n, t0, tend, y, dysub)
    ! Dummy arguments
    integer,  intent(in)                    :: n
    real(wp), intent(in)                    :: tend
    real(wp), intent(inout)                 :: t0
    real(wp), intent(inout),  dimension(n)  :: y
    real(wp),                 dimension(n)  :: dy

    interface
      subroutine dysub(n, t, y, dy)
        import wp
        integer,  intent(in)                  :: n
        real(wp), intent(in)                  :: t
        real(wp), intent(in),   dimension(n)  :: y
        real(wp), intent(out),  dimension(n)  :: dy
      end subroutine
    end interface

    ! Local arguments
    integer                               :: ii, jj, kk, m
    real(wp), parameter                   :: eps1 = 1d-8 ! sqrt(epsilon(1._wp))
    real(wp), parameter                   :: min_dt = 1d-10, max_dt = 1d-1
    real(wp)                              :: t, dt, error
    real(wp), dimension(n)                :: dummy_y, y_star
    real(wp), dimension(:),   allocatable :: c, dummy_b
    real(wp), dimension(:,:), allocatable :: b, a, k
    logical                               :: switch

    dt  = tend / 10._wp
    ! dt  = (tend-t0) / 10._wp
    dt  = minval([dt, max_dt])
    t   = t0
    switch = .true.
    kk = 1

    call dormand_prince_rk45(a, b, c, m)
    allocate(dummy_b(m))

    ! Begin loop to solve for y until end of timestep (tend)
    do while ( t <= tend + t0 )

      ! Begin loop to reach convergence from a specified point in time
      error = 1._wp
      do while ( error >= eps1 )

        ! I`m recalculating k(m,n) in each of these rk calls. Maybe that's a
        ! performance issue
        dummy_y = y
        y_star  = y
        dummy_b = b(1,:)
        call rk_explicit(n, t, dt, dummy_y, dysub, a, dummy_b, c, m)
        dummy_b = b(2,:)
        call rk_explicit(n, t, dt, y_star, dysub, a, dummy_b, c, m)

        error = maxval(abs(dummy_y - y_star))

        dt = dt * 0.75_wp
      end do

      y = dummy_y
      ! print 120, t, dt, tend, y

      if ( t + dt >= tend + t0) then
        if ( switch ) then
          dt = tend+t0 - t
          t = t + dt
          switch = .false.
        else
          exit
        end if
      else
        t = t + dt
        ! dt = init_dt
        dt = dt * 2._wp
        dt = minval([dt, max_dt])
      end if

      ! print 120, t, dt, tend, y
      kk = kk+1

    end do

    t0 = t
    ! print*, kk
    deallocate(a, b, c, dummy_b)
    120 format (1x,5(e15.7))

    return
  end subroutine rk_adaptive

  subroutine rk_explicit(n, t, dt, y, dysub, a, b, c, m)
    ! Dummy arguments
    integer,  intent(in)                      :: n, m
    real(wp), intent(in)                      :: t, dt
    real(wp), intent(inout),  dimension(n)    :: y
    real(wp), intent(in),     dimension(m)    :: b, c
    real(wp), intent(in),     dimension(m,m)  :: a

    interface
      subroutine dysub(n, t, y, dy)
        import wp
        integer,  intent(in)                  :: n
        real(wp), intent(in)                  :: t
        real(wp), intent(in),   dimension(n)  :: y
        real(wp), intent(out),  dimension(n)  :: dy
      end subroutine
    end interface

    ! Local arguments
    integer                   :: ii, jj
    real(wp)                  :: dummy_t
    real(wp), dimension(n)    :: dummy_y
    real(wp), dimension(m,n)  :: k
    real(wp), dimension(n)    :: dy

    k = 0._wp
    dummy_y = y

    do ii = 1, m

      dummy_t = t + dt*c(ii)
      do jj = 1, ii-1
        dummy_y = dummy_y + dt*a(ii,jj)*k(jj,:)
      end do

      call dysub(n, dummy_t, dummy_y, dy)
      k(ii, :) = dy

    end do

    y = y + [( dt*dot_product(b,k(:,jj)), jj = 1, n )]

    return
  end subroutine rk_explicit

  subroutine classic_rk4(a, b, c, m)
    integer,  intent(out)                               :: m
    real(wp), intent(out), dimension(:),    allocatable :: c, b
    real(wp), intent(out), dimension(:,:),  allocatable :: a

    m = 4
    allocate(a(m,m), b(m), c(m))

    c = [0._wp, 0.5_wp, 0.5_wp, 1._wp]
    b = [1._wp/6._wp, 1._wp/3._wp, 1._wp/3._wp, 1._wp/6._wp ]

    a = 0._wp
    a(2,1) = 0.5_wp
    a(3,2) = 0.5_wp
    a(4,3) = 1._wp

    return
  end subroutine classic_rk4

  subroutine three_eighths_rk4(a, b, c, m)
    integer,  intent(out)                               :: m
    real(wp), intent(out), dimension(:),    allocatable :: c, b
    real(wp), intent(out), dimension(:,:),  allocatable :: a

    m = 4
    allocate(a(m,m), b(m), c(m))

    c = [0._wp, 1._wp/3._wp, 2._wp/6._wp, 1._wp]
    b = [1._wp/8._wp, 3._wp/8._wp, 3._wp/8._wp, 1._wp/8._wp ]

    a = 0._wp
    a(2,1) = 0.5_wp
    a(3,1:2) = [-1._wp/3._wp, 1._wp]
    a(4,1:3) = [1._wp, -1._wp, 1._wp]

    return
  end subroutine three_eighths_rk4

  subroutine midpoint(a, b, c, m)
    integer,  intent(out)                               :: m
    real(wp)                                            :: alpha
    real(wp), intent(out), dimension(:),    allocatable :: c, b
    real(wp), intent(out), dimension(:,:),  allocatable :: a

    m = 2
    allocate(a(m,m), b(m), c(m))

    alpha = 0.5_wp
    call two_stage(a, b, c, m, alpha)

    return
  end subroutine midpoint

  subroutine heun(a, b, c, m)
    integer,  intent(out)                               :: m
    real(wp)                                            :: alpha
    real(wp), intent(out), dimension(:),    allocatable :: c, b
    real(wp), intent(out), dimension(:,:),  allocatable :: a

    m = 2
    allocate(a(m,m), b(m), c(m))

    alpha = 1._wp
    call two_stage(a, b, c, m, alpha)

    return
  end subroutine heun

  subroutine ralston(a, b, c, m)
    integer,  intent(out)                               :: m
    real(wp)                                            :: alpha
    real(wp), intent(out), dimension(:),    allocatable :: c, b
    real(wp), intent(out), dimension(:,:),  allocatable :: a

    m = 2
    allocate(a(m,m), b(m), c(m))

    alpha = 2._wp/3._wp
    call two_stage(a, b, c, m, alpha)

    return
  end subroutine ralston

  subroutine two_stage(a, b, c, m, alpha)
    integer,  intent(in)                  :: m
    real(wp), intent(in)                  :: alpha
    real(wp), intent(out), dimension(m)   :: c, b
    real(wp), intent(out), dimension(m,m) :: a

    c = [0._wp, alpha]
    b = [(1._wp - (1._wp/(2._wp*alpha))), (1._wp/(2._wp*alpha))]

    a = 0._wp
    a(2,1) = alpha

    return
  end subroutine two_stage

  subroutine heun_euler(a, b, c, m)
    integer,  intent(out)                               :: m
    real(wp)                                            :: alpha
    real(wp), intent(out), dimension(:),    allocatable :: c
    real(wp), intent(out), dimension(:,:),  allocatable :: a, b

    real(wp),              dimension(:),    allocatable :: dummy_b

    m = 2
    allocate(a(m,m), b(2,m), c(m), dummy_b(m))

    ! There is already a routine to calculate the heun coefficients - just
    ! reuse it for the first row of `b`

    alpha = 1._wp
    call two_stage(a, dummy_b, c, m, alpha)

    b(1,:) = dummy_b
    b(2,:) = [1._wp, 0._wp]

    deallocate(dummy_b)

    return
  end subroutine heun_euler

  subroutine fehlbery_rk12(a, b, c, m)
    integer,  intent(out)                               :: m
    real(wp), intent(out), dimension(:),    allocatable :: c
    real(wp), intent(out), dimension(:,:),  allocatable :: a, b

    m = 3
    allocate(a(m,m), b(2,m), c(m))

    c       = [0._wp, 0.5_wp, 1.0_wp]
    b(1,:)  = [1._wp/256._wp, 255._wp/256._wp, 0._wp]
    b(2,:)  = [1._wp/512._wp, 255._wp/256._wp, 1._wp/512._wp]

    a       = 0._wp
    a(2,1)  = 0.5_wp
    a(3,:)  = [1._wp/256._wp, 255._wp/256._wp, 0._wp]

    return
  end subroutine fehlbery_rk12

  subroutine bogacki_shampine(a, b, c, m)
    integer,  intent(out)                               :: m
    real(wp), intent(out), dimension(:),    allocatable :: c
    real(wp), intent(out), dimension(:,:),  allocatable :: a, b

    m = 4
    allocate(a(m,m), b(2,m), c(m))

    c         = [0._wp, 0.5_wp, 0.75_wp, 1.0_wp]
    b(1,:)    = [2._wp/9._wp, 1._wp/3._wp, 4._wp/9._wp, 0._wp]
    b(2,:)    = [7._wp/24._wp, 1._wp/4._wp, 1._wp/3._wp, 1._wp/8._wp]

    a         = 0._wp
    a(2,1)    = 0.5_wp
    a(3,2)    = 0.75_wp
    a(4,1:3)  = [2._wp/9._wp, 1._wp/3._wp, 4._wp/9._wp]

    return
  end subroutine bogacki_shampine

  subroutine fehlbery_rk45(a, b, c, m)
    integer,  intent(out)                               :: m
    real(wp), intent(out), dimension(:),    allocatable :: c
    real(wp), intent(out), dimension(:,:),  allocatable :: a, b

    m = 6
    allocate(a(m,m), b(2,m), c(m))

    c         = [0._wp, 0.25_wp, 3._wp/8._wp, 12._wp/13._wp, 1._wp, 0.5_wp]
    b(1,:)    = [16._wp/135._wp, 0._wp, 6656._wp/12825._wp, &
                  & 28561._wp/56430._wp, -9._wp/50._wp, 2._wp/55._wp]
    b(2,:)    = [25._wp/216._wp, 0._wp, 1408._wp/2565._wp, 2197._wp/4104._wp, &
                  & -1._wp/5._wp, 0._wp]

    a         = 0._wp
    a(2,1)    = 1._wp/4._wp
    a(3, 1:2) = [3._wp/32._wp, 9._wp/32._wp]
    a(4, 1:3) = [1932._wp/2197._wp, -7200._wp/2197._wp, 7296._wp/2197._wp]
    a(5, 1:4) = [439._wp/216._wp, -8._wp, 3680._wp/513._wp, -845._wp/4104._wp]
    a(6, 1:5) = [-8._wp/27._wp, 2._wp, -3544._wp/2565._wp, 1859._wp/4104._wp, &
                  & -11._wp/40._wp]

    return
  end subroutine fehlbery_rk45

  subroutine cash_karp_rk45(a, b, c, m)
    integer,  intent(out)                               :: m
    real(wp), intent(out), dimension(:),    allocatable :: c
    real(wp), intent(out), dimension(:,:),  allocatable :: a, b

    m = 6
    allocate(a(m,m), b(2,m), c(m))

    c         = [0._wp, 1._wp/5._wp, 3._wp/10._wp, 3._wp/5._wp, 1._wp, &
                  & 7._wp/8._wp]
    b(1,:)    = [37._wp/378._wp, 0._wp, 250._wp/621._wp, 125._wp/594._wp, &
                  & 0._wp, 512._wp/1771._wp]
    b(2,:)    = [2825._wp/27648._wp, 0._wp, 18575._wp/48384._wp, &
                  & 13525._wp/55296._wp, 277._wp/14336._wp, 1._wp/4._wp]

    a         = 0._wp
    a(2,1)    = 1._wp/5._wp
    a(3, 1:2) = [3._wp/40._wp, 9._wp/40._wp]
    a(4, 1:3) = [3._wp/10._wp, -9._wp/10._wp, 6._wp/5._wp]
    a(5, 1:4) = [-11._wp/54._wp, 5._wp/2._wp, -70._wp/27._wp, 35._wp/27._wp]
    a(6, 1:5) = [1631._wp/55296._wp, 175._wp/512._wp, 575._wp/13824._wp, &
                  & 44275._wp/110592._wp, 253._wp/4096._wp]

    return
  end subroutine cash_karp_rk45

  subroutine dormand_prince_rk45(a, b, c, m)
    integer,  intent(out)                               :: m
    real(wp), intent(out), dimension(:),    allocatable :: c
    real(wp), intent(out), dimension(:,:),  allocatable :: a, b

    m = 7
    allocate(a(m,m), b(2,m), c(m))

    c         = [0._wp, 1._wp/5._wp, 3._wp/10._wp, 4._wp/5._wp, 8._wp/9._wp, &
                  & 1._wp, 1._wp]
    b(1,:)    = [35._wp/384._wp, 0._wp, 500._wp/1113._wp, 125._wp/192._wp, &
                  & -2187._wp/6784._wp, 11._wp/84._wp, 0._wp]
    b(2,:)    = [5179._wp/57600._wp, 0._wp, 7571._wp/16695._wp, &
                  & 393._wp/640._wp, -92097._wp/339200._wp, 187._wp/2100._wp, &
                  & 1._wp/40._wp]

    a         = 0._wp
    a(2,1)    = 1._wp/5._wp
    a(3, 1:2) = [3._wp/40._wp, 9._wp/40._wp]
    a(4, 1:3) = [44._wp/45._wp, -56._wp/15._wp, 32._wp/9._wp]
    a(5, 1:4) = [19372._wp/6561._wp, -25360._wp/2187._wp, 64448._wp/6561._wp, &
                  & -212._wp/729._wp]
    a(6, 1:5) = [9017._wp/3168._wp, -355._wp/33._wp, 46732._wp/5247._wp, &
                  & 49._wp/176._wp, -5103._wp/18656._wp]
    a(7, 1:6) = [35._wp/384._wp, 0._wp, 500._wp/1113._wp, 125._wp/192._wp, &
                  & -2187._wp/6784._wp, 11._wp/84._wp]

    return
  end subroutine dormand_prince_rk45

end module runge_kutta
