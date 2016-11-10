module runge_kutta
  use iso_fortran_env, only: wp => real64
  use rk_constants, only: classic_rk4, &
                        & three_eighths_rk4, &
                        & midpoint, &
                        & heun, &
                        & ralston
  use rk_constants, only: heun_euler, &
                        & fehlbery_rk12, &
                        & bogacki_shampine, &
                        & fehlbery_rk45, &
                        & cash_karp_rk45, &
                        & dormand_prince_rk45, &
                        & fehlbery_rk78
  implicit none

  interface
    subroutine sub_interface(n, t, y, dy)
      import wp
      integer,  intent(in)                  :: n
      real(wp), intent(in)                  :: t
      real(wp), intent(in),   dimension(n)  :: y
      real(wp), intent(out),  dimension(n)  :: dy
    end subroutine
  end interface

contains

  subroutine rk_wrapper(n, num_t, t, y, dysub)
    ! Dummy arguments
    integer,  intent(in)                          :: n, num_t
    real(wp), intent(inout),  dimension(num_t)    :: t
    real(wp),                 dimension(num_t, n) :: y
    procedure(sub_interface)                      :: dysub

    ! Local arguments
    integer                                       :: ii, m, p
    real(wp)                                      :: dt
    real(wp), dimension(2)                        :: tspan
    real(wp), dimension(n)                        :: yy
    real(wp), dimension(:),   allocatable         :: c, b, bstar
    real(wp), dimension(:,:), allocatable         :: a

    ! Get the tableau coefficients associated with a runge kutta implementation
    ! call heun_euler(a, b, bstar, c, m, p)
    ! call fehlbery_rk12(a, b, bstar, c, m, p)
    ! call bogacki_shampine(a, b, bstar, c, m, p)
    ! call fehlbery_rk45(a, b, bstar, c, m, p)
    call dormand_prince_rk45(a, b, bstar, c, m, p)
    ! call fehlbery_rk78(a, b, bstar, c, m, p)

    ! call heun(a, b, c, m, p)
    ! call three_eighths_rk4(a, b, c, m, p)

    ! Initial guess for dt is just the difference between t(1:2), this way new
    ! versions of dt will be saved in each call to rk_adaptive
    !
    ! EDIT 10-11-2016: change dt to be the minimum value between t(2)-t(1) and
    ! 1d-2 because sometimes you only want a very sparse output of data, and
    ! that could overshoot a good initial guess. This should be smarter
    dt = minval([t(2) - t(1), 1d-1])

    do ii = 1, num_t-1
      yy = y(ii,:)
      tspan = t(ii:ii+1)

      call rk_adaptive(n, tspan, dt, yy, dysub, a, b, bstar, c, m, p)
      ! call rk_explicit(n, t(ii), t(ii+1)-t(ii), yy, dysub, a, b, c, m)

      y(ii+1, :) = yy

      ! if (ii == 1) stop
    end do

    if (allocated(a)) deallocate(a)
    if (allocated(b)) deallocate(b)
    if (allocated(bstar)) deallocate(bstar)
    if (allocated(c)) deallocate(c)

    return
  end subroutine rk_wrapper

  subroutine rk_adaptive(n, tspan, dt, y, dysub, a, b, bstar, c, m, p)
    ! Dummy arguments
    integer,  intent(in)                      :: n, m, p
    real(wp), intent(inout)                   :: dt
    real(wp), intent(in),     dimension(2)    :: tspan
    real(wp), intent(inout),  dimension(n)    :: y
    real(wp), intent(in),     dimension(m)    :: c, b, bstar
    real(wp), intent(in),     dimension(m,m)  :: a
    procedure(sub_interface)                  :: dysub

    ! Local arguments
    real(wp)                              :: t, error
    real(wp), parameter                   :: eps_abs = 1d-4 ! sqrt(epsilon(1._wp))
    real(wp), parameter                   :: eps_rel = 1d-4
    real(wp), dimension(n)                :: dummy_y, y_star

    integer                               :: eval, num_dt
    integer, parameter                    :: max_eval = 500
    real(wp)                              :: s, gamma, sc, maxy
    logical                               :: last

    num_dt = 1
    t = tspan(1)
    last = .false.

    ! Begin loop to solve for y until end of timestep (tend)
    do while ( t <= tspan(2) )

      ! Begin loop to reach convergence from a specified point in time
      do eval = 1, max_eval
        ! print*, 'eval     = ', eval
        ! print*, 'time     = ', t
        ! print*, 'y        = ', y
        ! print*, 'dt       = ', dt
        ! k(m,n) is recalculated in each of these rk loops. Maybe that's a
        ! performance issue
        dummy_y = y
        y_star  = y
        call rk_explicit(n, t, dt, dummy_y, dysub, a, b,      c, m)
        call rk_explicit(n, t, dt, y_star,  dysub, a, bstar,  c, m)

        ! print*, 'dummy_y  = ', dummy_y
        ! print*, 'y_star   = ', y_star

        ! error = maxval(abs(dummy_y - y_star))
        error = norm2(dummy_y - y_star)
        maxy = maxval([ abs(dummy_y), abs(y_star) ])

        gamma = 0.25_wp ** (1._wp / real(p,wp))
        ! gamma = 0.9_wp

        sc = eps_abs + maxy * eps_rel

        s = gamma * (sc/error) ** ( 1._wp/real(p,wp) )
        ! s = eps_abs*dt / &
        !     & (2._wp * (tspan(2) - tspan(1)) * error)
        ! print*, 'error    = ', error
        ! print*, 'maxy     = ', maxy
        ! print*, 'sc       = ', sc
        ! print*, 's        = ', s
        ! print*, 'sqrt(s)  = ', sqrt(s)
        ! print*,
        ! stop

        ! Adjust dt based on estimation of error
        if ( s >= 2._wp ) then
          t = t + dt
          y = y_star
          dt = dt * 2._wp
          ! print*, 'Adjusted dt with factor of 2'
          exit
        else if ( s < 2._wp .and. s >= 1._wp ) then
          t = t + dt
          y = y_star
          dt = dt * s
          ! print 121, 'Adjusted dt with factor of ', s
          exit
        else if ( s <= 1._wp ) then
          dt = dt * s
          ! print 121, 'Adjusted dt with factor of ', s
        end if

        if ( eval == max_eval ) then
          print*, "Warning, max eval iterations reached!"
        end if

      end do

      write(22,*) t, dt, y

      ! num_dt = num_dt + 1
      ! print 120, t, tspan, dt, y

      ! Make sure t does not overshoot tspan(2) (i.e. end of timestep)
      if ( last ) then
        exit
      end if

      if ( t + dt > tspan(2) ) then
        dt = tspan(2) - t
        last = .true.
      end if

    end do

    ! print*, num_dt
    ! print 120, t, tspan, dt, y
    ! stop

    ! deallocate(a, b, bstar, c)
    120 format (1x,6(e14.6))
    121 format (a, e14.6)

    return
  end subroutine rk_adaptive

  subroutine rk_explicit(n, t, dt, y, dysub, a, b, c, m)
    ! Dummy arguments
    integer,  intent(in)                      :: n, m
    real(wp), intent(in)                      :: t, dt
    real(wp), intent(inout),  dimension(n)    :: y
    real(wp), intent(in),     dimension(m)    :: b, c
    real(wp), intent(in),     dimension(m,m)  :: a
    procedure(sub_interface)                  :: dysub

    ! Local arguments
    integer                   :: ii
    real(wp), dimension(n)    :: dy
    real(wp), dimension(m,n)  :: k

    ! Initialize guess for the values of k, just k = dy0 = f(t0, y0)
    do ii = 1, m
      call dysub(n, t, y, dy)
      k(ii, :) = dy
    end do

    ! Calculate values of k
    call calc_k(n, t, dt, y, dysub, a, b, c, m, k)

    ! Calculate y(t+dt) using y(t) and k()'s
    y = y + [( dt*dot_product(b,k(:,ii)), ii = 1, n )]

    return
  end subroutine rk_explicit

  subroutine calc_k(n, t, dt, y, dysub, a, b, c, m, k)
    ! Dummy arguments
    integer,  intent(in)                      :: n, m
    real(wp), intent(in)                      :: t, dt
    real(wp), intent(inout),  dimension(n)    :: y
    real(wp), intent(in),     dimension(m)    :: b, c
    real(wp), intent(in),     dimension(m,m)  :: a
    real(wp), intent(inout),  dimension(m,n)  :: k
    procedure(sub_interface)                  :: dysub

    ! Local arguments
    integer                   :: ii, jj
    real(wp)                  :: dummy_t
    real(wp), dimension(n)    :: dummy_y
    real(wp), dimension(n)    :: dy

    dummy_y = y

    do ii = 1, m

      dummy_t = t + dt*c(ii)
      do jj = 1, ii-1
        dummy_y = dummy_y + dt*a(ii,jj)*k(jj,:)
      end do

      call dysub(n, dummy_t, dummy_y, dy)
      k(ii, :) = dy

    end do

    return
  end subroutine calc_k

end module runge_kutta
