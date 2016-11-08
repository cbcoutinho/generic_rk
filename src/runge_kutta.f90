module runge_kutta
  use iso_fortran_env, only: wp => real64
  use rk_constants, only: classic_rk4, three_eighths_rk4
  use rk_constants, only: midpoint, heun, ralston
  use rk_constants, only: heun_euler, fehlbery_rk12, bogacki_shampine
  use rk_constants, only: fehlbery_rk45, cash_karp_rk45, dormand_prince_rk45
  implicit none

contains

  subroutine rk_wrapper(n, num_t, t, y, dysub)
    ! Dummy arguments
    integer,  intent(in)                          :: n, num_t
    real(wp), intent(inout),  dimension(num_t)    :: t
    real(wp),                 dimension(num_t, n) :: y

    interface
      subroutine dysub(n, t, y, dy)
        import wp
        integer,  intent(in)                      :: n
        real(wp), intent(in)                      :: t
        real(wp), intent(in),   dimension(n)      :: y
        real(wp), intent(out),  dimension(n)      :: dy
      end subroutine
    end interface

    ! Local arguments
    integer                                       :: ii, m, p
    real(wp)                                      :: dt
    real(wp), dimension(2)                        :: tspan
    real(wp), dimension(n)                        :: yy
    real(wp), dimension(:),   allocatable         :: c, b, bstar
    real(wp), dimension(:,:), allocatable         :: a

    ! Get the tableau coefficients associated with a runge kutta implementation
    call dormand_prince_rk45(a, b, bstar, c, m, p)
    ! call heun(a, b, c, m, p)

    ! Initial guess for dt is just the difference between t(1:2), this way new
    ! versions of dt will be saved in each call to rk_adaptive
    dt = t(2) - t(1)

    do ii = 1, num_t-1
      yy = y(ii,:)
      tspan = t(ii:ii+1)

      call rk_adaptive(n, tspan, dt, yy, dysub, a, b, bstar, c, m, p)
      ! call rk_explicit(n, t(ii), t(ii+1)-t(ii), yy, dysub, a, b, c, m, p)

      y(ii+1, :) = yy

      ! if (ii == 1) stop
    end do

    deallocate(a, b, bstar, c)

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
    integer                               :: ii, jj, kk
    real(wp)                              :: t, error
    real(wp), parameter                   :: eps1 = sqrt(epsilon(1._wp))
    real(wp), dimension(n)                :: dummy_y, y_star

    integer                               :: eval, num_dt
    integer, parameter                    :: max_eval = 100
    real(wp)                              :: s, new_s
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

        s = eps1*dt / &
            & (2._wp * (tspan(2) - tspan(1)) * error)
        ! print*, 'error    = ', error
        ! print*, 's        = ', s
        ! print*, 'sqrt(s)  = ', sqrt(s)
        ! print*,

        new_s = 0.9_wp * (eps1/error) ** (1._wp/real(p, wp))
        ! print*, 'new_s    = ', new_s
        ! stop

        ! Adjust dt based on estimation of error
        if ( new_s >= 2._wp ) then
          t = t + dt
          y = y_star
          dt = dt * 2._wp
          ! print*, 'Adjusted dt with factor of 2'
          exit
        else if ( new_s < 2._wp .and. new_s >= 1._wp ) then
          t = t + dt
          y = y_star
          ! dt = dt * s
          ! print 121, 'Adjusted dt with factor of ', s
          ! dt = dt * new_s
          ! print 121, 'Adjusted dt with factor of ', new_s
          exit
        else if ( new_s <= 1._wp ) then
          ! dt = dt * s
          ! print 121, 'Adjusted dt with factor of ', s
          dt = dt * new_s
          ! print 121, 'Adjusted dt with factor of ', new_s
        end if

        if ( eval == max_eval ) then
          print*, "Warning, max eval iterations reached!"
        end if

      end do

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

end module runge_kutta
