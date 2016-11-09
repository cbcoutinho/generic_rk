module rk_constants
  use iso_fortran_env, only: wp => real64
  implicit none

  public :: classic_rk4, three_eighths_rk4
  public :: midpoint, heun, ralston
  public :: heun_euler, fehlbery_rk12, bogacki_shampine
  public :: fehlbery_rk45, cash_karp_rk45, dormand_prince_rk45
  public :: fehlbery_rk78

  private :: two_stage


contains

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



  subroutine classic_rk4(a, b, c, m)
    integer,  intent(out)                               :: m
    real(wp), intent(out), dimension(:),    allocatable :: c, b
    real(wp), intent(out), dimension(:,:),  allocatable :: a

    m = 4
    allocate(a(m,m), b(m), c(m))

    c = [0._wp, 0.5_wp, 0.5_wp, 1._wp]
    b = [1._wp/6._wp, 1._wp/3._wp, 1._wp/3._wp, 1._wp/6._wp]

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
    b = [1._wp/8._wp, 3._wp/8._wp, 3._wp/8._wp, 1._wp/8._wp]

    a = 0._wp
    a(2,1) = 0.5_wp
    a(3,1:2) = [-1._wp/3._wp, 1._wp]
    a(4,1:3) = [1._wp, -1._wp, 1._wp]

    return
  end subroutine three_eighths_rk4



  subroutine heun_euler(a, b, bstar, c, m, p)
    integer,  intent(out)                               :: m, p
    real(wp), intent(out), dimension(:),    allocatable :: c, b, bstar
    real(wp), intent(out), dimension(:,:),  allocatable :: a

    p = 2
    m = 2
    allocate(a(m,m), b(m), bstar(m), c(m))

    ! There is already a routine to calculate the heun coefficients - just
    ! reuse it for the first row of `b`, then set bstar

    call heun(a, b, c, m)
    bstar = [1._wp, 0._wp]

    return
  end subroutine heun_euler

  subroutine fehlbery_rk12(a, b, bstar, c, m, p)
    integer,  intent(out)                               :: m, p
    real(wp), intent(out), dimension(:),    allocatable :: c, b, bstar
    real(wp), intent(out), dimension(:,:),  allocatable :: a

    p = 2
    m = 3
    allocate(a(m,m), b(m), bstar(m), c(m))

    c       = [0._wp, 0.5_wp, 1.0_wp]
    b       = [1._wp/512._wp, 255._wp/256._wp, 1._wp/512._wp]
    bstar   = [1._wp/256._wp, 255._wp/256._wp, 0._wp]

    a       = 0._wp
    a(2,1)  = 0.5_wp
    a(3,:)  = [1._wp/256._wp, 255._wp/256._wp, 0._wp]

    return
  end subroutine fehlbery_rk12

  subroutine bogacki_shampine(a, b, bstar, c, m, p)
    integer,  intent(out)                               :: m, p
    real(wp), intent(out), dimension(:),    allocatable :: c, b, bstar
    real(wp), intent(out), dimension(:,:),  allocatable :: a

    p = 3
    m = 4
    allocate(a(m,m), b(m), bstar(m), c(m))

    c         = [0._wp, 0.5_wp, 0.75_wp, 1.0_wp]
    b         = [2._wp/9._wp, 1._wp/3._wp, 4._wp/9._wp, 0._wp]
    bstar     = [7._wp/24._wp, 1._wp/4._wp, 1._wp/3._wp, 1._wp/8._wp]

    a         = 0._wp
    a(2,1)    = 0.5_wp
    a(3,2)    = 0.75_wp
    a(4,1:3)  = [2._wp/9._wp, 1._wp/3._wp, 4._wp/9._wp]

    return
  end subroutine bogacki_shampine

  subroutine fehlbery_rk45(a, b, bstar, c, m, p)
    integer,  intent(out)                               :: m, p
    real(wp), intent(out), dimension(:),    allocatable :: c, b, bstar
    real(wp), intent(out), dimension(:,:),  allocatable :: a

    p = 5
    m = 6
    allocate(a(m,m), b(m), bstar(m), c(m))

    c         = [0._wp, 0.25_wp, 3._wp/8._wp, 12._wp/13._wp, 1._wp, 0.5_wp]
    b         = [16._wp/135._wp, 0._wp, 6656._wp/12825._wp, &
                  & 28561._wp/56430._wp, -9._wp/50._wp, 2._wp/55._wp]
    bstar     = [25._wp/216._wp, 0._wp, 1408._wp/2565._wp, 2197._wp/4104._wp, &
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

  subroutine cash_karp_rk45(a, b, bstar, c, m, p)
    integer,  intent(out)                               :: m, p
    real(wp), intent(out), dimension(:),    allocatable :: c, b, bstar
    real(wp), intent(out), dimension(:,:),  allocatable :: a

    p = 5
    m = 6
    allocate(a(m,m), b(m), bstar(m), c(m))

    c         = [0._wp, 1._wp/5._wp, 3._wp/10._wp, 3._wp/5._wp, 1._wp, &
                  & 7._wp/8._wp]
    b         = [37._wp/378._wp, 0._wp, 250._wp/621._wp, 125._wp/594._wp, &
                  & 0._wp, 512._wp/1771._wp]
    bstar     = [2825._wp/27648._wp, 0._wp, 18575._wp/48384._wp, &
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

  subroutine dormand_prince_rk45(a, b, bstar, c, m, p)
    integer,  intent(out)                               :: m, p
    real(wp), intent(out), dimension(:),    allocatable :: c, b, bstar
    real(wp), intent(out), dimension(:,:),  allocatable :: a

    p = 5
    m = 7
    allocate(a(m,m), b(m), bstar(m), c(m))

    c         = [0._wp, 1._wp/5._wp, 3._wp/10._wp, 4._wp/5._wp, 8._wp/9._wp, &
                  & 1._wp, 1._wp]
    b         = [35._wp/384._wp, 0._wp, 500._wp/1113._wp, 125._wp/192._wp, &
                  & -2187._wp/6784._wp, 11._wp/84._wp, 0._wp]
    bstar     = [5179._wp/57600._wp, 0._wp, 7571._wp/16695._wp, &
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

  subroutine fehlbery_rk78(a, b, bstar, c, m, p)
    integer,  intent(out)                               :: m, p
    real(wp), intent(out), dimension(:),    allocatable :: c, b, bstar
    real(wp), intent(out), dimension(:,:),  allocatable :: a

    p = 7
    m = 13
    allocate(a(m,m), b(m), bstar(m), c(m))

    c         = [0._wp, 2._wp/27._wp, 1._wp/9._wp, 1._wp/6._wp, 5._wp/12._wp, &
                & 0.5_wp, 5._wp/6._wp, 1._wp/6._wp, 2._wp/3._wp, &
                & 1._wp/3._wp, 1._wp, 0._wp, 1._wp]
    b         = [41._wp/840._wp, 0._wp, 0._wp, 0._wp, 0._wp, 34._wp/105._wp, &
                & 9._wp/35._wp, 9._wp/35._wp, 9._wp/280._wp, 9._wp/280._wp, &
                & 41._wp/840._wp, 0._wp, 0._wp]
    bstar     = [0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 34._wp/105._wp, &
                & 9._wp/35._wp, 9._wp/35._wp, 9._wp/280._wp, &
                & 9._wp/280._wp, 0._wp, 41._wp/840._wp, 41._wp/840._wp]

    a           = 0._wp
    a(2,1)      = 2._wp/27._wp
    a(3, 1:2)   = [1._wp/36._wp, 1._wp/12._wp]
    a(4, 1:3)   = [1._wp/24._wp, 0._wp, 1._wp/8._wp]
    a(5, 1:4)   = [5._wp/12._wp, 0._wp, -25._wp/16._wp, 25._wp/16._wp]
    a(6, 1:5)   = [1._wp/20._wp, 0._wp, 0._wp, 1._wp/4._wp, 1._wp/5._wp]
    a(7, 1:6)   = [-25._wp/108._wp, 0._wp, 0._wp, 125._wp/108._wp, &
                  & -65._wp/27._wp, 125._wp/54._wp]
    a(8, 1:7)   = [31._wp/300._wp, 0._wp, 0._wp, 0._wp, 61._wp/225._wp, &
                  & -2._wp/9._wp, 13._wp/900._wp]
    a(9, 1:8)   = [2._wp, 0._wp, 0._wp, -53._wp/6._wp, 704._wp/45._wp, &
                  & -107._wp/9._wp, 67._wp/90._wp, 3._wp]
    a(10, 1:9)  = [-91._wp/108._wp, 0._wp, 0._wp, 23._wp/108._wp, &
                  & -976._wp/135._wp, 311._wp/54._wp, -19._wp/60._wp, &
                                & 17._wp/6._wp, -1._wp/12._wp]
    a(11, 1:10) = [2383._wp/4100._wp, 0._wp, 0._wp, -341._wp/164._wp, &
                  & 4496._wp/1025._wp, -301._wp/82._wp, 2133._wp/4100._wp, &
                  & 45._wp/82._wp, 45._wp/164._wp, 18._wp/41._wp]
    a(12, 1:11) = [3._wp/205._wp, 0._wp, 0._wp, 0._wp, 0._wp, -6._wp/41._wp, &
                  & -3._wp/205._wp, -3._wp/41._wp, 3._wp/41._wp, 6._wp/41._wp, &
                  & 0._wp]
    a(13, 1:12) = [-1777._wp/4100._wp, 0._wp, 0._wp, -341._wp/164._wp, &
                  & 4496._wp/1025._wp, -289._wp/82._wp, 2193._wp/4100._wp, &
                  & 51._wp/82._wp, 33._wp/164, 19._wp/41._wp, 0._wp, 1._wp]

    return
  end subroutine fehlbery_rk78

end module rk_constants
