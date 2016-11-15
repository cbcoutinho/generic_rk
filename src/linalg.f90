module linalg
  use iso_fortran_env, only: wp => real64
  implicit none

  public  :: linsolve_quick, linsolve

contains
  subroutine linsolve_quick(n, a, nrhs, b, x)

    ! Quick wrapper around linsolve

    integer,  intent(in)                          :: n, nrhs
    real(wp), intent(in),     dimension(n, n)     :: a
    real(wp), intent(in),     dimension(n, nrhs)  :: b
    real(wp), intent(out),    dimension(n, nrhs)  :: x

    integer,                  dimension(n)        :: P
    real(wp),                 dimension(n, n)     :: LU

    call linsolve (n, a, nrhs, b, x, LU, P, .False.)
    return
  end subroutine linsolve_quick

  subroutine linsolve (n, a, nrhs, b, x, LU, P, toggle)
    ! This routine is a wrapper dgesv, splitting it into its two primary
    ! components:
    !             dgetrf - Decomposes A into P*L*U
    !             dgetrs - Uses P*L*U to solve for x (Ax=b => (P*L*U)x=b)
    !
    ! Splitting these two up like this allows you to save the decomposed
    ! version of 'a' so that you don't have to do it again. If 'toggle' is
    ! equal to true, then the decomposition has already occured and LU can be
    ! trusted - OK to skip dgetrf

    ! Dummy variables
    integer,  intent(in)                          :: n, nrhs
    integer,  intent(inout),  dimension(n)        :: P
    real(wp), intent(in),     dimension(n, nrhs)  :: b
    real(wp), intent(in),     dimension(n, n)     :: a
    real(wp), intent(out),    dimension(n, nrhs)  :: x
    real(wp), intent(inout),  dimension(n, n)     :: LU
    logical,  intent(in)                          :: toggle

    ! Local variables
    integer                                       :: info
    integer,                  dimension(n)        :: my_P
    real(wp),                 dimension(n, n)     :: my_a
    real(wp),                 dimension(n, nrhs)  :: my_b

    my_a = a
    my_b = b

    if ( .not. toggle ) then
      call dgetrf (n, n, my_a, n, my_P, info)

      if ( info /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DGETRF'
        write ( *, '(a,i8)' ) '  Factorization failed, INFO = ', info
        stop
      end if

      LU  = my_a
      P   = my_P

    end if

    call dgetrs ('No transpose', n, nrhs, LU, n, P, my_b, n, info)

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DGETRS'
      write ( *, '(a,i8)' ) '  Back substitution failed, INFO = ', info
      stop
    end if

    x = my_b

    return
  end subroutine linsolve
end module linalg
