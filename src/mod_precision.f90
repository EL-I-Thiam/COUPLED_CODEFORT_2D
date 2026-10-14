module mod_precision
  use, intrinsic :: iso_fortran_env, only: real64, int32, int64
  implicit none
  public

  ! Précision double (64 bits, ~15-17 chiffres significatifs)
  integer, parameter :: dp = real64
  integer, parameter :: ip = int32
  integer, parameter :: lp = int64

  ! Constantes numériques fondamentales
  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: one  = 1.0_dp
  real(dp), parameter :: two  = 2.0_dp
  real(dp), parameter :: half = 0.5_dp
  real(dp), parameter :: pi   = 3.141592653589793238462643383279502884197_dp

end module mod_precision
