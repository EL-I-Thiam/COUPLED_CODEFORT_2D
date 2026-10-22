! ==============================================================================
! MODULE : mod_viscosity
! Solveur 2D Tout-Couple (FUGU-2D) - Fortran 2008 / OpenMPI
!
! Calcul des proprietes physiques effectives et decomposition du tenseur visqueux:
!   - Masse volumique rho(C) : centres et faces de vitesse (moyenne arithmetique)
!   - Viscosite d'elongation kappa = 2*mu(C) au centre des mailles
!   - Viscosites de cisaillement zeta = 2*mu et de rotation eta = mu aux sommets
!     calculees par moyenne harmonique d'apres la these de M. El Ouafa (Ch. 4.4.2)
! ==============================================================================
module mod_viscosity
  use mod_precision, only: dp, ip, zero, one, half, two
  use mod_params,    only: rho1, rho2, mu1, mu2
  implicit none
  private

  public :: compute_mixture_densities
  public :: compute_mixture_viscosities
  public :: harmonic_viscosity

contains

  ! ----------------------------------------------------------------------------
  ! Fonction pure : moyenne harmonique de deux viscosites
  ! mu(C) = (mu1 * mu2) / (mu2 * C + mu1 * (1 - C))
  ! Assure l'ordre 2 au passage de l'interface liquide-gaz ou liquide-liquide
  ! ----------------------------------------------------------------------------
  pure function harmonic_viscosity(c_val) result(mu_eff)
    real(dp), intent(in) :: c_val
    real(dp) :: mu_eff
    real(dp) :: c_clamped, denom

    ! Bornage de securite pour eviter toute extrapolation non-physique hors [0, 1]
    c_clamped = max(zero, min(one, c_val))

    if (abs(mu1 - mu2) < 1.0e-14_dp) then
      mu_eff = mu1
    else
      denom = mu2 * c_clamped + mu1 * (one - c_clamped)
      if (denom > 1.0e-15_dp) then
        mu_eff = (mu1 * mu2) / denom
      else
        mu_eff = min(mu1, mu2)
      end if
    end if
  end function harmonic_viscosity

  ! ----------------------------------------------------------------------------
  ! Calcul des masses volumiques effectives
  !   - rho_c au centre des mailles : rho_c = rho1 * C + rho2 * (1 - C)
  !   - rho_u sur les faces verticales (decalage en x)
  !   - rho_v sur les faces horizontales (decalage en y)
  ! ----------------------------------------------------------------------------
  subroutine compute_mixture_densities(c_field, rho_c, rho_u, rho_v, nx, ny, nh)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: c_field
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: rho_c
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: rho_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(inout) :: rho_v

    integer(ip) :: i, j
    real(dp) :: c_val

    ! 1. Masse volumique centree
    do j = 1-nh, ny+nh
      do i = 1-nh, nx+nh
        c_val = max(zero, min(one, c_field(i, j)))
        rho_c(i, j) = rho1 * c_val + rho2 * (one - c_val)
      end do
    end do

    ! 2. Masse volumique sur les faces u (i+1/2, j) : moyenne arithmetique
    !    Les faces i utilisent les cellules i et i+1, valides pour i dans [1-nh, nx+nh-1]
    do j = 1-nh, ny+nh
      do i = 1-nh, nx+nh-1
        rho_u(i, j) = half * (rho_c(i, j) + rho_c(i+1, j))
      end do
      ! Extension sur les bords extrêmes
      rho_u(0-nh, j)  = rho_u(1-nh, j)
      rho_u(nx+nh, j) = rho_u(nx+nh-1, j)
    end do

    ! 3. Masse volumique sur les faces v (i, j+1/2) : moyenne arithmetique
    !    Les faces j utilisent les cellules j et j+1, valides pour j dans [1-nh, ny+nh-1]
    do i = 1-nh, nx+nh
      do j = 1-nh, ny+nh-1
        rho_v(i, j) = half * (rho_c(i, j) + rho_c(i, j+1))
      end do
      rho_v(i, 0-nh)  = rho_v(i, 1-nh)
      rho_v(i, ny+nh) = rho_v(i, ny+nh-1)
    end do

  end subroutine compute_mixture_densities

  ! ----------------------------------------------------------------------------
  ! Calcul des viscosites decomposees (Caltagirone & Vincent / El Ouafa Ch. 4)
  !   - kappa = 2 * mu_c(i, j) aux centres des mailles (elongation)
  !   - zeta  = 2 * mu_vert aux sommets des mailles (cisaillement)
  !   - eta   = mu_vert aux sommets des mailles (rotation)
  !
  ! Indication d'indexation des sommets :
  !   Le sommet (i, j) correspond geometriquement au coin (i+1/2, j+1/2)
  !   entoure par les 4 centres (i, j), (i+1, j), (i, j+1), (i+1, j+1).
  ! ----------------------------------------------------------------------------
  subroutine compute_mixture_viscosities(c_field, kappa, zeta, eta, nx, ny, nh)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: c_field
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: kappa
    real(dp), dimension(0-nh:nx+nh, 0-nh:ny+nh), intent(inout) :: zeta
    real(dp), dimension(0-nh:nx+nh, 0-nh:ny+nh), intent(inout) :: eta

    integer(ip) :: i, j
    real(dp) :: c_vert, mu_node, mu_vert

    ! 1. Viscosite d'elongation kappa = 2*mu au centre des cellules
    do j = 1-nh, ny+nh
      do i = 1-nh, nx+nh
        mu_node = harmonic_viscosity(c_field(i, j))
        kappa(i, j) = two * mu_node
      end do
    end do

    ! 2. Viscosites de cisaillement zeta et de rotation eta aux sommets
    !    Les 4 centres (i, j), (i+1, j), (i, j+1), (i+1, j+1) existent
    !    pour i dans [1-nh, nx+nh-1] et j dans [1-nh, ny+nh-1]
    do j = 1-nh, ny+nh-1
      do i = 1-nh, nx+nh-1
        c_vert = 0.25_dp * (c_field(i, j)   + c_field(i+1, j) + &
                            c_field(i, j+1) + c_field(i+1, j+1))
        mu_vert = harmonic_viscosity(c_vert)
        zeta(i, j) = two * mu_vert
        eta(i, j)  = mu_vert
      end do
    end do

    ! Remplissage des bordures extremes
    do j = 1-nh, ny+nh-1
      zeta(0-nh, j)  = zeta(1-nh, j)
      zeta(nx+nh, j) = zeta(nx+nh-1, j)
      eta(0-nh, j)   = eta(1-nh, j)
      eta(nx+nh, j)  = eta(nx+nh-1, j)
    end do

    do i = 0-nh, nx+nh
      zeta(i, 0-nh)  = zeta(i, 1-nh)
      zeta(i, ny+nh) = zeta(i, ny+nh-1)
      eta(i, 0-nh)   = eta(i, 1-nh)
      eta(i, ny+nh)  = eta(i, ny+nh-1)
    end do

  end subroutine compute_mixture_viscosities

end module mod_viscosity
