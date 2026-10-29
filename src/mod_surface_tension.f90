! ==============================================================================
! MODULE : mod_surface_tension
! Solveur 2D Tout-Couplé (FUGU-2D) - Fortran 2008 / OpenMPI
!
! Calcul de la normale unitaire, de la courbure kappa_Gamma et de la force
! volumique de tension superficielle CSF de Brackbill et al. (1992)
! sur maillage décalé (MAC).
! ==============================================================================
module mod_surface_tension
  use mod_precision, only: dp, ip, zero, one, half, two
  use mod_params,    only: sigma
  use mod_mpi,       only: dx_grid, dy_grid, coords, dims, periods, comm_cart, &
                           mpi_exchange_halos_scalar
  implicit none
  private

  public :: compute_interface_curvature
  public :: compute_csf_force

contains

  ! ----------------------------------------------------------------------------
  ! 1. Calcul de la normale unitaire et de la courbure kappa_Gamma = - div(n)
  ! ----------------------------------------------------------------------------
  subroutine compute_interface_curvature(C_smooth, kappa_curv, nx, ny, nh)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)  :: C_smooth
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(out) :: kappa_curv

    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh) :: nx_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh) :: ny_v
    real(dp) :: grad_x, grad_y, norm_grad, eps_reg
    real(dp) :: inv_dx, inv_dy
    integer(ip) :: i, j

    inv_dx = one / dx_grid
    inv_dy = one / dy_grid
    eps_reg = 1.0e-12_dp

    ! a) Normale n_x sur les faces u : (i+1/2, j)
    do j = 1, ny
      do i = 0, nx
        grad_x = (C_smooth(i+1, j) - C_smooth(i, j)) * inv_dx
        grad_y = 0.25_dp * (C_smooth(i, j+1) + C_smooth(i+1, j+1) - &
                            C_smooth(i, j-1) - C_smooth(i+1, j-1)) * inv_dy
        norm_grad = sqrt(grad_x**2 + grad_y**2 + eps_reg)
        nx_u(i, j) = grad_x / norm_grad
      end do
    end do

    ! b) Normale n_y sur les faces v : (i, j+1/2)
    do j = 0, ny
      do i = 1, nx
        grad_x = 0.25_dp * (C_smooth(i+1, j) + C_smooth(i+1, j+1) - &
                            C_smooth(i-1, j) - C_smooth(i-1, j+1)) * inv_dx
        grad_y = (C_smooth(i, j+1) - C_smooth(i, j)) * inv_dy
        norm_grad = sqrt(grad_x**2 + grad_y**2 + eps_reg)
        ny_v(i, j) = grad_y / norm_grad
      end do
    end do

    ! c) Courbure kappa = - div(n) aux centres des mailles
    kappa_curv = zero
    do j = 1, ny
      do i = 1, nx
        kappa_curv(i, j) = - ( (nx_u(i, j) - nx_u(i-1, j)) * inv_dx + &
                               (ny_v(i, j) - ny_v(i, j-1)) * inv_dy )
      end do
    end do

    call mpi_exchange_halos_scalar(kappa_curv, nx, ny, nh)
  end subroutine compute_interface_curvature

  ! ----------------------------------------------------------------------------
  ! 2. Force volumique CSF : f_sigma = sigma * kappa_Gamma * grad(C)
  !    Discrétisée directement sur les faces de vitesse MAC décalées
  ! ----------------------------------------------------------------------------
  subroutine compute_csf_force(C_raw, kappa_curv, f_sigma_u, f_sigma_v, nx, ny, nh)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)  :: C_raw
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)  :: kappa_curv
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(out) :: f_sigma_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(out) :: f_sigma_v

    real(dp) :: inv_dx, inv_dy, kappa_face
    integer(ip) :: i, j

    inv_dx = one / dx_grid
    inv_dy = one / dy_grid

    f_sigma_u = zero
    f_sigma_v = zero

    ! Force sur les faces u
    do j = 1, ny
      do i = 0, nx
        kappa_face = half * (kappa_curv(i, j) + kappa_curv(i+1, j))
        f_sigma_u(i, j) = sigma * kappa_face * (C_raw(i+1, j) - C_raw(i, j)) * inv_dx
      end do
    end do

    ! Force sur les faces v
    do j = 0, ny
      do i = 1, nx
        kappa_face = half * (kappa_curv(i, j) + kappa_curv(i, j+1))
        f_sigma_v(i, j) = sigma * kappa_face * (C_raw(i, j+1) - C_raw(i, j)) * inv_dy
      end do
    end do
  end subroutine compute_csf_force

end module mod_surface_tension
