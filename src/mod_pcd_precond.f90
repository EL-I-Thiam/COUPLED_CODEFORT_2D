module mod_pcd_precond
  use mod_precision,  only: dp, ip, zero, one, half, two
  use mod_mpi,        only: dx_grid, dy_grid, &
                            mpi_exchange_halos_scalar, mpi_exchange_halos_u, &
                            mpi_exchange_halos_v
  use mod_subsolvers, only: solve_poisson_pressure, solve_helmholtz_u, solve_helmholtz_v
  implicit none
  private

  public :: apply_pcd_precond

contains

  subroutine apply_pcd_precond(z_u, z_v, z_p, r_u, r_v, r_p, &
                               u_adv, v_adv, rho_u, rho_v, &
                               kappa, zeta, eta, dt_step, nx, ny, nh)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)  :: z_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)  :: z_v
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)  :: z_p
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(out) :: r_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(out) :: r_v
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(out) :: r_p
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)  :: u_adv
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)  :: v_adv
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)  :: rho_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)  :: rho_v
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)  :: kappa
    real(dp), dimension(0-nh:nx+nh, 0-nh:ny+nh), intent(in)  :: zeta, eta
    real(dp), intent(in) :: dt_step

    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh) :: phi_p
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh) :: z_u_work
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh) :: z_v_work
    real(dp) :: inv_dx, inv_dy, cross_term, u_avg, v_avg, conv_p
    integer(ip) :: i, j

    inv_dx = one / dx_grid
    inv_dy = one / dy_grid

    r_u = zero
    r_v = zero
    r_p = zero
    phi_p = zero

    ! 1. Résolution de Poisson pour le potentiel de pression : A_p * phi_p = z_p
    call solve_poisson_pressure(phi_p, z_p, rho_u, rho_v, nx, ny, nh)
    call mpi_exchange_halos_scalar(phi_p, nx, ny, nh)

    ! Schur de pression PCD : r_p = - ( (1/dt) * phi_p + N_p(phi_p) )
    do j = 1, ny
      do i = 1, nx
        u_avg = half * (u_adv(i, j) + u_adv(i-1, j))
        v_avg = half * (v_adv(i, j) + v_adv(i, j-1))
        conv_p = u_avg * (phi_p(i+1, j) - phi_p(i-1, j)) * (half * inv_dx) + &
                 v_avg * (phi_p(i, j+1) - phi_p(i, j-1)) * (half * inv_dy)
        r_p(i, j) = - ((one / dt_step) * phi_p(i, j) + conv_p)
      end do
    end do
    call mpi_exchange_halos_scalar(r_p, nx, ny, nh)

    ! 2. Mise à jour du résidu de vitesse : z - B_p^T * r_p
    z_u_work = zero
    do j = 1, ny
      do i = 1, nx-1
        z_u_work(i, j) = z_u(i, j) - (r_p(i+1, j) - r_p(i, j)) * inv_dx
      end do
      z_u_work(0, j)  = z_u(0, j)
      z_u_work(nx, j) = z_u(nx, j)
    end do

    z_v_work = zero
    do j = 1, ny-1
      do i = 1, nx
        z_v_work(i, j) = z_v(i, j) - (r_p(i, j+1) - r_p(i, j)) * inv_dy
      end do
    end do
    do i = 1, nx
      z_v_work(i, 0)  = z_v(i, 0)
      z_v_work(i, ny) = z_v(i, ny)
    end do

    ! 3. Inversion des blocs de vitesse
    call solve_helmholtz_v(r_v, z_v_work, rho_v, kappa, zeta, eta, dt_step, nx, ny, nh, n_sweeps_in=3)
    call mpi_exchange_halos_v(r_v, nx, ny, nh)

    do j = 1, ny
      do i = 1, nx-1
        cross_term = - inv_dy * ( &
          eta(i, j)   * (r_v(i+1, j)   - r_v(i, j))   * inv_dx - &
          eta(i, j-1) * (r_v(i+1, j-1) - r_v(i, j-1)) * inv_dx )
        z_u_work(i, j) = z_u_work(i, j) - cross_term
      end do
    end do

    call solve_helmholtz_u(r_u, z_u_work, rho_u, kappa, zeta, eta, dt_step, nx, ny, nh, n_sweeps_in=3)
    call mpi_exchange_halos_u(r_u, nx, ny, nh)

  end subroutine apply_pcd_precond

end module mod_pcd_precond
