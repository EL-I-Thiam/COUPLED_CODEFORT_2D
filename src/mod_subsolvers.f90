! ==============================================================================
! MODULE : mod_subsolvers
! Solveur 2D Tout-Couple (FUGU-2D) - Fortran 2008 / OpenMPI
!
! Sous-solveurs itératifs scalaires pour les blocs du préconditionneur :
!   1. Inversion du Laplacien variable de pression : A_p^(1/rho) phi = rhs
!      (Gradient Conjugué avec préconditionneur Jacobi et projection moyenne nulle)
!   2. Lissage / Inversion des blocs de vitesse diagonaux : F_uu u = rhs, F_vv v = rhs
!      (Relaxation par Jacobi / Gauss-Seidel amorti)
! ==============================================================================
module mod_subsolvers
  use mod_precision, only: dp, ip, zero, one, half, two
  use mod_params,    only: Nx_global, Ny_global
  use mod_mpi,       only: dx_grid, dy_grid, coords, dims, periods, comm_cart, &
                           mpi_exchange_halos_scalar, mpi_exchange_halos_u, &
                           mpi_exchange_halos_v, mpi_global_dot_product
  use mpi
  implicit none
  private

  public :: solve_poisson_pressure
  public :: solve_helmholtz_u
  public :: solve_helmholtz_v

contains

  ! ----------------------------------------------------------------------------
  ! 1. Solveur de Poisson à masse volumique variable pour la pression
  !    A_p^(1/rho) phi = - div( (1/rho) grad phi ) = rhs
  !    Résolu par Gradient Conjugué (PCG) avec préconditionneur diagonal
  ! ----------------------------------------------------------------------------
  subroutine solve_poisson_pressure(phi, rhs, rho_u, rho_v, nx, ny, nh, &
                                    tol_in, max_iter_in, iter_out, res_out)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: phi
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: rhs
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: rho_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)    :: rho_v
    real(dp), optional, intent(in)  :: tol_in
    integer(ip), optional, intent(in) :: max_iter_in
    integer(ip), optional, intent(out) :: iter_out
    real(dp), optional, intent(out) :: res_out

    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh) :: r, z, p, Ap, diag_inv, rhs_work
    real(dp) :: tol, inv_dx2, inv_dy2
    real(dp) :: alpha, beta, rz_old, rz_new, pAp
    real(dp) :: res_norm, rhs_norm, d_val
    real(dp) :: c_east, c_west, c_north, c_south
    real(dp) :: loc_dot, glob_dot, loc_sum, glob_sum, mean_val
    integer(ip) :: max_iter, iter, i, j, ierr

    res_norm = zero
    tol = 1.0e-5_dp
    if (present(tol_in)) tol = tol_in
    max_iter = 150
    if (present(max_iter_in)) max_iter = max_iter_in

    inv_dx2 = one / (dx_grid * dx_grid)
    inv_dy2 = one / (dy_grid * dy_grid)

    ! a) Calcul de la diagonale inverse pour préconditionnement Jacobi
    diag_inv = one
    do j = 1, ny
      do i = 1, nx
        c_east  = inv_dx2 / rho_u(i, j)
        c_west  = inv_dx2 / rho_u(i-1, j)
        c_north = inv_dy2 / rho_v(i, j)
        c_south = inv_dy2 / rho_v(i, j-1)

        ! Ajustement des parois étanches (flux nul si frontière physique non-périodique)
        if (.not. periods(1)) then
          if (coords(1) == 0 .and. i == 1) c_west = zero
          if (coords(1) == dims(1)-1 .and. i == nx) c_east = zero
        end if
        if (.not. periods(2)) then
          if (coords(2) == 0 .and. j == 1) c_south = zero
          if (coords(2) == dims(2)-1 .and. j == ny) c_north = zero
        end if

        d_val = c_east + c_west + c_north + c_south
        if (d_val > 1.0e-14_dp) then
          diag_inv(i, j) = one / d_val
        else
          diag_inv(i, j) = one
        end if
      end do
    end do

    ! b) Projection compatibilité Neumann (moyenne nulle) et calcul de norme du second membre
    rhs_work = rhs
    if (.not. periods(1) .and. .not. periods(2)) then
      loc_sum = zero
      do j = 1, ny
        do i = 1, nx
          loc_sum = loc_sum + rhs(i, j)
        end do
      end do
      call MPI_Allreduce(loc_sum, glob_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)
      mean_val = glob_sum / real(Nx_global * Ny_global, dp)
      do j = 1, ny
        do i = 1, nx
          rhs_work(i, j) = rhs(i, j) - mean_val
        end do
      end do
    end if

    loc_dot = zero
    do j = 1, ny
      do i = 1, nx
        loc_dot = loc_dot + rhs_work(i, j)**2
      end do
    end do
    call MPI_Allreduce(loc_dot, glob_dot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)
    rhs_norm = sqrt(glob_dot)
    if (rhs_norm < 1.0e-14_dp) then
      phi = zero
      if (present(iter_out)) iter_out = 0
      if (present(res_out)) res_out = zero
      return
    end if

    ! c) Résidu initial r = rhs_work - A * phi
    call mpi_exchange_halos_scalar(phi, nx, ny, nh)
    call apply_neumann_bc_phi(phi, nx, ny, nh)

    do j = 1, ny
      do i = 1, nx
        c_east  = inv_dx2 / rho_u(i, j)
        c_west  = inv_dx2 / rho_u(i-1, j)
        c_north = inv_dy2 / rho_v(i, j)
        c_south = inv_dy2 / rho_v(i, j-1)

        if (.not. periods(1)) then
          if (coords(1) == 0 .and. i == 1) c_west = zero
          if (coords(1) == dims(1)-1 .and. i == nx) c_east = zero
        end if
        if (.not. periods(2)) then
          if (coords(2) == 0 .and. j == 1) c_south = zero
          if (coords(2) == dims(2)-1 .and. j == ny) c_north = zero
        end if

        Ap(i, j) = (c_east + c_west + c_north + c_south) * phi(i, j) - &
                   (c_east * phi(i+1, j) + c_west * phi(i-1, j) + &
                    c_north * phi(i, j+1) + c_south * phi(i, j-1))
        r(i, j) = rhs_work(i, j) - Ap(i, j)
        z(i, j) = diag_inv(i, j) * r(i, j)
        p(i, j) = z(i, j)
      end do
    end do

    loc_dot = zero
    do j = 1, ny
      do i = 1, nx
        loc_dot = loc_dot + r(i, j) * z(i, j)
      end do
    end do
    call MPI_Allreduce(loc_dot, rz_old, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)

    ! d) Boucle d'itérations Gradient Conjugué (PCG)
    do iter = 1, max_iter
      call mpi_exchange_halos_scalar(p, nx, ny, nh)
      call apply_neumann_bc_phi(p, nx, ny, nh)

      ! Mat-vec Ap = A * p
      loc_dot = zero
      do j = 1, ny
        do i = 1, nx
          c_east  = inv_dx2 / rho_u(i, j)
          c_west  = inv_dx2 / rho_u(i-1, j)
          c_north = inv_dy2 / rho_v(i, j)
          c_south = inv_dy2 / rho_v(i, j-1)

          if (.not. periods(1)) then
            if (coords(1) == 0 .and. i == 1) c_west = zero
            if (coords(1) == dims(1)-1 .and. i == nx) c_east = zero
          end if
          if (.not. periods(2)) then
            if (coords(2) == 0 .and. j == 1) c_south = zero
            if (coords(2) == dims(2)-1 .and. j == ny) c_north = zero
          end if

          Ap(i, j) = (c_east + c_west + c_north + c_south) * p(i, j) - &
                     (c_east * p(i+1, j) + c_west * p(i-1, j) + &
                      c_north * p(i, j+1) + c_south * p(i, j-1))
          loc_dot = loc_dot + p(i, j) * Ap(i, j)
        end do
      end do

      call MPI_Allreduce(loc_dot, pAp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)
      if (abs(pAp) < 1.0e-30_dp) exit

      alpha = rz_old / pAp

      loc_dot = zero
      do j = 1, ny
        do i = 1, nx
          phi(i, j) = phi(i, j) + alpha * p(i, j)
          r(i, j)   = r(i, j) - alpha * Ap(i, j)
          z(i, j)   = diag_inv(i, j) * r(i, j)
          loc_dot   = loc_dot + r(i, j)**2
        end do
      end do

      call MPI_Allreduce(loc_dot, glob_dot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)
      res_norm = sqrt(glob_dot)

      if (res_norm / rhs_norm < tol) exit

      loc_dot = zero
      do j = 1, ny
        do i = 1, nx
          loc_dot = loc_dot + r(i, j) * z(i, j)
        end do
      end do
      call MPI_Allreduce(loc_dot, rz_new, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)

      beta = rz_new / rz_old
      rz_old = rz_new

      do j = 1, ny
        do i = 1, nx
          p(i, j) = z(i, j) + beta * p(i, j)
        end do
      end do
    end do

    if (.not. periods(1) .and. .not. periods(2)) then
      loc_sum = zero
      do j = 1, ny
        do i = 1, nx
          loc_sum = loc_sum + phi(i, j)
        end do
      end do
      call MPI_Allreduce(loc_sum, glob_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)
      mean_val = glob_sum / real(Nx_global * Ny_global, dp)
      do j = 1, ny
        do i = 1, nx
          phi(i, j) = phi(i, j) - mean_val
        end do
      end do
    end if

    if (present(iter_out)) iter_out = min(iter, max_iter)
    if (present(res_out))  res_out  = res_norm / rhs_norm

  end subroutine solve_poisson_pressure

  ! ----------------------------------------------------------------------------
  ! 2. Solveur / Lisseur scalaire pour le bloc de vitesse U : F_uu u = rhs
  !    Par Jacobi amorti (relaxation sous-relaxée pour assurer la contraction)
  ! ----------------------------------------------------------------------------
  subroutine solve_helmholtz_u(u_sol, rhs_u, rho_u, kappa, zeta, eta, dt_step, &
                               nx, ny, nh, n_sweeps_in)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: u_sol
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: rhs_u
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: rho_u
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: kappa
    real(dp), dimension(0-nh:nx+nh, 0-nh:ny+nh), intent(in)    :: zeta, eta
    real(dp), intent(in) :: dt_step
    integer(ip), optional, intent(in) :: n_sweeps_in

    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh) :: u_old
    real(dp) :: inv_dx2, inv_dy2, diag_val, off_diag
    real(dp) :: c_E, c_W, c_N, c_S
    real(dp) :: omega_relax
    integer(ip) :: n_sweeps, sweep, i, j

    n_sweeps = 5
    if (present(n_sweeps_in)) n_sweeps = n_sweeps_in
    omega_relax = 0.8_dp

    inv_dx2 = one / (dx_grid * dx_grid)
    inv_dy2 = one / (dy_grid * dy_grid)

    do sweep = 1, n_sweeps
      call mpi_exchange_halos_u(u_sol, nx, ny, nh)
      u_old = u_sol

      do j = 1, ny
        do i = 1, nx-1
          ! Coefficients diagonaux et hors-diagonaux de la diffusion
          c_E = kappa(i+1, j) * inv_dx2
          c_W = kappa(i, j)   * inv_dx2
          c_N = (zeta(i, j)   - eta(i, j))   * inv_dy2
          c_S = (zeta(i, j-1) - eta(i, j-1)) * inv_dy2

          diag_val = (rho_u(i, j) / dt_step) + c_E + c_W + c_N + c_S
          off_diag = c_E * u_old(i+1, j) + c_W * u_old(i-1, j) + &
                     c_N * u_old(i, j+1) + c_S * u_old(i, j-1)

          if (diag_val > 1.0e-15_dp) then
            u_sol(i, j) = (one - omega_relax) * u_old(i, j) + &
                          omega_relax * ((rhs_u(i, j) + off_diag) / diag_val)
          end if
        end do
        if (.not. periods(1)) then
          if (coords(1) == 0) u_sol(0, j) = rhs_u(0, j)
          if (coords(1) == dims(1)-1) u_sol(nx, j) = rhs_u(nx, j)
        end if
      end do
    end do

  end subroutine solve_helmholtz_u

  ! ----------------------------------------------------------------------------
  ! 3. Solveur / Lisseur scalaire pour le bloc de vitesse V : F_vv v = rhs
  ! ----------------------------------------------------------------------------
  subroutine solve_helmholtz_v(v_sol, rhs_v, rho_v, kappa, zeta, eta, dt_step, &
                               nx, ny, nh, n_sweeps_in)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(inout) :: v_sol
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)    :: rhs_v
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)    :: rho_v
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: kappa
    real(dp), dimension(0-nh:nx+nh, 0-nh:ny+nh), intent(in)    :: zeta, eta
    real(dp), intent(in) :: dt_step
    integer(ip), optional, intent(in) :: n_sweeps_in

    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh) :: v_old
    real(dp) :: inv_dx2, inv_dy2, diag_val, off_diag
    real(dp) :: c_E, c_W, c_N, c_S
    real(dp) :: omega_relax
    integer(ip) :: n_sweeps, sweep, i, j

    n_sweeps = 5
    if (present(n_sweeps_in)) n_sweeps = n_sweeps_in
    omega_relax = 0.8_dp

    inv_dx2 = one / (dx_grid * dx_grid)
    inv_dy2 = one / (dy_grid * dy_grid)

    do sweep = 1, n_sweeps
      call mpi_exchange_halos_v(v_sol, nx, ny, nh)
      v_old = v_sol

      do j = 1, ny-1
        do i = 1, nx
          c_E = (zeta(i, j)   - eta(i, j))   * inv_dx2
          c_W = (zeta(i-1, j) - eta(i-1, j)) * inv_dx2
          c_N = kappa(i, j+1) * inv_dy2
          c_S = kappa(i, j)   * inv_dy2

          diag_val = (rho_v(i, j) / dt_step) + c_E + c_W + c_N + c_S
          off_diag = c_E * v_old(i+1, j) + c_W * v_old(i-1, j) + &
                     c_N * v_old(i, j+1) + c_S * v_old(i, j-1)

          if (diag_val > 1.0e-15_dp) then
            v_sol(i, j) = (one - omega_relax) * v_old(i, j) + &
                          omega_relax * ((rhs_v(i, j) + off_diag) / diag_val)
          end if
        end do
      end do
      if (.not. periods(2)) then
        do i = 1, nx
          if (coords(2) == 0) v_sol(i, 0) = rhs_v(i, 0)
          if (coords(2) == dims(2)-1) v_sol(i, ny) = rhs_v(i, ny)
        end do
      end if
    end do

  end subroutine solve_helmholtz_v

  ! ----------------------------------------------------------------------------
  ! Routine utilitaire interne : conditions aux limites Neumann pour la pression
  ! (d phi / dn = 0 sur les frontières externes non périodiques)
  ! ----------------------------------------------------------------------------
  subroutine apply_neumann_bc_phi(phi, nx, ny, nh)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: phi
    integer(ip) :: i, j

    if (.not. periods(1)) then
      if (coords(1) == 0) then
        do j = 1-nh, ny+nh
          do i = 1, nh
            phi(1-i, j) = phi(i, j)
          end do
        end do
      end if
      if (coords(1) == dims(1)-1) then
        do j = 1-nh, ny+nh
          do i = 1, nh
            phi(nx+i, j) = phi(nx-i+1, j)
          end do
        end do
      end if
    end if

    if (.not. periods(2)) then
      if (coords(2) == 0) then
        do i = 1-nh, nx+nh
          do j = 1, nh
            phi(i, 1-j) = phi(i, j)
          end do
        end do
      end if
      if (coords(2) == dims(2)-1) then
        do i = 1-nh, nx+nh
          do j = 1, nh
            phi(i, ny+j) = phi(i, ny-j+1)
          end do
        end do
      end if
    end if
  end subroutine apply_neumann_bc_phi

end module mod_subsolvers
