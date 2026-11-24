! ==============================================================================
! PROGRAMME DE TEST : test_poiseuille_twophase.f90
! Validation analytique de l'écoulement de Poiseuille diphasique (Semaine 6)
! Vérifie la continuité de la vitesse et le saut de pente à l'interface
! grâce à la moyenne harmonique de la viscosité sur maillage MAC.
! Projet Solveur 2D Tout-Couplé (FUGU-2D) - Fortran 2008 / OpenMPI
! ==============================================================================
program test_poiseuille_twophase
  use mod_precision,      only: dp, ip, zero, one, half, two
  use mod_params,         only: Nx_global, Ny_global, Lx, Ly, rho1, rho2, &
                                mu1, mu2, sigma, gx, gy, time_order, &
                                tol_bicgstab2, max_iter_bicgstab2
  use mod_mpi,            only: my_rank, n_procs, dims, nx_loc, ny_loc, &
                                j_start, dx_grid, dy_grid, comm_cart, &
                                mpi_init_env, mpi_finalize_env, &
                                mpi_setup_cartesian_mesh, mpi_exchange_halos_scalar
  use mod_viscosity,      only: compute_mixture_densities, compute_mixture_viscosities
  use mod_solver_coupled, only: advance_coupled_step
  use mpi
  implicit none

  integer(ip), parameter :: nh = 2
  real(dp), allocatable :: u_vel(:,:), v_vel(:,:), u_prev(:,:), v_prev(:,:)
  real(dp), allocatable :: p_field(:,:), c_field(:,:)
  real(dp), allocatable :: rho_c(:,:), rho_u(:,:), rho_v(:,:)
  real(dp), allocatable :: kappa(:,:), zeta(:,:), eta(:,:)
  real(dp), allocatable :: kappa_curv(:,:), f_sigma_u(:,:), f_sigma_v(:,:)

  real(dp) :: y_c, y_int, G_press, A1, A2, B2, u_exact, dt_step, res_sol
  real(dp) :: err_l2_loc, err_l2_glob, norm_l2_loc, norm_l2_glob
  real(dp) :: err_inf_loc, err_inf_glob, rel_l2
  integer(ip) :: i_loc, j_loc, j_glob, step, n_steps, iter_sol, ierr

  call mpi_init_env()

  Nx_global = 16
  Ny_global = 32
  Lx        = one
  Ly        = one
  rho1      = 1.0_dp
  rho2      = 1.0_dp
  mu1       = 1.0_dp
  mu2       = 0.25_dp
  sigma     = zero
  G_press   = 1.0_dp
  gx        = G_press
  gy        = zero
  time_order = 2
  tol_bicgstab2 = 1.0e-7_dp
  max_iter_bicgstab2 = 200

  ! Périodique en X, parois adhérentes en Y (y=0 et y=Ly)
  call mpi_setup_cartesian_mesh(Nx_global, Ny_global, Lx, Ly, periodic_x=.true., periodic_y=.false.)

  if (my_rank == 0) then
    write(*,'(A)') "================================================================"
    write(*,'(A)') " TEST SEMAINE 6 : Poiseuille Diphasique Analytique (mu1/mu2=4)  "
    write(*,'(A)') "================================================================"
    write(*,'(A,I4,A,I3,A,I3)') " Processus MPI : ", n_procs, " | Grille : ", dims(1), " x ", dims(2)
    write(*,'(A,I4,A,I4,A,F5.2,A,F5.2)') " Maillage : ", Nx_global, " x ", Ny_global, &
                                         " | mu1 = ", mu1, " | mu2 = ", mu2
  end if

  allocate(u_vel(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), v_vel(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(u_prev(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), v_prev(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(p_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh), c_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(rho_c(1-nh:nx_loc+nh, 1-nh:ny_loc+nh), rho_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), rho_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(kappa(1-nh:nx_loc+nh, 1-nh:ny_loc+nh), zeta(0-nh:nx_loc+nh, 0-nh:ny_loc+nh), eta(0-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(kappa_curv(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(f_sigma_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), f_sigma_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))

  u_vel   = zero
  v_vel   = zero
  u_prev  = zero
  v_prev  = zero
  p_field = zero

  y_int = half * Ly
  do j_loc = 1-nh, ny_loc+nh
    j_glob = j_start + j_loc - 1
    y_c = (real(j_glob, dp) - half) * dy_grid
    do i_loc = 1-nh, nx_loc+nh
      if (y_c <= y_int) then
        c_field(i_loc, j_loc) = one
      else
        c_field(i_loc, j_loc) = zero
      end if
    end do
  end do
  call mpi_exchange_halos_scalar(c_field, nx_loc, ny_loc, nh)

  call compute_mixture_densities(c_field, rho_c, rho_u, rho_v, nx_loc, ny_loc, nh)
  call compute_mixture_viscosities(c_field, kappa, zeta, eta, nx_loc, ny_loc, nh)

  ! Intégration jusqu'au régime permanent avec BDF2
  dt_step = 0.05_dp
  n_steps = 80
  do step = 1, n_steps
    call advance_coupled_step(u_vel, v_vel, p_field, c_field, u_prev, v_prev, &
                              rho_c, rho_u, rho_v, kappa, zeta, eta, &
                              kappa_curv, f_sigma_u, f_sigma_v, &
                              dt_step, step, nx_loc, ny_loc, nh, &
                              iter_out=iter_sol, res_out=res_sol, use_bicgstab2=.false.)
  end do

  ! Coefficients de la solution analytique exacte par morceaux
  A1 = G_press * (y_int**2 / (two * mu1) + (Ly**2 - y_int**2) / (two * mu2)) / &
       (y_int + (mu1 / mu2) * (Ly - y_int))
  A2 = (mu1 / mu2) * A1
  B2 = G_press * Ly**2 / (two * mu2) - A2 * Ly

  err_l2_loc  = zero
  norm_l2_loc = zero
  err_inf_loc = zero

  do j_loc = 1, ny_loc
    j_glob = j_start + j_loc - 1
    y_c = (real(j_glob, dp) - half) * dy_grid
    if (y_c <= y_int) then
      u_exact = - (G_press / (two * mu1)) * y_c**2 + A1 * y_c
    else
      u_exact = - (G_press / (two * mu2)) * y_c**2 + A2 * y_c + B2
    end if
    do i_loc = 1, nx_loc
      err_l2_loc  = err_l2_loc  + (u_vel(i_loc, j_loc) - u_exact)**2 * dx_grid * dy_grid
      norm_l2_loc = norm_l2_loc + u_exact**2 * dx_grid * dy_grid
      err_inf_loc = max(err_inf_loc, abs(u_vel(i_loc, j_loc) - u_exact))
    end do
  end do

  call MPI_Allreduce(err_l2_loc,  err_l2_glob,  1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)
  call MPI_Allreduce(norm_l2_loc, norm_l2_glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)
  call MPI_Allreduce(err_inf_loc, err_inf_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_cart, ierr)

  rel_l2 = sqrt(err_l2_glob / max(norm_l2_glob, 1.0e-14_dp))

  if (my_rank == 0) then
    write(*,'(A)') "----------------------------------------------------------------"
    write(*,'(A,ES12.4)') " Erreur relative L2 vs solution analytique : ", rel_l2
    write(*,'(A,ES12.4)') " Erreur max L_inf                          : ", err_inf_glob
    write(*,'(A,I4,A,ES10.3)') " Dernier pas Krylov : iter = ", iter_sol, " | res = ", res_sol
    write(*,'(A)') "----------------------------------------------------------------"
    if (rel_l2 < 1.0e-2_dp .and. err_inf_glob < 5.0e-3_dp) then
      write(*,'(A)') " --> [SUCCÈS] Poiseuille diphasique conforme à la solution exacte !"
    else
      write(*,'(A)') " --> [ATTENTION] Écart par rapport au profil analytique."
    end if
    write(*,'(A)') "================================================================"
  end if

  deallocate(u_vel, v_vel, u_prev, v_prev, p_field, c_field)
  deallocate(rho_c, rho_u, rho_v, kappa, zeta, eta)
  deallocate(kappa_curv, f_sigma_u, f_sigma_v)

  call mpi_finalize_env()
end program test_poiseuille_twophase
