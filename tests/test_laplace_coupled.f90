! ==============================================================================
! PROGRAMME DE TEST : test_laplace_coupled.f90
! Validation du couplage complet Navier-Stokes + VOF + CSF (Semaine 5)
! Cas-test : Goutte statique 2D — Loi de Laplace & Courants parasites
! Projet Solveur 2D Tout-Couplé (FUGU-2D) - Fortran 2008 / OpenMPI
! ==============================================================================
program test_laplace_coupled
  use mod_precision,      only: dp, ip, zero, one, half
  use mod_params,         only: Nx_global, Ny_global, Lx, Ly, rho1, rho2, &
                                mu1, mu2, sigma, gx, gy, cfl_target, &
                                tol_bicgstab2, max_iter_bicgstab2
  use mod_mpi,            only: my_rank, n_procs, dims, nx_loc, ny_loc, &
                                i_start, j_start, dx_grid, dy_grid, comm_cart, &
                                mpi_init_env, mpi_finalize_env, &
                                mpi_setup_cartesian_mesh, mpi_exchange_halos_scalar
  use mod_viscosity,      only: compute_mixture_densities, compute_mixture_viscosities
  use mod_operators,      only: op_divergence
  use mod_solver_coupled, only: compute_adaptive_dt, advance_coupled_step
  use mpi
  implicit none

  integer(ip), parameter :: nh = 2
  real(dp), allocatable :: u_vel(:,:), v_vel(:,:), p_field(:,:), c_field(:,:)
  real(dp), allocatable :: u_prev(:,:), v_prev(:,:)
  real(dp), allocatable :: rho_c(:,:), rho_u(:,:), rho_v(:,:)
  real(dp), allocatable :: kappa(:,:), zeta(:,:), eta(:,:)
  real(dp), allocatable :: kappa_curv(:,:), f_sigma_u(:,:), f_sigma_v(:,:)
  real(dp), allocatable :: div_field(:,:)

  real(dp) :: x_c, y_c, r_dist, R_drop, dt_step, dt_max, res_bicg
  real(dp) :: mass_init_loc, mass_init_glob, mass_curr_loc, mass_curr_glob, mass_err
  real(dp) :: p_in_sum_loc, p_out_sum_loc, p_in_sum_glob, p_out_sum_glob
  real(dp) :: p_in_mean, p_out_mean, delta_p_calc, delta_p_exact, err_laplace
  real(dp) :: u_max_loc, u_max_glob, div_max_loc, div_max_glob
  integer(ip) :: n_in_loc, n_out_loc, n_in_glob, n_out_glob
  integer(ip) :: i_loc, j_loc, i_glob, j_glob, step, n_steps, iter_bicg, ierr

  call mpi_init_env()

  Nx_global = 64
  Ny_global = 64
  Lx        = one
  Ly        = one
  rho1      = 1.0_dp
  rho2      = 1.0_dp
  mu1       = 0.1_dp
  mu2       = 0.1_dp
  sigma     = 1.0_dp
  gx        = zero
  gy        = zero
  R_drop    = 0.25_dp
  cfl_target = 0.25_dp
  tol_bicgstab2 = 1.0e-5_dp
  max_iter_bicgstab2 = 150

  call mpi_setup_cartesian_mesh(Nx_global, Ny_global, Lx, Ly, periodic_x=.false., periodic_y=.false.)

  if (my_rank == 0) then
    write(*,'(A)') "================================================================"
    write(*,'(A)') " TEST SEMAINE 5 : Goutte Statique Couplée (NS + VOF + CSF)      "
    write(*,'(A)') "================================================================"
    write(*,'(A,I4,A,I3,A,I3)') " Processus MPI : ", n_procs, " | Grille : ", dims(1), " x ", dims(2)
    write(*,'(A,I4,A,I4,A,F6.3,A,F6.3)') " Maillage : ", Nx_global, " x ", Ny_global, &
                                         " | R0 = ", R_drop, " | sigma = ", sigma
  end if

  allocate(u_vel(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), v_vel(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(u_prev(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), v_prev(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(p_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh), c_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(rho_c(1-nh:nx_loc+nh, 1-nh:ny_loc+nh), rho_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), rho_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(kappa(1-nh:nx_loc+nh, 1-nh:ny_loc+nh), zeta(0-nh:nx_loc+nh, 0-nh:ny_loc+nh), eta(0-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(kappa_curv(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(f_sigma_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), f_sigma_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(div_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))

  u_vel   = zero
  v_vel   = zero
  u_prev  = zero
  v_prev  = zero
  p_field = zero

  ! Initialisation de la goutte circulaire au centre (0.5, 0.5)
  mass_init_loc = zero
  do j_loc = 1-nh, ny_loc+nh
    j_glob = j_start + j_loc - 1
    y_c = (real(j_glob, dp) - half) * dy_grid
    do i_loc = 1-nh, nx_loc+nh
      i_glob = i_start + i_loc - 1
      x_c = (real(i_glob, dp) - half) * dx_grid
      r_dist = sqrt((x_c - 0.5_dp)**2 + (y_c - 0.5_dp)**2)
      if (r_dist <= R_drop) then
        c_field(i_loc, j_loc) = one
      else
        c_field(i_loc, j_loc) = zero
      end if
      if (i_loc >= 1 .and. i_loc <= nx_loc .and. j_loc >= 1 .and. j_loc <= ny_loc) then
        mass_init_loc = mass_init_loc + c_field(i_loc, j_loc) * dx_grid * dy_grid
      end if
    end do
  end do
  call mpi_exchange_halos_scalar(c_field, nx_loc, ny_loc, nh)
  call MPI_Allreduce(mass_init_loc, mass_init_glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)

  call compute_mixture_densities(c_field, rho_c, rho_u, rho_v, nx_loc, ny_loc, nh)
  call compute_mixture_viscosities(c_field, kappa, zeta, eta, nx_loc, ny_loc, nh)

  dt_max  = 5.0e-4_dp
  n_steps = 20

  do step = 1, n_steps
    call compute_adaptive_dt(u_vel, v_vel, dt_max, dt_step, nx_loc, ny_loc, nh)
    call advance_coupled_step(u_vel, v_vel, p_field, c_field, u_prev, v_prev, &
                              rho_c, rho_u, rho_v, kappa, zeta, eta, &
                              kappa_curv, f_sigma_u, f_sigma_v, &
                              dt_step, step, nx_loc, ny_loc, nh, &
                              iter_out=iter_bicg, res_out=res_bicg, use_bicgstab2=.false.)
  end do

  ! Diagnostic 1 : Conservation de la masse
  mass_curr_loc = zero
  do j_loc = 1, ny_loc
    do i_loc = 1, nx_loc
      mass_curr_loc = mass_curr_loc + c_field(i_loc, j_loc) * dx_grid * dy_grid
    end do
  end do
  call MPI_Allreduce(mass_curr_loc, mass_curr_glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)
  mass_err = abs(mass_curr_glob - mass_init_glob) / max(mass_init_glob, 1.0e-14_dp)

  ! Diagnostic 2 : Saut de pression de Laplace (r < 0.7*R0 vs r > 1.3*R0)
  p_in_sum_loc  = zero
  p_out_sum_loc = zero
  n_in_loc      = 0
  n_out_loc     = 0
  do j_loc = 1, ny_loc
    j_glob = j_start + j_loc - 1
    y_c = (real(j_glob, dp) - half) * dy_grid
    do i_loc = 1, nx_loc
      i_glob = i_start + i_loc - 1
      x_c = (real(i_glob, dp) - half) * dx_grid
      r_dist = sqrt((x_c - 0.5_dp)**2 + (y_c - 0.5_dp)**2)
      if (r_dist < 0.7_dp * R_drop) then
        p_in_sum_loc = p_in_sum_loc + p_field(i_loc, j_loc)
        n_in_loc = n_in_loc + 1
      else if (r_dist > 1.3_dp * R_drop) then
        p_out_sum_loc = p_out_sum_loc + p_field(i_loc, j_loc)
        n_out_loc = n_out_loc + 1
      end if
    end do
  end do
  call MPI_Allreduce(p_in_sum_loc,  p_in_sum_glob,  1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)
  call MPI_Allreduce(p_out_sum_loc, p_out_sum_glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)
  call MPI_Allreduce(n_in_loc,  n_in_glob,  1, MPI_INTEGER, MPI_SUM, comm_cart, ierr)
  call MPI_Allreduce(n_out_loc, n_out_glob, 1, MPI_INTEGER, MPI_SUM, comm_cart, ierr)

  p_in_mean    = p_in_sum_glob  / real(max(n_in_glob, 1), dp)
  p_out_mean   = p_out_sum_glob / real(max(n_out_glob, 1), dp)
  delta_p_calc = p_in_mean - p_out_mean
  delta_p_exact = sigma / R_drop
  err_laplace  = abs(delta_p_calc - delta_p_exact) / delta_p_exact * 100.0_dp

  ! Diagnostic 3 : Courants parasites ||u||_inf et incompressibilité ||div u||_inf
  u_max_loc = zero
  do j_loc = 1, ny_loc
    do i_loc = 1, nx_loc
      u_max_loc = max(u_max_loc, abs(u_vel(i_loc, j_loc)), abs(v_vel(i_loc, j_loc)))
    end do
  end do
  call MPI_Allreduce(u_max_loc, u_max_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_cart, ierr)

  call op_divergence(u_vel, v_vel, div_field, nx_loc, ny_loc, nh)
  div_max_loc = zero
  do j_loc = 1, ny_loc
    do i_loc = 1, nx_loc
      div_max_loc = max(div_max_loc, abs(div_field(i_loc, j_loc)))
    end do
  end do
  call MPI_Allreduce(div_max_loc, div_max_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_cart, ierr)

  if (my_rank == 0) then
    write(*,'(A)') "----------------------------------------------------------------"
    write(*,'(A,I3,A,ES10.3,A,I3,A,ES10.3)') &
      " Après ", n_steps, " pas (dt = ", dt_step, ") | Dernier iter = ", iter_bicg, " | Res = ", res_bicg
    write(*,'(A,F10.4,A)') " Saut de pression Laplace exact : ", delta_p_exact, " Pa"
    write(*,'(A,F10.4,A)') " Saut de pression calculé       : ", delta_p_calc,  " Pa"
    write(*,'(A,F8.3,A)')  " Erreur relative Laplace        : ", err_laplace,   " %"
    write(*,'(A,ES11.4,A)') " Courants parasites ||u||_inf   : ", u_max_glob,    " m/s"
    write(*,'(A,ES11.4)')   " Erreur conservation masse VOF  : ", mass_err
    write(*,'(A,ES11.4)')   " Incompressibilité ||div u||_inf: ", div_max_glob
    write(*,'(A)') "----------------------------------------------------------------"
    if (err_laplace < 5.0_dp .and. mass_err < 1.0e-6_dp .and. div_max_glob < 1.0e-4_dp) then
      write(*,'(A)') " --> [SUCCÈS] Couplage complet NS + VOF + CSF validé !"
    else
      write(*,'(A)') " --> [ATTENTION] Vérifier le saut de pression ou la conservation."
    end if
    write(*,'(A)') "================================================================"
  end if

  deallocate(u_vel, v_vel, u_prev, v_prev, p_field, c_field)
  deallocate(rho_c, rho_u, rho_v, kappa, zeta, eta)
  deallocate(kappa_curv, f_sigma_u, f_sigma_v, div_field)

  call mpi_finalize_env()
end program test_laplace_coupled
