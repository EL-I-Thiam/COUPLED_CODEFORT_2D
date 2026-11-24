! ==============================================================================
! PROGRAMME DE TEST : test_week6_advanced.f90
! Validation des fonctionnalités avancées de la Semaine 6 :
!   1. Pénalisation de Darcy (cylindre solide immergé, alpha = 1e12) + BDF2
!   2. Transport consistant masse-quantité de mouvement (AMP) à fort contraste
!      (Hysing Cas 2 : rho1/rho2 = 1000, mu1/mu2 = 100)
!   3. Sauvegarde et reprise binaire exacte (Checkpoint / Restart)
! Projet Solveur 2D Tout-Couplé (FUGU-2D) - Fortran 2008 / OpenMPI
! ==============================================================================
program test_week6_advanced
  use mod_precision,      only: dp, ip, zero, one, half
  use mod_params,         only: Nx_global, Ny_global, Lx, Ly, rho1, rho2, &
                                mu1, mu2, sigma, gx, gy, time_order, &
                                alpha_penalty, tol_bicgstab2, max_iter_bicgstab2
  use mod_mpi,            only: my_rank, nx_loc, ny_loc, &
                                i_start, j_start, dx_grid, dy_grid, comm_cart, &
                                mpi_init_env, mpi_finalize_env, &
                                mpi_setup_cartesian_mesh, mpi_exchange_halos_scalar
  use mod_viscosity,      only: compute_mixture_densities, compute_mixture_viscosities
  use mod_operators,      only: op_divergence
  use mod_solver_coupled, only: build_cylinder_mask, advance_coupled_step
  use mod_io,             only: write_checkpoint, read_checkpoint
  use mpi
  implicit none

  integer(ip), parameter :: nh = 2
  real(dp), allocatable :: u_vel(:,:), v_vel(:,:), u_prev(:,:), v_prev(:,:)
  real(dp), allocatable :: p_field(:,:), c_field(:,:)
  real(dp), allocatable :: rho_c(:,:), rho_u(:,:), rho_v(:,:)
  real(dp), allocatable :: kappa(:,:), zeta(:,:), eta(:,:)
  real(dp), allocatable :: kappa_curv(:,:), f_sigma_u(:,:), f_sigma_v(:,:)
  real(dp), allocatable :: phi_s_u(:,:), phi_s_v(:,:), div_field(:,:)
  real(dp), allocatable :: u_ref(:,:), v_ref(:,:)

  real(dp) :: dt_step, res_sol, u_sol_max_loc, u_sol_max_glob
  real(dp) :: div_max_loc, div_max_glob, t_chk, diff_loc, diff_glob
  real(dp) :: x_c, y_c, r_dist, vol_init, vol_curr_loc, vol_curr_glob, mass_err
  real(dp) :: y_sum_loc, y_sum_glob, y_cm
  integer(ip) :: i_loc, j_loc, i_glob, j_glob, step, step_chk, iter_sol, ierr
  logical :: ok_chk, pass_pen, pass_amp, pass_chk

  call mpi_init_env()

  if (my_rank == 0) then
    call execute_command_line("mkdir -p output")
    write(*,'(A)') "================================================================"
    write(*,'(A)') " TEST SEMAINE 6 : BDF2, Pénalisation Darcy, AMP & Checkpoint    "
    write(*,'(A)') "================================================================"
  end if

  ! =========================================================================
  ! PARTIE 1 : Pénalisation de Darcy (Cylindre immergé) + BDF2
  ! =========================================================================
  Nx_global = 32
  Ny_global = 32
  Lx = one
  Ly = one
  rho1 = 1.0_dp
  rho2 = 1.0_dp
  mu1  = 0.1_dp
  mu2  = 0.1_dp
  sigma = zero
  gx = 1.0_dp
  gy = zero
  time_order = 2
  alpha_penalty = 1.0e12_dp
  tol_bicgstab2 = 1.0e-5_dp
  max_iter_bicgstab2 = 150

  call mpi_setup_cartesian_mesh(Nx_global, Ny_global, Lx, Ly, periodic_x=.true., periodic_y=.false.)

  allocate(u_vel(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), v_vel(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(u_prev(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), v_prev(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(p_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh), c_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(rho_c(1-nh:nx_loc+nh, 1-nh:ny_loc+nh), rho_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), rho_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(kappa(1-nh:nx_loc+nh, 1-nh:ny_loc+nh), zeta(0-nh:nx_loc+nh, 0-nh:ny_loc+nh), eta(0-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(kappa_curv(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(f_sigma_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), f_sigma_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(phi_s_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), phi_s_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(div_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))

  u_vel = zero; v_vel = zero; u_prev = zero; v_prev = zero
  p_field = zero; c_field = one
  call compute_mixture_densities(c_field, rho_c, rho_u, rho_v, nx_loc, ny_loc, nh)
  call compute_mixture_viscosities(c_field, kappa, zeta, eta, nx_loc, ny_loc, nh)
  call build_cylinder_mask(0.5_dp, 0.5_dp, 0.15_dp, phi_s_u, phi_s_v, nx_loc, ny_loc, nh)

  dt_step = 0.01_dp
  do step = 1, 10
    call advance_coupled_step(u_vel, v_vel, p_field, c_field, u_prev, v_prev, &
                              rho_c, rho_u, rho_v, kappa, zeta, eta, &
                              kappa_curv, f_sigma_u, f_sigma_v, &
                              dt_step, step, nx_loc, ny_loc, nh, &
                              iter_out=iter_sol, res_out=res_sol, use_bicgstab2=.false., &
                              phi_s_u=phi_s_u, phi_s_v=phi_s_v)
  end do

  u_sol_max_loc = zero
  do j_loc = 1, ny_loc
    do i_loc = 1, nx_loc
      if (phi_s_u(i_loc, j_loc) > half) u_sol_max_loc = max(u_sol_max_loc, abs(u_vel(i_loc, j_loc)))
      if (phi_s_v(i_loc, j_loc) > half) u_sol_max_loc = max(u_sol_max_loc, abs(v_vel(i_loc, j_loc)))
    end do
  end do
  call MPI_Allreduce(u_sol_max_loc, u_sol_max_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_cart, ierr)

  call op_divergence(u_vel, v_vel, div_field, nx_loc, ny_loc, nh)
  div_max_loc = zero
  do j_loc = 1, ny_loc
    do i_loc = 1, nx_loc
      div_max_loc = max(div_max_loc, abs(div_field(i_loc, j_loc)))
    end do
  end do
  call MPI_Allreduce(div_max_loc, div_max_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_cart, ierr)

  pass_pen = (u_sol_max_glob < 1.0e-8_dp .and. div_max_glob < 1.0e-4_dp)
  if (my_rank == 0) then
    write(*,'(A,ES11.4,A,ES11.4)') " [TEST A: Darcy+BDF2] ||u||_solid = ", u_sol_max_glob, &
                                    " | ||div u||_inf = ", div_max_glob
  end if

  ! =========================================================================
  ! PARTIE 2 : Checkpoint / Restart bit-a-bit
  ! =========================================================================
  allocate(u_ref(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), v_ref(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  call write_checkpoint("output/test_chk", 10_ip, 0.1_dp, u_vel, v_vel, &
                        u_prev, v_prev, p_field, c_field, nx_loc, ny_loc, nh)

  ! Avancement de référence de 2 pas supplémentaires
  do step = 11, 12
    call advance_coupled_step(u_vel, v_vel, p_field, c_field, u_prev, v_prev, &
                              rho_c, rho_u, rho_v, kappa, zeta, eta, &
                              kappa_curv, f_sigma_u, f_sigma_v, &
                              dt_step, step, nx_loc, ny_loc, nh, &
                              use_bicgstab2=.false., phi_s_u=phi_s_u, phi_s_v=phi_s_v)
  end do
  u_ref = u_vel; v_ref = v_vel

  ! Reprise depuis le checkpoint et ré-exécution des 2 pas
  call read_checkpoint("output/test_chk", step_chk, t_chk, u_vel, v_vel, &
                       u_prev, v_prev, p_field, c_field, nx_loc, ny_loc, nh, ok_chk)
  do step = step_chk + 1, step_chk + 2
    call advance_coupled_step(u_vel, v_vel, p_field, c_field, u_prev, v_prev, &
                              rho_c, rho_u, rho_v, kappa, zeta, eta, &
                              kappa_curv, f_sigma_u, f_sigma_v, &
                              dt_step, step, nx_loc, ny_loc, nh, &
                              use_bicgstab2=.false., phi_s_u=phi_s_u, phi_s_v=phi_s_v)
  end do

  diff_loc = max(maxval(abs(u_vel - u_ref)), maxval(abs(v_vel - v_ref)))
  call MPI_Allreduce(diff_loc, diff_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_cart, ierr)
  pass_chk = (ok_chk .and. diff_glob < 1.0e-14_dp)
  if (my_rank == 0) then
    write(*,'(A,L2,A,ES11.4)') " [TEST B: Checkpoint] Lecture OK = ", ok_chk, &
                               " | Ecart reprise = ", diff_glob
  end if

  deallocate(u_vel, v_vel, u_prev, v_prev, p_field, c_field, u_ref, v_ref)
  deallocate(rho_c, rho_u, rho_v, kappa, zeta, eta)
  deallocate(kappa_curv, f_sigma_u, f_sigma_v, phi_s_u, phi_s_v, div_field)

  ! =========================================================================
  ! PARTIE 3 : Hysing Cas 2 (rho1/rho2 = 1000, mu1/mu2 = 100) + AMP + BDF2
  ! =========================================================================
  Nx_global = 32
  Ny_global = 64
  Lx = 1.0_dp
  Ly = 2.0_dp
  rho1 = 1000.0_dp
  rho2 = 1.0_dp
  mu1  = 10.0_dp
  mu2  = 0.1_dp
  sigma = 1.96_dp
  gx = zero
  gy = -0.98_dp
  time_order = 2

  call mpi_setup_cartesian_mesh(Nx_global, Ny_global, Lx, Ly, periodic_x=.false., periodic_y=.false.)

  allocate(u_vel(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), v_vel(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(u_prev(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), v_prev(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(p_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh), c_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(rho_c(1-nh:nx_loc+nh, 1-nh:ny_loc+nh), rho_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), rho_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(kappa(1-nh:nx_loc+nh, 1-nh:ny_loc+nh), zeta(0-nh:nx_loc+nh, 0-nh:ny_loc+nh), eta(0-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(kappa_curv(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(f_sigma_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), f_sigma_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))

  u_vel = zero; v_vel = zero; u_prev = zero; v_prev = zero; p_field = zero

  vol_curr_loc = zero
  do j_loc = 1-nh, ny_loc+nh
    j_glob = j_start + j_loc - 1
    y_c = (real(j_glob, dp) - half) * dy_grid
    do i_loc = 1-nh, nx_loc+nh
      i_glob = i_start + i_loc - 1
      x_c = (real(i_glob, dp) - half) * dx_grid
      r_dist = sqrt((x_c - 0.5_dp)**2 + (y_c - 0.5_dp)**2)
      if (r_dist <= 0.25_dp) then
        c_field(i_loc, j_loc) = zero
      else
        c_field(i_loc, j_loc) = one
      end if
      if (i_loc >= 1 .and. i_loc <= nx_loc .and. j_loc >= 1 .and. j_loc <= ny_loc) then
        vol_curr_loc = vol_curr_loc + (one - c_field(i_loc, j_loc)) * dx_grid * dy_grid
      end if
    end do
  end do
  call mpi_exchange_halos_scalar(c_field, nx_loc, ny_loc, nh)
  call MPI_Allreduce(vol_curr_loc, vol_init, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)

  dt_step = 1.0e-3_dp
  do step = 1, 20
    call advance_coupled_step(u_vel, v_vel, p_field, c_field, u_prev, v_prev, &
                              rho_c, rho_u, rho_v, kappa, zeta, eta, &
                              kappa_curv, f_sigma_u, f_sigma_v, &
                              dt_step, step, nx_loc, ny_loc, nh, &
                              iter_out=iter_sol, res_out=res_sol, &
                              use_bicgstab2=.false., consistent_amp=.true.)
  end do

  vol_curr_loc = zero
  y_sum_loc    = zero
  do j_loc = 1, ny_loc
    j_glob = j_start + j_loc - 1
    y_c = (real(j_glob, dp) - half) * dy_grid
    do i_loc = 1, nx_loc
      vol_curr_loc = vol_curr_loc + (one - c_field(i_loc, j_loc)) * dx_grid * dy_grid
      y_sum_loc    = y_sum_loc    + (one - c_field(i_loc, j_loc)) * y_c * dx_grid * dy_grid
    end do
  end do
  call MPI_Allreduce(vol_curr_loc, vol_curr_glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)
  call MPI_Allreduce(y_sum_loc,    y_sum_glob,    1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)
  y_cm     = y_sum_glob / max(vol_curr_glob, 1.0e-14_dp)
  mass_err = abs(vol_curr_glob - vol_init) / max(vol_init, 1.0e-14_dp)

  pass_amp = (y_cm > 0.5_dp .and. mass_err < 1.0e-5_dp)
  if (my_rank == 0) then
    write(*,'(A,F8.5,A,ES11.4,A,I3)') " [TEST C: Hysing2+AMP] y_cm = ", y_cm, &
                                       " | err_masse = ", mass_err, " | iter = ", iter_sol
    write(*,'(A)') "----------------------------------------------------------------"
    if (pass_pen .and. pass_chk .and. pass_amp) then
      write(*,'(A)') " --> [SUCCÈS] Toutes les fonctionnalités avancées Semaine 6 validées !"
    else
      write(*,'(A)') " --> [ATTENTION] Un des sous-tests Semaine 6 nécessite vérification."
    end if
    write(*,'(A)') "================================================================"
  end if

  deallocate(u_vel, v_vel, u_prev, v_prev, p_field, c_field)
  deallocate(rho_c, rho_u, rho_v, kappa, zeta, eta)
  deallocate(kappa_curv, f_sigma_u, f_sigma_v)

  call mpi_finalize_env()
end program test_week6_advanced
