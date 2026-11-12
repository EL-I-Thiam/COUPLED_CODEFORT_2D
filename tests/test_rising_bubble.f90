! ==============================================================================
! PROGRAMME DE TEST : test_rising_bubble.f90
! Benchmark diphasique gravitaire + capillaire (Hysing et al., 2009 - Cas 1)
! Bulle légère ascendante en domaine [0, 1] x [0, 2]
! Projet Solveur 2D Tout-Couplé (FUGU-2D) - Fortran 2008 / OpenMPI
! ==============================================================================
program test_rising_bubble
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
  use mod_io,             only: write_vtk_fields
  use mpi
  implicit none

  integer(ip), parameter :: nh = 2
  real(dp), allocatable :: u_vel(:,:), v_vel(:,:), p_field(:,:), c_field(:,:)
  real(dp), allocatable :: u_prev(:,:), v_prev(:,:)
  real(dp), allocatable :: rho_c(:,:), rho_u(:,:), rho_v(:,:)
  real(dp), allocatable :: kappa(:,:), zeta(:,:), eta(:,:)
  real(dp), allocatable :: kappa_curv(:,:), f_sigma_u(:,:), f_sigma_v(:,:)
  real(dp), allocatable :: div_field(:,:)

  real(dp) :: x_c, y_c, r_dist, R_bub, dt_step, dt_max, t_curr, t_end, res_bicg
  real(dp) :: vol_bub_loc, vol_bub_glob, vol_init, y_sum_loc, y_sum_glob
  real(dp) :: v_sum_loc, v_sum_glob, y_cm, v_rise, mass_err, c_bub, v_cell
  real(dp) :: div_max_loc, div_max_glob
  integer(ip) :: i_loc, j_loc, i_glob, j_glob, step, max_steps, iter_bicg, ierr
  character(len=128) :: vtk_name

  call mpi_init_env()

  ! Paramètres Hysing et al. (2009) - Cas 1
  Nx_global = 32
  Ny_global = 64
  Lx        = 1.0_dp
  Ly        = 2.0_dp
  rho1      = 1000.0_dp   ! Fluide porteur lourd (C = 1)
  rho2      = 100.0_dp    ! Bulle légère (C = 0)
  mu1       = 10.0_dp
  mu2       = 1.0_dp
  sigma     = 24.5_dp
  gx        = zero
  gy        = -0.98_dp
  R_bub     = 0.25_dp
  cfl_target = 0.25_dp
  tol_bicgstab2 = 1.0e-5_dp
  max_iter_bicgstab2 = 150

  call mpi_setup_cartesian_mesh(Nx_global, Ny_global, Lx, Ly, periodic_x=.false., periodic_y=.false.)

  if (my_rank == 0) then
    call execute_command_line("mkdir -p output")
    write(*,'(A)') "================================================================"
    write(*,'(A)') " TEST SEMAINE 5 : Bulle Ascendante 2D (Hysing et al. Cas 1)     "
    write(*,'(A)') "================================================================"
    write(*,'(A,I4,A,I3,A,I3)') " Processus MPI : ", n_procs, " | Grille : ", dims(1), " x ", dims(2)
    write(*,'(A,I4,A,I4,A,F5.2,A,F5.2)') " Maillage : ", Nx_global, " x ", Ny_global, &
                                         " | Domaine : ", Lx, " x ", Ly
    write(*,'(A,F7.1,A,F6.1,A,F6.2)') " rho1/rho2 = ", rho1/rho2, " | mu1/mu2 = ", mu1/mu2, &
                                      " | sigma = ", sigma
    write(*,'(A)') "----------------------------------------------------------------"
  end if

  allocate(u_vel(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), v_vel(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(u_prev(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), v_prev(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(p_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh), c_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(rho_c(1-nh:nx_loc+nh, 1-nh:ny_loc+nh), rho_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), rho_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(kappa(1-nh:nx_loc+nh, 1-nh:ny_loc+nh), zeta(0-nh:nx_loc+nh, 0-nh:ny_loc+nh), eta(0-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(kappa_curv(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(f_sigma_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), f_sigma_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(div_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))

  u_vel      = zero
  v_vel      = zero
  u_prev     = zero
  v_prev     = zero
  p_field    = zero
  kappa_curv = zero

  ! Initialisation de la bulle légère (C = 0) centrée en (0.5, 0.5)
  vol_bub_loc = zero
  do j_loc = 1-nh, ny_loc+nh
    j_glob = j_start + j_loc - 1
    y_c = (real(j_glob, dp) - half) * dy_grid
    do i_loc = 1-nh, nx_loc+nh
      i_glob = i_start + i_loc - 1
      x_c = (real(i_glob, dp) - half) * dx_grid
      r_dist = sqrt((x_c - 0.5_dp)**2 + (y_c - 0.5_dp)**2)
      if (r_dist <= R_bub) then
        c_field(i_loc, j_loc) = zero
      else
        c_field(i_loc, j_loc) = one
      end if
      if (i_loc >= 1 .and. i_loc <= nx_loc .and. j_loc >= 1 .and. j_loc <= ny_loc) then
        vol_bub_loc = vol_bub_loc + (one - c_field(i_loc, j_loc)) * dx_grid * dy_grid
      end if
    end do
  end do
  call mpi_exchange_halos_scalar(c_field, nx_loc, ny_loc, nh)
  call MPI_Allreduce(vol_bub_loc, vol_init, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)

  call compute_mixture_densities(c_field, rho_c, rho_u, rho_v, nx_loc, ny_loc, nh)
  call compute_mixture_viscosities(c_field, kappa, zeta, eta, nx_loc, ny_loc, nh)

  call write_vtk_fields("output/rising_bubble_0000.vtk", p_field, c_field, rho_c, &
                        kappa_curv, u_vel, v_vel, nx_loc, ny_loc, nh)

  t_curr    = zero
  t_end     = 0.5_dp
  dt_max    = 2.0e-3_dp
  max_steps = 40

  do step = 1, max_steps
    if (t_curr >= t_end) exit
    call compute_adaptive_dt(u_vel, v_vel, dt_max, dt_step, nx_loc, ny_loc, nh)
    if (t_curr + dt_step > t_end) dt_step = t_end - t_curr

    call advance_coupled_step(u_vel, v_vel, p_field, c_field, u_prev, v_prev, &
                              rho_c, rho_u, rho_v, kappa, zeta, eta, &
                              kappa_curv, f_sigma_u, f_sigma_v, &
                              dt_step, step, nx_loc, ny_loc, nh, &
                              iter_out=iter_bicg, res_out=res_bicg, use_bicgstab2=.false.)
    t_curr = t_curr + dt_step

    ! Suivi du centre de masse et de la vitesse d'ascension de la bulle
    vol_bub_loc = zero
    y_sum_loc   = zero
    v_sum_loc   = zero
    do j_loc = 1, ny_loc
      j_glob = j_start + j_loc - 1
      y_c = (real(j_glob, dp) - half) * dy_grid
      do i_loc = 1, nx_loc
        c_bub  = one - c_field(i_loc, j_loc)
        v_cell = half * (v_vel(i_loc, j_loc-1) + v_vel(i_loc, j_loc))
        vol_bub_loc = vol_bub_loc + c_bub * dx_grid * dy_grid
        y_sum_loc   = y_sum_loc   + c_bub * y_c * dx_grid * dy_grid
        v_sum_loc   = v_sum_loc   + c_bub * v_cell * dx_grid * dy_grid
      end do
    end do

    call MPI_Allreduce(vol_bub_loc, vol_bub_glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)
    call MPI_Allreduce(y_sum_loc,   y_sum_glob,   1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)
    call MPI_Allreduce(v_sum_loc,   v_sum_glob,   1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)

    y_cm     = y_sum_glob / max(vol_bub_glob, 1.0e-14_dp)
    v_rise   = v_sum_glob / max(vol_bub_glob, 1.0e-14_dp)
    mass_err = abs(vol_bub_glob - vol_init) / max(vol_init, 1.0e-14_dp)

    if (my_rank == 0 .and. (step == 1 .or. mod(step, 10) == 0 .or. step == max_steps)) then
      write(*,'(A,I4,A,F7.4,A,F7.4,A,F8.4,A,I3,A,ES9.2)') &
        " Pas ", step, " | t = ", t_curr, " | y_cm = ", y_cm, &
        " | v_rise = ", v_rise, " | iter = ", iter_bicg, " | res = ", res_bicg
    end if

    if (mod(step, 20) == 0 .or. step == max_steps) then
      write(vtk_name, '(A,I4.4,A)') "output/rising_bubble_", step, ".vtk"
      call write_vtk_fields(trim(vtk_name), p_field, c_field, rho_c, &
                            kappa_curv, u_vel, v_vel, nx_loc, ny_loc, nh)
    end if
  end do

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
    write(*,'(A,F8.5,A)') " Position finale centre de masse y_cm : ", y_cm, " m (init : 0.50000 m)"
    write(*,'(A,F8.5,A)') " Vitesse d'ascension v_rise           : ", v_rise, " m/s"
    write(*,'(A,ES11.4)')  " Erreur conservation volume bulle     : ", mass_err
    write(*,'(A,ES11.4)')  " Incompressibilite finale ||div u||   : ", div_max_glob
    write(*,'(A)') " Fichiers VTK exportés dans output/rising_bubble_*.vtk"
    write(*,'(A)') "----------------------------------------------------------------"
    if (y_cm > 0.5_dp .and. v_rise > zero .and. mass_err < 1.0e-6_dp .and. div_max_glob < 1.0e-4_dp) then
      write(*,'(A)') " --> [SUCCÈS] Ascension de bulle 2D validée (flottabilité + VOF + CSF) !"
    else
      write(*,'(A)') " --> [ATTENTION] Vérifier la dynamique d'ascension ou la conservation."
    end if
    write(*,'(A)') "================================================================"
  end if

  deallocate(u_vel, v_vel, u_prev, v_prev, p_field, c_field)
  deallocate(rho_c, rho_u, rho_v, kappa, zeta, eta)
  deallocate(kappa_curv, f_sigma_u, f_sigma_v, div_field)

  call mpi_finalize_env()
end program test_rising_bubble
