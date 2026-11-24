! ==============================================================================
! PROGRAMME PRINCIPAL : main.f90 (FUGU-2D)
! Solveur 2D Tout-Couplé Navier-Stokes / VOF-SVOF / CSF / Pénalisation Darcy
! Piloté par fichier de configuration Fortran namelist (input.nml)
! avec support sauvegarde/reprise (checkpoint/restart) et export VTK.
! ==============================================================================
program main_coupled
  use mod_precision,      only: dp, ip, zero, one, half
  use mod_params,         only: Lx, Ly, Nx_global, Ny_global, rho1, rho2, &
                                mu1, mu2, sigma, gx, gy, dt, t_final, cfl_target, &
                                max_time_steps, time_order, tol_bicgstab2, &
                                max_iter_bicgstab2, alpha_penalty
  use mod_mpi,            only: my_rank, n_procs, nx_loc, ny_loc, &
                                i_start, j_start, dx_grid, dy_grid, &
                                mpi_init_env, mpi_finalize_env, &
                                mpi_setup_cartesian_mesh, mpi_exchange_halos_scalar
  use mod_viscosity,      only: compute_mixture_densities, compute_mixture_viscosities
  use mod_solver_coupled, only: compute_adaptive_dt, build_cylinder_mask, advance_coupled_step
  use mod_io,             only: write_vtk_fields, write_checkpoint, read_checkpoint
  implicit none

  integer(ip), parameter :: nh = 2

  ! Paramètres namelist
  character(len=32)  :: case_type
  logical            :: periodic_x, periodic_y, adaptive_dt, consistent_amp
  logical            :: use_bicgstab2, restart
  integer(ip)        :: vtk_freq, chk_freq
  real(dp)           :: u_lid, cyl_xc, cyl_yc, cyl_rad, drop_rad
  character(len=128) :: chk_prefix, nml_file

  namelist /sim_params/ case_type, Lx, Ly, Nx_global, Ny_global, &
                        rho1, rho2, mu1, mu2, sigma, gx, gy, &
                        dt, t_final, cfl_target, max_time_steps, time_order, &
                        tol_bicgstab2, max_iter_bicgstab2, alpha_penalty, &
                        periodic_x, periodic_y, adaptive_dt, consistent_amp, &
                        use_bicgstab2, restart, vtk_freq, chk_freq, &
                        u_lid, cyl_xc, cyl_yc, cyl_rad, drop_rad, chk_prefix

  real(dp), allocatable :: u_vel(:,:), v_vel(:,:), u_prev(:,:), v_prev(:,:)
  real(dp), allocatable :: p_field(:,:), c_field(:,:)
  real(dp), allocatable :: rho_c(:,:), rho_u(:,:), rho_v(:,:)
  real(dp), allocatable :: kappa(:,:), zeta(:,:), eta(:,:)
  real(dp), allocatable :: kappa_curv(:,:), f_sigma_u(:,:), f_sigma_v(:,:)
  real(dp), allocatable :: phi_s_u(:,:), phi_s_v(:,:)

  integer(ip) :: unit_nml, step, start_step, iter_sol, i_loc, j_loc, i_glob, j_glob
  real(dp) :: t_curr, dt_step, res_sol, x_c, y_c, r_dist
  logical :: exists_nml, ok_chk, has_cyl
  character(len=256) :: vtk_name

  call mpi_init_env()

  ! Valeurs par défaut
  case_type      = 'rising_bubble'
  periodic_x     = .false.
  periodic_y     = .false.
  adaptive_dt    = .true.
  consistent_amp = .true.
  use_bicgstab2  = .false.
  restart        = .false.
  vtk_freq       = 20
  chk_freq       = 50
  u_lid          = one
  cyl_xc         = 0.5_dp
  cyl_yc         = 0.5_dp
  cyl_rad        = 0.125_dp
  drop_rad       = 0.25_dp
  chk_prefix     = 'output/chk_fugu2d'
  nml_file       = 'input.nml'

  if (command_argument_count() >= 1) then
    call get_command_argument(1, nml_file)
  end if

  inquire(file=trim(nml_file), exist=exists_nml)
  if (exists_nml) then
    open(newunit=unit_nml, file=trim(nml_file), status='old', action='read')
    read(unit_nml, nml=sim_params)
    close(unit_nml)
  end if

  call mpi_setup_cartesian_mesh(Nx_global, Ny_global, Lx, Ly, &
                                periodic_x=periodic_x, periodic_y=periodic_y)

  if (my_rank == 0) then
    call execute_command_line("mkdir -p output")
    write(*,'(A)') "================================================================"
    write(*,'(A)') " FUGU-2D : Solveur Monolithique Tout-Couplé (Navier-Stokes/VOF) "
    write(*,'(A)') "================================================================"
    write(*,'(A,A,A,I3,A,I3,A,I4)') " Cas : ", trim(case_type), &
      " | Grille : ", Nx_global, " x ", Ny_global, " | Procs MPI : ", n_procs
    write(*,'(A,I2,A,L2,A,L2)') " Ordre temporel : BDF", time_order, &
      " | AMP consistant : ", consistent_amp, " | Reprise : ", restart
  end if

  allocate(u_vel(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), v_vel(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(u_prev(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), v_prev(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(p_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh), c_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(rho_c(1-nh:nx_loc+nh, 1-nh:ny_loc+nh), rho_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), rho_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(kappa(1-nh:nx_loc+nh, 1-nh:ny_loc+nh), zeta(0-nh:nx_loc+nh, 0-nh:ny_loc+nh), eta(0-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(kappa_curv(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(f_sigma_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), f_sigma_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))

  u_vel      = zero
  v_vel      = zero
  u_prev     = zero
  v_prev     = zero
  p_field    = zero
  c_field    = one
  kappa_curv = zero
  has_cyl    = (trim(case_type) == 'cylinder_flow')

  if (has_cyl) then
    allocate(phi_s_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh), phi_s_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
    call build_cylinder_mask(cyl_xc, cyl_yc, cyl_rad, phi_s_u, phi_s_v, nx_loc, ny_loc, nh)
  end if

  start_step = 1
  t_curr     = zero

  if (restart) then
    call read_checkpoint(trim(chk_prefix), start_step, t_curr, u_vel, v_vel, &
                         u_prev, v_prev, p_field, c_field, nx_loc, ny_loc, nh, ok_chk)
    if (ok_chk) then
      start_step = start_step + 1
      if (my_rank == 0) write(*,'(A,I6,A,F10.5)') " Reprise depuis le pas ", start_step - 1, " | t = ", t_curr
    else
      if (my_rank == 0) write(*,'(A)') " [INFO] Aucun checkpoint valide trouvé, initialisation à t = 0."
      start_step = 1
      t_curr     = zero
    end if
  end if

  if (.not. restart .or. start_step == 1) then
    select case (trim(case_type))
    case ('rising_bubble')
      do j_loc = 1-nh, ny_loc+nh
        j_glob = j_start + j_loc - 1
        y_c = (real(j_glob, dp) - half) * dy_grid
        do i_loc = 1-nh, nx_loc+nh
          i_glob = i_start + i_loc - 1
          x_c = (real(i_glob, dp) - half) * dx_grid
          r_dist = sqrt((x_c - 0.5_dp * Lx)**2 + (y_c - 0.5_dp)**2)
          if (r_dist <= drop_rad) then
            c_field(i_loc, j_loc) = zero
          else
            c_field(i_loc, j_loc) = one
          end if
        end do
      end do
    case ('laplace')
      do j_loc = 1-nh, ny_loc+nh
        j_glob = j_start + j_loc - 1
        y_c = (real(j_glob, dp) - half) * dy_grid
        do i_loc = 1-nh, nx_loc+nh
          i_glob = i_start + i_loc - 1
          x_c = (real(i_glob, dp) - half) * dx_grid
          r_dist = sqrt((x_c - 0.5_dp * Lx)**2 + (y_c - 0.5_dp * Ly)**2)
          if (r_dist <= drop_rad) then
            c_field(i_loc, j_loc) = one
          else
            c_field(i_loc, j_loc) = zero
          end if
        end do
      end do
    case ('poiseuille')
      do j_loc = 1-nh, ny_loc+nh
        j_glob = j_start + j_loc - 1
        y_c = (real(j_glob, dp) - half) * dy_grid
        do i_loc = 1-nh, nx_loc+nh
          if (y_c <= 0.5_dp * Ly) then
            c_field(i_loc, j_loc) = one
          else
            c_field(i_loc, j_loc) = zero
          end if
        end do
      end do
    case default
      c_field = one
    end select
    call mpi_exchange_halos_scalar(c_field, nx_loc, ny_loc, nh)
  end if

  call compute_mixture_densities(c_field, rho_c, rho_u, rho_v, nx_loc, ny_loc, nh)
  call compute_mixture_viscosities(c_field, kappa, zeta, eta, nx_loc, ny_loc, nh)

  if (vtk_freq > 0 .and. start_step == 1) then
    call write_vtk_fields("output/fugu2d_0000.vtk", p_field, c_field, rho_c, &
                          kappa_curv, u_vel, v_vel, nx_loc, ny_loc, nh)
  end if

  do step = start_step, max_time_steps
    if (t_curr >= t_final) exit

    if (adaptive_dt) then
      call compute_adaptive_dt(u_vel, v_vel, dt, dt_step, nx_loc, ny_loc, nh)
    else
      dt_step = dt
    end if
    if (t_curr + dt_step > t_final) dt_step = t_final - t_curr

    if (trim(case_type) == 'cavity') then
      call advance_coupled_step(u_vel, v_vel, p_field, c_field, u_prev, v_prev, &
                                rho_c, rho_u, rho_v, kappa, zeta, eta, &
                                kappa_curv, f_sigma_u, f_sigma_v, &
                                dt_step, step, nx_loc, ny_loc, nh, &
                                iter_out=iter_sol, res_out=res_sol, &
                                use_bicgstab2=use_bicgstab2, u_top_lid=u_lid, &
                                consistent_amp=consistent_amp)
    else if (has_cyl) then
      call advance_coupled_step(u_vel, v_vel, p_field, c_field, u_prev, v_prev, &
                                rho_c, rho_u, rho_v, kappa, zeta, eta, &
                                kappa_curv, f_sigma_u, f_sigma_v, &
                                dt_step, step, nx_loc, ny_loc, nh, &
                                iter_out=iter_sol, res_out=res_sol, &
                                use_bicgstab2=use_bicgstab2, &
                                consistent_amp=consistent_amp, &
                                phi_s_u=phi_s_u, phi_s_v=phi_s_v)
    else
      call advance_coupled_step(u_vel, v_vel, p_field, c_field, u_prev, v_prev, &
                                rho_c, rho_u, rho_v, kappa, zeta, eta, &
                                kappa_curv, f_sigma_u, f_sigma_v, &
                                dt_step, step, nx_loc, ny_loc, nh, &
                                iter_out=iter_sol, res_out=res_sol, &
                                use_bicgstab2=use_bicgstab2, &
                                consistent_amp=consistent_amp)
    end if

    t_curr = t_curr + dt_step

    if (my_rank == 0 .and. (step == start_step .or. mod(step, 10) == 0)) then
      write(*,'(A,I5,A,F8.4,A,ES9.2,A,I3,A,ES9.2)') &
        " Pas ", step, " | t = ", t_curr, " | dt = ", dt_step, &
        " | iter = ", iter_sol, " | res = ", res_sol
    end if

    if (vtk_freq > 0 .and. mod(step, vtk_freq) == 0) then
      write(vtk_name, '(A,I4.4,A)') "output/fugu2d_", step, ".vtk"
      call write_vtk_fields(trim(vtk_name), p_field, c_field, rho_c, &
                            kappa_curv, u_vel, v_vel, nx_loc, ny_loc, nh)
    end if

    if (chk_freq > 0 .and. mod(step, chk_freq) == 0) then
      call write_checkpoint(trim(chk_prefix), step, t_curr, u_vel, v_vel, &
                            u_prev, v_prev, p_field, c_field, nx_loc, ny_loc, nh)
    end if
  end do

  if (my_rank == 0) then
    write(*,'(A)') "================================================================"
    write(*,'(A,F10.5)') " Simulation terminée à t = ", t_curr
    write(*,'(A)') "================================================================"
  end if

  deallocate(u_vel, v_vel, u_prev, v_prev, p_field, c_field)
  deallocate(rho_c, rho_u, rho_v, kappa, zeta, eta)
  deallocate(kappa_curv, f_sigma_u, f_sigma_v)
  if (has_cyl) deallocate(phi_s_u, phi_s_v)

  call mpi_finalize_env()
end program main_coupled
