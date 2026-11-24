! ==============================================================================
! MODULE : mod_solver_coupled
! Solveur 2D Tout-Couplé (FUGU-2D) - Fortran 2008 / OpenMPI
!
! Moteur d'intégration temporelle diphasique complet :
!   1. Calcul du pas de temps adaptatif (CFL convective + contrainte capillaire)
!   2. Advection conservative de l'interface VOF (Van Leer TVD)
!   3. Transport consistant masse-quantité de mouvement (AMP) en option
!   4. Mise à jour des propriétés effectives (rho, kappa, zeta, eta)
!   5. Lissage SVOF Helmholtz + courbure + force capillaire CSF
!   6. Schéma temporel Euler ordre 1 ou BDF2/Gear ordre 2 (time_order = 2)
!   7. Pénalisation de Darcy pour obstacles solides immergés (alpha_penalty)
!   8. Résolution monolithique Navier-Stokes par BiCGStab / BiCGStab(2) + PCD
! ==============================================================================
module mod_solver_coupled
  use mod_precision,       only: dp, ip, zero, one, half, two, pi
  use mod_params,          only: rho1, rho2, sigma, gx, gy, cfl_target, &
                                 tol_bicgstab2, max_iter_bicgstab2, &
                                 time_order, alpha_penalty
  use mod_mpi,             only: dx_grid, dy_grid, coords, dims, periods, &
                                 i_start, j_start, &
                                 mpi_exchange_halos_scalar, mpi_exchange_halos_u, &
                                 mpi_exchange_halos_v, mpi_global_min_dt
  use mod_viscosity,       only: compute_mixture_densities, compute_mixture_viscosities
  use mod_operators,       only: apply_velocity_bcs
  use mod_vof,             only: vof_advection_tvd, svof_helmholtz_filter
  use mod_surface_tension, only: compute_interface_curvature, compute_csf_force
  use mod_bicgstab2,       only: bicgstab_solve, bicgstab2_solve
  implicit none
  private

  public :: compute_adaptive_dt
  public :: build_cylinder_mask
  public :: advance_coupled_step

contains

  ! ----------------------------------------------------------------------------
  ! 1. Calcul du pas de temps adaptatif global (CFL convective + Capillaire)
  ! ----------------------------------------------------------------------------
  subroutine compute_adaptive_dt(u_vel, v_vel, dt_max, dt_out, nx, ny, nh)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)  :: u_vel
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)  :: v_vel
    real(dp), intent(in)  :: dt_max
    real(dp), intent(out) :: dt_out

    real(dp) :: u_max_loc, v_max_loc, h_min
    real(dp) :: dt_cfl_loc, dt_cap, dt_grav, g_mag, dt_loc
    integer(ip) :: i, j

    h_min = min(dx_grid, dy_grid)

    u_max_loc = zero
    do j = 1, ny
      do i = 0, nx
        u_max_loc = max(u_max_loc, abs(u_vel(i, j)))
      end do
    end do

    v_max_loc = zero
    do j = 0, ny
      do i = 1, nx
        v_max_loc = max(v_max_loc, abs(v_vel(i, j)))
      end do
    end do

    ! Contrainte CFL convective
    if (u_max_loc / dx_grid + v_max_loc / dy_grid > 1.0e-12_dp) then
      dt_cfl_loc = cfl_target / (u_max_loc / dx_grid + v_max_loc / dy_grid)
    else
      dt_cfl_loc = dt_max
    end if

    dt_loc = min(dt_max, dt_cfl_loc)

    ! Contrainte capillaire de Brackbill et al. (1992)
    if (sigma > 1.0e-14_dp) then
      dt_cap = 0.5_dp * sqrt((rho1 + rho2) * (h_min**3) / (4.0_dp * pi * sigma))
      dt_loc = min(dt_loc, dt_cap)
    end if

    ! Contrainte gravitaire
    g_mag = sqrt(gx**2 + gy**2)
    if (g_mag > 1.0e-14_dp) then
      dt_grav = 0.5_dp * sqrt(h_min / g_mag)
      dt_loc = min(dt_loc, dt_grav)
    end if

    dt_out = mpi_global_min_dt(dt_loc)
  end subroutine compute_adaptive_dt

  ! ----------------------------------------------------------------------------
  ! 2. Construction du masque solide (0 ou 1) d'un cylindre immergé sur faces MAC
  ! ----------------------------------------------------------------------------
  subroutine build_cylinder_mask(xc, yc, rad, phi_s_u, phi_s_v, nx, ny, nh)
    real(dp), intent(in) :: xc, yc, rad
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(out) :: phi_s_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(out) :: phi_s_v

    integer(ip) :: i_loc, j_loc, i_glob, j_glob
    real(dp) :: x_pos, y_pos, r_dist

    phi_s_u = zero
    phi_s_v = zero

    do j_loc = 1, ny
      j_glob = j_start + j_loc - 1
      y_pos = (real(j_glob, dp) - half) * dy_grid
      do i_loc = 0, nx
        i_glob = i_start + i_loc - 1
        x_pos = real(i_glob, dp) * dx_grid
        r_dist = sqrt((x_pos - xc)**2 + (y_pos - yc)**2)
        if (r_dist <= rad) phi_s_u(i_loc, j_loc) = one
      end do
    end do

    do j_loc = 0, ny
      j_glob = j_start + j_loc - 1
      y_pos = real(j_glob, dp) * dy_grid
      do i_loc = 1, nx
        i_glob = i_start + i_loc - 1
        x_pos = (real(i_glob, dp) - half) * dx_grid
        r_dist = sqrt((x_pos - xc)**2 + (y_pos - yc)**2)
        if (r_dist <= rad) phi_s_v(i_loc, j_loc) = one
      end do
    end do

    call mpi_exchange_halos_u(phi_s_u, nx, ny, nh)
    call mpi_exchange_halos_v(phi_s_v, nx, ny, nh)
  end subroutine build_cylinder_mask

  ! ----------------------------------------------------------------------------
  ! 3. Avancement d'un pas de temps complet du solveur diphasique tout-couplé
  ! ----------------------------------------------------------------------------
  subroutine advance_coupled_step(u_vel, v_vel, p_field, c_field, &
                                  u_prev, v_prev, &
                                  rho_c, rho_u, rho_v, kappa, zeta, eta, &
                                  kappa_curv, f_sigma_u, f_sigma_v, &
                                  dt_step, step_idx, nx, ny, nh, &
                                  iter_out, res_out, use_bicgstab2, u_top_lid, &
                                  tol_in, max_iter_in, &
                                  consistent_amp, phi_s_u, phi_s_v, u_solid_in, v_solid_in)
    integer(ip), intent(in) :: step_idx, nx, ny, nh
    real(dp), intent(in) :: dt_step
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: u_vel, u_prev
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(inout) :: v_vel, v_prev
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: p_field, c_field
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: rho_c, kappa, kappa_curv
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: rho_u, f_sigma_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(inout) :: rho_v, f_sigma_v
    real(dp), dimension(0-nh:nx+nh, 0-nh:ny+nh), intent(inout) :: zeta, eta
    integer(ip), optional, intent(out) :: iter_out
    real(dp), optional, intent(out)    :: res_out
    logical, optional, intent(in)      :: use_bicgstab2, consistent_amp
    real(dp), optional, intent(in)     :: u_top_lid, tol_in, u_solid_in, v_solid_in
    integer(ip), optional, intent(in)  :: max_iter_in
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), optional, intent(in) :: phi_s_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), optional, intent(in) :: phi_s_v

    real(dp), allocatable :: u_n(:,:), v_n(:,:), u_nm1(:,:), v_nm1(:,:)
    real(dp), allocatable :: u_adv(:,:), v_adv(:,:)
    real(dp), allocatable :: c_new(:,:), c_smooth(:,:)
    real(dp), allocatable :: flux_c_x(:,:), flux_c_y(:,:)
    real(dp), allocatable :: mdot_u_cons(:,:), mdot_v_cons(:,:)
    real(dp), allocatable :: pen_u(:,:), pen_v(:,:)
    real(dp), allocatable :: B_u(:,:), B_v(:,:), B_p(:,:)

    logical :: use_l2, use_amp, has_penalty, use_bdf2
    real(dp) :: tol_sol, u_lid, delta_filter, res_loc, dt_eff
    real(dp) :: u_sol_val, v_sol_val
    integer(ip) :: max_it, it_loc, i, j, i_u_max, j_v_max

    use_l2 = .true.
    if (present(use_bicgstab2)) use_l2 = use_bicgstab2

    use_amp = .false.
    if (present(consistent_amp)) use_amp = consistent_amp

    has_penalty = present(phi_s_u) .and. present(phi_s_v)
    u_sol_val = zero
    v_sol_val = zero
    if (present(u_solid_in)) u_sol_val = u_solid_in
    if (present(v_solid_in)) v_sol_val = v_solid_in

    tol_sol = tol_bicgstab2
    if (present(tol_in)) tol_sol = tol_in

    max_it = max_iter_bicgstab2
    if (present(max_iter_in)) max_it = max_iter_in

    u_lid = zero
    if (present(u_top_lid)) u_lid = u_top_lid

    use_bdf2 = (time_order == 2 .and. step_idx >= 2)
    if (use_bdf2) then
      dt_eff = dt_step / 1.5_dp
    else
      dt_eff = dt_step
    end if

    allocate(u_n(0-nh:nx+nh, 1-nh:ny+nh), v_n(1-nh:nx+nh, 0-nh:ny+nh))
    allocate(u_nm1(0-nh:nx+nh, 1-nh:ny+nh), v_nm1(1-nh:nx+nh, 0-nh:ny+nh))
    allocate(u_adv(0-nh:nx+nh, 1-nh:ny+nh), v_adv(1-nh:nx+nh, 0-nh:ny+nh))
    allocate(c_new(1-nh:nx+nh, 1-nh:ny+nh), c_smooth(1-nh:nx+nh, 1-nh:ny+nh))
    allocate(B_u(0-nh:nx+nh, 1-nh:ny+nh), B_v(1-nh:nx+nh, 0-nh:ny+nh), B_p(1-nh:nx+nh, 1-nh:ny+nh))

    ! 1. Sauvegarde de u^n, u^(n-1) et extrapolation d'Adams-Bashforth u_adv = 2 u^n - u^(n-1)
    u_n   = u_vel
    v_n   = v_vel
    u_nm1 = u_prev
    v_nm1 = v_prev
    if (step_idx <= 1) then
      u_adv = u_n
      v_adv = v_n
    else
      u_adv = two * u_n - u_nm1
      v_adv = two * v_n - v_nm1
    end if
    u_prev = u_n
    v_prev = v_n

    ! 2. Transport conservatif VOF (C^n -> C^(n+1)) et flux de masse AMP si demandé
    if (use_amp) then
      allocate(flux_c_x(0-nh:nx+nh, 1-nh:ny+nh), flux_c_y(1-nh:nx+nh, 0-nh:ny+nh))
      allocate(mdot_u_cons(0-nh:nx+nh, 1-nh:ny+nh), mdot_v_cons(1-nh:nx+nh, 0-nh:ny+nh))
      call vof_advection_tvd(c_field, c_new, u_n, v_n, dt_step, nx, ny, nh, &
                             flux_x_out=flux_c_x, flux_y_out=flux_c_y)
      mdot_u_cons = (rho1 - rho2) * flux_c_x + rho2 * u_adv
      mdot_v_cons = (rho1 - rho2) * flux_c_y + rho2 * v_adv
    else
      call vof_advection_tvd(c_field, c_new, u_n, v_n, dt_step, nx, ny, nh)
    end if
    c_field = c_new

    ! 3. Mise à jour des propriétés de mélange à t^(n+1)
    call compute_mixture_densities(c_field, rho_c, rho_u, rho_v, nx, ny, nh)
    call compute_mixture_viscosities(c_field, kappa, zeta, eta, nx, ny, nh)

    ! 4. Lissage SVOF et force capillaire CSF
    if (sigma > 1.0e-14_dp) then
      delta_filter = 1.5_dp * min(dx_grid, dy_grid)
      call svof_helmholtz_filter(c_field, c_smooth, delta_filter, nx, ny, nh, &
                                 max_iter=150, tol=1.0e-7_dp)
      call compute_interface_curvature(c_smooth, kappa_curv, nx, ny, nh)
      call compute_csf_force(c_field, kappa_curv, f_sigma_u, f_sigma_v, nx, ny, nh)
    else
      kappa_curv = zero
      f_sigma_u  = zero
      f_sigma_v  = zero
    end if

    ! 5. Préparation de la pénalisation de Darcy si présente
    if (has_penalty) then
      allocate(pen_u(0-nh:nx+nh, 1-nh:ny+nh), pen_v(1-nh:nx+nh, 0-nh:ny+nh))
      pen_u = alpha_penalty * phi_s_u
      pen_v = alpha_penalty * phi_s_v
    end if

    ! 6. Assemblage du second membre monolithique (B_u, B_v, B_p)
    i_u_max = nx
    if (.not. periods(1) .and. coords(1) == dims(1) - 1) i_u_max = nx - 1

    j_v_max = ny
    if (.not. periods(2) .and. coords(2) == dims(2) - 1) j_v_max = ny - 1

    B_u = zero
    do j = 1, ny
      do i = 1, i_u_max
        if (use_bdf2) then
          B_u(i, j) = (rho_u(i, j) / dt_step) * (two * u_n(i, j) - half * u_nm1(i, j)) + &
                      rho_u(i, j) * gx + f_sigma_u(i, j)
        else
          B_u(i, j) = (rho_u(i, j) / dt_step) * u_n(i, j) + &
                      rho_u(i, j) * gx + f_sigma_u(i, j)
        end if
        if (has_penalty) B_u(i, j) = B_u(i, j) + pen_u(i, j) * u_sol_val
      end do
    end do

    if (present(u_top_lid) .and. .not. periods(2) .and. coords(2) == dims(2) - 1) then
      do i = 1, i_u_max
        B_u(i, ny) = B_u(i, ny) + &
          (two * (zeta(i, ny) - eta(i, ny)) / (dy_grid * dy_grid)) * u_lid
      end do
    end if

    B_v = zero
    do j = 1, j_v_max
      do i = 1, nx
        if (use_bdf2) then
          B_v(i, j) = (rho_v(i, j) / dt_step) * (two * v_n(i, j) - half * v_nm1(i, j)) + &
                      rho_v(i, j) * gy + f_sigma_v(i, j)
        else
          B_v(i, j) = (rho_v(i, j) / dt_step) * v_n(i, j) + &
                      rho_v(i, j) * gy + f_sigma_v(i, j)
        end if
        if (has_penalty) B_v(i, j) = B_v(i, j) + pen_v(i, j) * v_sol_val
      end do
    end do

    B_p = zero

    ! 7. Résolution du système point-selle monolithique
    if (use_l2) then
      if (has_penalty .and. use_amp) then
        call bicgstab2_solve(u_vel, v_vel, p_field, B_u, B_v, B_p, &
                             u_adv, v_adv, rho_c, rho_u, rho_v, &
                             kappa, zeta, eta, dt_eff, nx, ny, nh, &
                             tol_in=tol_sol, max_iter_in=max_it, &
                             iter_out=it_loc, res_out=res_loc, u_top_lid=u_lid, &
                             penalty_u=pen_u, penalty_v=pen_v, &
                             mdot_u_cons=mdot_u_cons, mdot_v_cons=mdot_v_cons)
      else if (has_penalty) then
        call bicgstab2_solve(u_vel, v_vel, p_field, B_u, B_v, B_p, &
                             u_adv, v_adv, rho_c, rho_u, rho_v, &
                             kappa, zeta, eta, dt_eff, nx, ny, nh, &
                             tol_in=tol_sol, max_iter_in=max_it, &
                             iter_out=it_loc, res_out=res_loc, u_top_lid=u_lid, &
                             penalty_u=pen_u, penalty_v=pen_v)
      else if (use_amp) then
        call bicgstab2_solve(u_vel, v_vel, p_field, B_u, B_v, B_p, &
                             u_adv, v_adv, rho_c, rho_u, rho_v, &
                             kappa, zeta, eta, dt_eff, nx, ny, nh, &
                             tol_in=tol_sol, max_iter_in=max_it, &
                             iter_out=it_loc, res_out=res_loc, u_top_lid=u_lid, &
                             mdot_u_cons=mdot_u_cons, mdot_v_cons=mdot_v_cons)
      else
        call bicgstab2_solve(u_vel, v_vel, p_field, B_u, B_v, B_p, &
                             u_adv, v_adv, rho_c, rho_u, rho_v, &
                             kappa, zeta, eta, dt_eff, nx, ny, nh, &
                             tol_in=tol_sol, max_iter_in=max_it, &
                             iter_out=it_loc, res_out=res_loc, u_top_lid=u_lid)
      end if
    else
      if (has_penalty .and. use_amp) then
        call bicgstab_solve(u_vel, v_vel, p_field, B_u, B_v, B_p, &
                            u_adv, v_adv, rho_c, rho_u, rho_v, &
                            kappa, zeta, eta, dt_eff, nx, ny, nh, &
                            tol_in=tol_sol, max_iter_in=max_it, &
                            iter_out=it_loc, res_out=res_loc, u_top_lid=u_lid, &
                            penalty_u=pen_u, penalty_v=pen_v, &
                            mdot_u_cons=mdot_u_cons, mdot_v_cons=mdot_v_cons)
      else if (has_penalty) then
        call bicgstab_solve(u_vel, v_vel, p_field, B_u, B_v, B_p, &
                            u_adv, v_adv, rho_c, rho_u, rho_v, &
                            kappa, zeta, eta, dt_eff, nx, ny, nh, &
                            tol_in=tol_sol, max_iter_in=max_it, &
                            iter_out=it_loc, res_out=res_loc, u_top_lid=u_lid, &
                            penalty_u=pen_u, penalty_v=pen_v)
      else if (use_amp) then
        call bicgstab_solve(u_vel, v_vel, p_field, B_u, B_v, B_p, &
                            u_adv, v_adv, rho_c, rho_u, rho_v, &
                            kappa, zeta, eta, dt_eff, nx, ny, nh, &
                            tol_in=tol_sol, max_iter_in=max_it, &
                            iter_out=it_loc, res_out=res_loc, u_top_lid=u_lid, &
                            mdot_u_cons=mdot_u_cons, mdot_v_cons=mdot_v_cons)
      else
        call bicgstab_solve(u_vel, v_vel, p_field, B_u, B_v, B_p, &
                            u_adv, v_adv, rho_c, rho_u, rho_v, &
                            kappa, zeta, eta, dt_eff, nx, ny, nh, &
                            tol_in=tol_sol, max_iter_in=max_it, &
                            iter_out=it_loc, res_out=res_loc, u_top_lid=u_lid)
      end if
    end if

    ! 8. Synchronisation finale des halos
    call mpi_exchange_halos_u(u_vel, nx, ny, nh)
    call mpi_exchange_halos_v(v_vel, nx, ny, nh)
    call mpi_exchange_halos_scalar(p_field, nx, ny, nh)
    call apply_velocity_bcs(u_vel, v_vel, nx, ny, nh, u_top_lid=u_lid)

    if (present(iter_out)) iter_out = it_loc
    if (present(res_out))  res_out  = res_loc

    deallocate(u_n, v_n, u_nm1, v_nm1, u_adv, v_adv, c_new, c_smooth, B_u, B_v, B_p)
    if (use_amp) deallocate(flux_c_x, flux_c_y, mdot_u_cons, mdot_v_cons)
    if (has_penalty) deallocate(pen_u, pen_v)
  end subroutine advance_coupled_step

end module mod_solver_coupled
