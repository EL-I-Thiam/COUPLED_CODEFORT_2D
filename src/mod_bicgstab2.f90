! ==============================================================================
! MODULE : mod_bicgstab2
! Solveur 2D Tout-Couple (FUGU-2D) - Fortran 2008 / OpenMPI
!
! Famille de méthodes de Krylov préconditionnées à gauche pour le système selle 2D :
!   1. bicgstab_solve  : BiCGStab standard (ordre l = 1, van der Vorst 1992)
!   2. bicgstab2_solve : BiCGStab(2) (ordre l = 2, Sleijpen & van der Vorst 1994,
!                        Listing 4.2 de la thèse de Mohamed El Ouafa)
! ==============================================================================
module mod_bicgstab2
  use mod_precision,  only: dp, ip, zero, one, two
  use mod_mpi,        only: dx_grid, dy_grid, comm_cart, my_rank, &
                            mpi_exchange_halos_scalar, mpi_exchange_halos_u, &
                            mpi_exchange_halos_v
  use mod_operators,  only: matvec_coupled, apply_velocity_bcs
  use mod_pcd_precond, only: apply_pcd_precond
  use mpi
  implicit none
  private

  public :: bicgstab_solve
  public :: bicgstab2_solve
  public :: coupled_dot

contains

  ! ----------------------------------------------------------------------------
  ! Produit scalaire global L2 sur l'espace produit décalé H = U x V x P
  ! (u, v) = sum( u_i * v_i * vol ) avec réduction MPI
  ! ----------------------------------------------------------------------------
  function coupled_dot(u1, v1, p1, u2, v2, p2, nx, ny, nh) result(dot_val)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in) :: u1, u2
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in) :: v1, v2
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in) :: p1, p2
    real(dp) :: dot_val, loc_dot
    real(dp) :: vol
    integer(ip) :: i, j, ierr

    vol = dx_grid * dy_grid
    loc_dot = zero

    ! Faces u internes et frontières physiques
    do j = 1, ny
      do i = 0, nx
        loc_dot = loc_dot + u1(i, j) * u2(i, j) * vol
      end do
    end do

    ! Faces v internes et frontières physiques
    do j = 0, ny
      do i = 1, nx
        loc_dot = loc_dot + v1(i, j) * v2(i, j) * vol
      end do
    end do

    ! Centres pour la pression
    do j = 1, ny
      do i = 1, nx
        loc_dot = loc_dot + p1(i, j) * p2(i, j) * vol
      end do
    end do

    call MPI_Allreduce(loc_dot, dot_val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)
  end function coupled_dot

  ! ----------------------------------------------------------------------------
  ! 1. Solveur BiCGStab standard préconditionné à gauche (ordre l = 1)
  ! ----------------------------------------------------------------------------
  subroutine bicgstab_solve(X_u, X_v, X_p, B_u, B_v, B_p, &
                            u_adv, v_adv, rho_c, rho_u, rho_v, &
                            kappa, zeta, eta, dt_step, nx, ny, nh, &
                            tol_in, max_iter_in, iter_out, res_out, u_top_lid)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: X_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(inout) :: X_v
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: X_p
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: B_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)    :: B_v
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: B_p
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: u_adv
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)    :: v_adv
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: rho_c
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: rho_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)    :: rho_v
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: kappa
    real(dp), dimension(0-nh:nx+nh, 0-nh:ny+nh), intent(in)    :: zeta, eta
    real(dp), intent(in) :: dt_step
    real(dp), optional, intent(in)     :: tol_in
    integer(ip), optional, intent(in)  :: max_iter_in
    integer(ip), optional, intent(out) :: iter_out
    real(dp), optional, intent(out)    :: res_out
    real(dp), optional, intent(in)     :: u_top_lid

    real(dp), allocatable :: b_hat_u(:,:), b_hat_v(:,:), b_hat_p(:,:)
    real(dp), allocatable :: r_u(:,:), r_v(:,:), r_p(:,:)
    real(dp), allocatable :: r0_u(:,:), r0_v(:,:), r0_p(:,:)
    real(dp), allocatable :: p_u(:,:), p_v(:,:), p_p(:,:)
    real(dp), allocatable :: v_u(:,:), v_v(:,:), v_p(:,:)
    real(dp), allocatable :: s_u(:,:), s_v(:,:), s_p(:,:)
    real(dp), allocatable :: t_u(:,:), t_v(:,:), t_p(:,:)
    real(dp), allocatable :: z_u(:,:), z_v(:,:), z_p(:,:)

    real(dp) :: tol, norm_b, norm_r, rel_res
    real(dp) :: rho_old, rho_new, alpha, beta, omega_val
    real(dp) :: rv_dot, tt_dot, ts_dot
    real(dp) :: u_lid
    integer(ip) :: max_iter, k

    tol = 1.0e-5_dp
    if (present(tol_in)) tol = tol_in
    max_iter = 100
    if (present(max_iter_in)) max_iter = max_iter_in

    u_lid = zero
    if (present(u_top_lid)) u_lid = u_top_lid

    allocate(b_hat_u(0-nh:nx+nh, 1-nh:ny+nh), b_hat_v(1-nh:nx+nh, 0-nh:ny+nh), b_hat_p(1-nh:nx+nh, 1-nh:ny+nh))
    allocate(r_u(0-nh:nx+nh, 1-nh:ny+nh), r_v(1-nh:nx+nh, 0-nh:ny+nh), r_p(1-nh:nx+nh, 1-nh:ny+nh))
    allocate(r0_u(0-nh:nx+nh, 1-nh:ny+nh), r0_v(1-nh:nx+nh, 0-nh:ny+nh), r0_p(1-nh:nx+nh, 1-nh:ny+nh))
    allocate(p_u(0-nh:nx+nh, 1-nh:ny+nh), p_v(1-nh:nx+nh, 0-nh:ny+nh), p_p(1-nh:nx+nh, 1-nh:ny+nh))
    allocate(v_u(0-nh:nx+nh, 1-nh:ny+nh), v_v(1-nh:nx+nh, 0-nh:ny+nh), v_p(1-nh:nx+nh, 1-nh:ny+nh))
    allocate(s_u(0-nh:nx+nh, 1-nh:ny+nh), s_v(1-nh:nx+nh, 0-nh:ny+nh), s_p(1-nh:nx+nh, 1-nh:ny+nh))
    allocate(t_u(0-nh:nx+nh, 1-nh:ny+nh), t_v(1-nh:nx+nh, 0-nh:ny+nh), t_p(1-nh:nx+nh, 1-nh:ny+nh))
    allocate(z_u(0-nh:nx+nh, 1-nh:ny+nh), z_v(1-nh:nx+nh, 0-nh:ny+nh), z_p(1-nh:nx+nh, 1-nh:ny+nh))

    call apply_velocity_bcs(X_u, X_v, nx, ny, nh)

    ! b_hat = P^(-1) * b
    call apply_pcd_precond(B_u, B_v, B_p, b_hat_u, b_hat_v, b_hat_p, &
                           u_adv, v_adv, rho_u, rho_v, &
                           kappa, zeta, eta, dt_step, nx, ny, nh)

    norm_b = sqrt(coupled_dot(b_hat_u, b_hat_v, b_hat_p, b_hat_u, b_hat_v, b_hat_p, nx, ny, nh))
    if (norm_b < 1.0e-14_dp) norm_b = one

    ! r^(0) = b_hat - P^(-1) * A * x^(0)
    call matvec_coupled(X_u, X_v, X_p, z_u, z_v, z_p, &
                        u_adv, v_adv, rho_c, rho_u, rho_v, &
                        kappa, zeta, eta, dt_step, nx, ny, nh)
    call apply_pcd_precond(z_u, z_v, z_p, v_u, v_v, v_p, &
                           u_adv, v_adv, rho_u, rho_v, &
                           kappa, zeta, eta, dt_step, nx, ny, nh)

    r_u = b_hat_u - v_u
    r_v = b_hat_v - v_v
    r_p = b_hat_p - v_p

    r0_u = r_u
    r0_v = r_v
    r0_p = r_p

    p_u = r_u
    p_v = r_v
    p_p = r_p

    rho_old = coupled_dot(r0_u, r0_v, r0_p, r_u, r_v, r_p, nx, ny, nh)
    norm_r  = sqrt(coupled_dot(r_u, r_v, r_p, r_u, r_v, r_p, nx, ny, nh))
    rel_res = norm_r / norm_b

    alpha = zero
    omega_val = one

    do k = 1, max_iter
      if (rel_res < tol) exit

      ! v = P^(-1) * A * p
      call apply_velocity_bcs(p_u, p_v, nx, ny, nh)
      call matvec_coupled(p_u, p_v, p_p, z_u, z_v, z_p, &
                          u_adv, v_adv, rho_c, rho_u, rho_v, &
                          kappa, zeta, eta, dt_step, nx, ny, nh)
      call apply_pcd_precond(z_u, z_v, z_p, v_u, v_v, v_p, &
                             u_adv, v_adv, rho_u, rho_v, &
                             kappa, zeta, eta, dt_step, nx, ny, nh)
      call apply_velocity_bcs(v_u, v_v, nx, ny, nh)

      rv_dot = coupled_dot(r0_u, r0_v, r0_p, v_u, v_v, v_p, nx, ny, nh)
      if (abs(rv_dot) < 1.0e-30_dp) exit
      alpha = rho_old / rv_dot

      s_u = r_u - alpha * v_u
      s_v = r_v - alpha * v_v
      s_p = r_p - alpha * v_p

      ! t = P^(-1) * A * s
      call apply_velocity_bcs(s_u, s_v, nx, ny, nh)
      call matvec_coupled(s_u, s_v, s_p, z_u, z_v, z_p, &
                          u_adv, v_adv, rho_c, rho_u, rho_v, &
                          kappa, zeta, eta, dt_step, nx, ny, nh)
      call apply_pcd_precond(z_u, z_v, z_p, t_u, t_v, t_p, &
                             u_adv, v_adv, rho_u, rho_v, &
                             kappa, zeta, eta, dt_step, nx, ny, nh)
      call apply_velocity_bcs(t_u, t_v, nx, ny, nh)

      tt_dot = coupled_dot(t_u, t_v, t_p, t_u, t_v, t_p, nx, ny, nh)
      ts_dot = coupled_dot(t_u, t_v, t_p, s_u, s_v, s_p, nx, ny, nh)
      if (tt_dot > 1.0e-30_dp) then
        omega_val = ts_dot / tt_dot
      else
        omega_val = zero
      end if

      X_u = X_u + alpha * p_u + omega_val * s_u
      X_v = X_v + alpha * p_v + omega_val * s_v
      X_p = X_p + alpha * p_p + omega_val * s_p

      r_u = s_u - omega_val * t_u
      r_v = s_v - omega_val * t_v
      r_p = s_p - omega_val * t_p

      norm_r  = sqrt(coupled_dot(r_u, r_v, r_p, r_u, r_v, r_p, nx, ny, nh))
      rel_res = norm_r / norm_b

      if (rel_res < tol) exit
      if (abs(omega_val) < 1.0e-30_dp) exit

      rho_new = coupled_dot(r0_u, r0_v, r0_p, r_u, r_v, r_p, nx, ny, nh)
      if (abs(rho_old) < 1.0e-30_dp) exit
      beta = (rho_new / rho_old) * (alpha / omega_val)
      rho_old = rho_new

      p_u = r_u + beta * (p_u - omega_val * v_u)
      p_v = r_v + beta * (p_v - omega_val * v_v)
      p_p = r_p + beta * (p_p - omega_val * v_p)

      call apply_velocity_bcs(X_u, X_v, nx, ny, nh)
    end do

    if (present(u_top_lid)) then
      call apply_velocity_bcs(X_u, X_v, nx, ny, nh, u_top_lid=u_lid)
    end if

    if (present(iter_out)) iter_out = min(k, max_iter)
    if (present(res_out))  res_out  = rel_res

    deallocate(b_hat_u, b_hat_v, b_hat_p)
    deallocate(r_u, r_v, r_p, r0_u, r0_v, r0_p)
    deallocate(p_u, p_v, p_p, v_u, v_v, v_p)
    deallocate(s_u, s_v, s_p, t_u, t_v, t_p)
    deallocate(z_u, z_v, z_p)

  end subroutine bicgstab_solve

  ! ----------------------------------------------------------------------------
  ! 2. Solveur BiCGStab(2) préconditionné à gauche (ordre l = 2)
  !    (Listing 4.2 de la thèse de Mohamed El Ouafa)
  ! ----------------------------------------------------------------------------
  subroutine bicgstab2_solve(X_u, X_v, X_p, B_u, B_v, B_p, &
                             u_adv, v_adv, rho_c, rho_u, rho_v, &
                             kappa, zeta, eta, dt_step, nx, ny, nh, &
                             tol_in, max_iter_in, iter_out, res_out, u_top_lid)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: X_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(inout) :: X_v
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: X_p
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: B_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)    :: B_v
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: B_p
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: u_adv
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)    :: v_adv
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: rho_c
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: rho_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)    :: rho_v
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: kappa
    real(dp), dimension(0-nh:nx+nh, 0-nh:ny+nh), intent(in)    :: zeta, eta
    real(dp), intent(in) :: dt_step
    real(dp), optional, intent(in)     :: tol_in
    integer(ip), optional, intent(in)  :: max_iter_in
    integer(ip), optional, intent(out) :: iter_out
    real(dp), optional, intent(out)    :: res_out
    real(dp), optional, intent(in)     :: u_top_lid

    ! Vecteurs de travail Krylov pour le système couplé
    real(dp), allocatable :: b_hat_u(:,:), b_hat_v(:,:), b_hat_p(:,:)
    real(dp), allocatable :: r_u(:,:), r_v(:,:), r_p(:,:)
    real(dp), allocatable :: r_tilde_u(:,:), r_tilde_v(:,:), r_tilde_p(:,:)
    real(dp), allocatable :: u_vec_u(:,:), u_vec_v(:,:), u_vec_p(:,:)
    real(dp), allocatable :: v_vec_u(:,:), v_vec_v(:,:), v_vec_p(:,:)
    real(dp), allocatable :: w_vec_u(:,:), w_vec_v(:,:), w_vec_p(:,:)
    real(dp), allocatable :: s_vec_u(:,:), s_vec_v(:,:), s_vec_p(:,:)
    real(dp), allocatable :: t_vec_u(:,:), t_vec_v(:,:), t_vec_p(:,:)
    real(dp), allocatable :: z_u(:,:), z_v(:,:), z_p(:,:)
    real(dp), allocatable :: q_u(:,:), q_v(:,:), q_p(:,:)

    real(dp) :: tol, norm_b, norm_r, rel_res
    real(dp) :: xi0, xi1, alpha, beta, gamma_val
    real(dp) :: omega1, omega2, lambda_val, nu_val, tau_val, det_val, c1_dot, c2_dot
    real(dp) :: u_lid
    integer(ip) :: max_iter, k

    tol = 1.0e-5_dp
    if (present(tol_in)) tol = tol_in
    max_iter = 100
    if (present(max_iter_in)) max_iter = max_iter_in

    u_lid = zero
    if (present(u_top_lid)) u_lid = u_top_lid

    allocate(b_hat_u(0-nh:nx+nh, 1-nh:ny+nh), b_hat_v(1-nh:nx+nh, 0-nh:ny+nh), b_hat_p(1-nh:nx+nh, 1-nh:ny+nh))
    allocate(r_u(0-nh:nx+nh, 1-nh:ny+nh), r_v(1-nh:nx+nh, 0-nh:ny+nh), r_p(1-nh:nx+nh, 1-nh:ny+nh))
    allocate(r_tilde_u(0-nh:nx+nh, 1-nh:ny+nh), r_tilde_v(1-nh:nx+nh, 0-nh:ny+nh), r_tilde_p(1-nh:nx+nh, 1-nh:ny+nh))
    allocate(u_vec_u(0-nh:nx+nh, 1-nh:ny+nh), u_vec_v(1-nh:nx+nh, 0-nh:ny+nh), u_vec_p(1-nh:nx+nh, 1-nh:ny+nh))
    allocate(v_vec_u(0-nh:nx+nh, 1-nh:ny+nh), v_vec_v(1-nh:nx+nh, 0-nh:ny+nh), v_vec_p(1-nh:nx+nh, 1-nh:ny+nh))
    allocate(w_vec_u(0-nh:nx+nh, 1-nh:ny+nh), w_vec_v(1-nh:nx+nh, 0-nh:ny+nh), w_vec_p(1-nh:nx+nh, 1-nh:ny+nh))
    allocate(s_vec_u(0-nh:nx+nh, 1-nh:ny+nh), s_vec_v(1-nh:nx+nh, 0-nh:ny+nh), s_vec_p(1-nh:nx+nh, 1-nh:ny+nh))
    allocate(t_vec_u(0-nh:nx+nh, 1-nh:ny+nh), t_vec_v(1-nh:nx+nh, 0-nh:ny+nh), t_vec_p(1-nh:nx+nh, 1-nh:ny+nh))
    allocate(z_u(0-nh:nx+nh, 1-nh:ny+nh), z_v(1-nh:nx+nh, 0-nh:ny+nh), z_p(1-nh:nx+nh, 1-nh:ny+nh))
    allocate(q_u(0-nh:nx+nh, 1-nh:ny+nh), q_v(1-nh:nx+nh, 0-nh:ny+nh), q_p(1-nh:nx+nh, 1-nh:ny+nh))

    call apply_velocity_bcs(X_u, X_v, nx, ny, nh)

    ! 1. Préconditionnement à gauche du second membre : b_hat = P^(-1) * b
    call apply_pcd_precond(B_u, B_v, B_p, b_hat_u, b_hat_v, b_hat_p, &
                           u_adv, v_adv, rho_u, rho_v, &
                           kappa, zeta, eta, dt_step, nx, ny, nh)

    norm_b = sqrt(coupled_dot(b_hat_u, b_hat_v, b_hat_p, b_hat_u, b_hat_v, b_hat_p, nx, ny, nh))
    if (norm_b < 1.0e-14_dp) norm_b = one

    ! 2. Résidu initial : r^(0) = b_hat - P^(-1) * A * x^(0)
    call matvec_coupled(X_u, X_v, X_p, z_u, z_v, z_p, &
                        u_adv, v_adv, rho_c, rho_u, rho_v, &
                        kappa, zeta, eta, dt_step, nx, ny, nh)
    call apply_pcd_precond(z_u, z_v, z_p, q_u, q_v, q_p, &
                           u_adv, v_adv, rho_u, rho_v, &
                           kappa, zeta, eta, dt_step, nx, ny, nh)

    r_u = b_hat_u - q_u
    r_v = b_hat_v - q_v
    r_p = b_hat_p - q_p

    r_tilde_u = r_u
    r_tilde_v = r_v
    r_tilde_p = r_p

    u_vec_u = zero
    u_vec_v = zero
    u_vec_p = zero

    v_vec_u = zero
    v_vec_v = zero
    v_vec_p = zero

    xi0 = one
    omega1 = zero
    omega2 = one
    alpha = zero

    norm_r = sqrt(coupled_dot(r_u, r_v, r_p, r_u, r_v, r_p, nx, ny, nh))
    rel_res = norm_r / norm_b

    ! 3. Boucle principale BiCGStab(2)
    do k = 1, max_iter
      if (rel_res < tol) exit

      if (abs(omega2) > 1.0e-14_dp) then
        xi0 = - omega2 * xi0
      else if (abs(omega1) > 1.0e-14_dp) then
        xi0 = - omega1 * xi0
      else
        xi0 = - xi0
      end if

      ! =======================================================================
      ! ÉTAPE PAIRE
      ! =======================================================================
      xi1 = coupled_dot(r_tilde_u, r_tilde_v, r_tilde_p, r_u, r_v, r_p, nx, ny, nh)
      if (abs(xi0) < 1.0e-30_dp) exit
      beta = alpha * (xi1 / xi0)
      xi0 = xi1

      u_vec_u = r_u - beta * u_vec_u
      u_vec_v = r_v - beta * u_vec_v
      u_vec_p = r_p - beta * u_vec_p
      call apply_velocity_bcs(u_vec_u, u_vec_v, nx, ny, nh)

      ! z = A * u
      call matvec_coupled(u_vec_u, u_vec_v, u_vec_p, z_u, z_v, z_p, &
                          u_adv, v_adv, rho_c, rho_u, rho_v, &
                          kappa, zeta, eta, dt_step, nx, ny, nh)
      ! v = P^(-1) * z
      call apply_pcd_precond(z_u, z_v, z_p, v_vec_u, v_vec_v, v_vec_p, &
                             u_adv, v_adv, rho_u, rho_v, &
                             kappa, zeta, eta, dt_step, nx, ny, nh)
      call apply_velocity_bcs(v_vec_u, v_vec_v, nx, ny, nh)

      gamma_val = coupled_dot(r_tilde_u, r_tilde_v, r_tilde_p, v_vec_u, v_vec_v, v_vec_p, nx, ny, nh)
      if (abs(gamma_val) < 1.0e-30_dp) exit
      alpha = xi0 / gamma_val

      r_u = r_u - alpha * v_vec_u
      r_v = r_v - alpha * v_vec_v
      r_p = r_p - alpha * v_vec_p

      ! z = A * r
      call matvec_coupled(r_u, r_v, r_p, z_u, z_v, z_p, &
                          u_adv, v_adv, rho_c, rho_u, rho_v, &
                          kappa, zeta, eta, dt_step, nx, ny, nh)
      ! s = P^(-1) * z
      call apply_pcd_precond(z_u, z_v, z_p, s_vec_u, s_vec_v, s_vec_p, &
                             u_adv, v_adv, rho_u, rho_v, &
                             kappa, zeta, eta, dt_step, nx, ny, nh)
      call apply_velocity_bcs(s_vec_u, s_vec_v, nx, ny, nh)

      X_u = X_u + alpha * u_vec_u
      X_v = X_v + alpha * u_vec_v
      X_p = X_p + alpha * u_vec_p

      ! =======================================================================
      ! ÉTAPE IMPAIRE
      ! =======================================================================
      xi1 = coupled_dot(r_tilde_u, r_tilde_v, r_tilde_p, s_vec_u, s_vec_v, s_vec_p, nx, ny, nh)
      if (abs(xi0) < 1.0e-30_dp) exit
      beta = alpha * (xi1 / xi0)
      xi0 = xi1

      v_vec_u = s_vec_u - beta * v_vec_u
      v_vec_v = s_vec_v - beta * v_vec_v
      v_vec_p = s_vec_p - beta * v_vec_p
      call apply_velocity_bcs(v_vec_u, v_vec_v, nx, ny, nh)

      ! z = A * v
      call matvec_coupled(v_vec_u, v_vec_v, v_vec_p, z_u, z_v, z_p, &
                          u_adv, v_adv, rho_c, rho_u, rho_v, &
                          kappa, zeta, eta, dt_step, nx, ny, nh)
      ! w = P^(-1) * z
      call apply_pcd_precond(z_u, z_v, z_p, w_vec_u, w_vec_v, w_vec_p, &
                             u_adv, v_adv, rho_u, rho_v, &
                             kappa, zeta, eta, dt_step, nx, ny, nh)
      call apply_velocity_bcs(w_vec_u, w_vec_v, nx, ny, nh)

      gamma_val = coupled_dot(r_tilde_u, r_tilde_v, r_tilde_p, w_vec_u, w_vec_v, w_vec_p, nx, ny, nh)
      if (abs(gamma_val) < 1.0e-30_dp) exit
      alpha = xi0 / gamma_val

      u_vec_u = r_u - beta * u_vec_u
      u_vec_v = r_v - beta * u_vec_v
      u_vec_p = r_p - beta * u_vec_p

      r_u = r_u - alpha * v_vec_u
      r_v = r_v - alpha * v_vec_v
      r_p = r_p - alpha * v_vec_p

      s_vec_u = s_vec_u - alpha * w_vec_u
      s_vec_v = s_vec_v - alpha * w_vec_v
      s_vec_p = s_vec_p - alpha * w_vec_p

      ! z = A * s
      call matvec_coupled(s_vec_u, s_vec_v, s_vec_p, z_u, z_v, z_p, &
                          u_adv, v_adv, rho_c, rho_u, rho_v, &
                          kappa, zeta, eta, dt_step, nx, ny, nh)
      ! t = P^(-1) * z
      call apply_pcd_precond(z_u, z_v, z_p, t_vec_u, t_vec_v, t_vec_p, &
                             u_adv, v_adv, rho_u, rho_v, &
                             kappa, zeta, eta, dt_step, nx, ny, nh)
      call apply_velocity_bcs(t_vec_u, t_vec_v, nx, ny, nh)

      ! =======================================================================
      ! PARTIE GCR(2) : Minimisation de résidu d'ordre 2
      ! =======================================================================
      omega1     = coupled_dot(r_u, r_v, r_p, s_vec_u, s_vec_v, s_vec_p, nx, ny, nh)
      lambda_val = coupled_dot(s_vec_u, s_vec_v, s_vec_p, s_vec_u, s_vec_v, s_vec_p, nx, ny, nh)
      nu_val     = coupled_dot(s_vec_u, s_vec_v, s_vec_p, t_vec_u, t_vec_v, t_vec_p, nx, ny, nh)
      tau_val    = coupled_dot(t_vec_u, t_vec_v, t_vec_p, t_vec_u, t_vec_v, t_vec_p, nx, ny, nh)
      omega2     = coupled_dot(r_u, r_v, r_p, t_vec_u, t_vec_v, t_vec_p, nx, ny, nh)

      if (lambda_val * tau_val > 1.0e-30_dp) then
        det_val = lambda_val * tau_val - nu_val**2
        if (det_val > 1.0e-14_dp * (lambda_val * tau_val)) then
          c1_dot = omega1
          c2_dot = omega2
          omega1 = (c1_dot * tau_val - c2_dot * nu_val) / det_val
          omega2 = (c2_dot * lambda_val - c1_dot * nu_val) / det_val
        else if (lambda_val > 1.0e-30_dp) then
          omega1 = omega1 / lambda_val
          omega2 = zero
        else
          omega1 = zero
          omega2 = zero
        end if
      else if (lambda_val > 1.0e-30_dp) then
        omega1 = omega1 / lambda_val
        omega2 = zero
      else
        omega1 = zero
        omega2 = zero
      end if

      ! Mise à jour finale de l'estimation de solution
      X_u = X_u + omega1 * r_u + omega2 * s_vec_u + alpha * u_vec_u
      X_v = X_v + omega1 * r_v + omega2 * s_vec_v + alpha * u_vec_v
      X_p = X_p + omega1 * r_p + omega2 * s_vec_p + alpha * u_vec_p

      ! Mise à jour finale du résidu
      r_u = r_u - omega1 * s_vec_u - omega2 * t_vec_u
      r_v = r_v - omega1 * s_vec_v - omega2 * t_vec_v
      r_p = r_p - omega1 * s_vec_p - omega2 * t_vec_p

      u_vec_u = u_vec_u - omega1 * v_vec_u - omega2 * w_vec_u
      u_vec_v = u_vec_v - omega1 * v_vec_v - omega2 * w_vec_v
      u_vec_p = u_vec_p - omega1 * v_vec_p - omega2 * w_vec_p

      norm_r = sqrt(coupled_dot(r_u, r_v, r_p, r_u, r_v, r_p, nx, ny, nh))
      rel_res = norm_r / norm_b

      call apply_velocity_bcs(X_u, X_v, nx, ny, nh)
    end do

    if (present(u_top_lid)) then
      call apply_velocity_bcs(X_u, X_v, nx, ny, nh, u_top_lid=u_lid)
    end if

    if (present(iter_out)) iter_out = min(k, max_iter)
    if (present(res_out))  res_out  = rel_res

    deallocate(b_hat_u, b_hat_v, b_hat_p)
    deallocate(r_u, r_v, r_p)
    deallocate(r_tilde_u, r_tilde_v, r_tilde_p)
    deallocate(u_vec_u, u_vec_v, u_vec_p)
    deallocate(v_vec_u, v_vec_v, v_vec_p)
    deallocate(w_vec_u, w_vec_v, w_vec_p)
    deallocate(s_vec_u, s_vec_v, s_vec_p)
    deallocate(t_vec_u, t_vec_v, t_vec_p)
    deallocate(z_u, z_v, z_p)
    deallocate(q_u, q_v, q_p)

  end subroutine bicgstab2_solve

end module mod_bicgstab2
