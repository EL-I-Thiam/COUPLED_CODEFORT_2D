! ==============================================================================
! MODULE : mod_vof
! Solveur 2D Tout-Couplé (FUGU-2D) - Fortran 2008 / OpenMPI
!
! 1. Advection conservative de la fraction volumique C (Van Leer TVD)
! 2. Régularisation SVOF par équation de Helmholtz : C^s - delta_s^2 * Lapl(C^s) = C
! ==============================================================================
module mod_vof
  use mod_precision, only: dp, ip, zero, one, half, two
  use mod_mpi,       only: dx_grid, dy_grid, coords, dims, periods, comm_cart, &
                           mpi_exchange_halos_scalar, mpi_exchange_halos_u, &
                           mpi_exchange_halos_v
  use mpi
  implicit none
  private

  public :: vof_advection_tvd
  public :: svof_helmholtz_filter

contains

  ! ----------------------------------------------------------------------------
  ! Limiteur de flux de Van Leer : psi(r) = (r + |r|) / (1 + |r|)
  ! ----------------------------------------------------------------------------
  pure function van_leer_limiter(r) result(psi)
    real(dp), intent(in) :: r
    real(dp) :: psi
    if (r <= zero) then
      psi = zero
    else
      psi = (two * r) / (one + r)
    end if
  end function van_leer_limiter

  ! ----------------------------------------------------------------------------
  ! 1. Advection conservative 2D de la fraction volumique C avec limiteur TVD
  ! ----------------------------------------------------------------------------
  subroutine vof_advection_tvd(C_in, C_out, u_vel, v_vel, dt_step, nx, ny, nh, &
                               flux_x_out, flux_y_out)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)  :: C_in
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(out) :: C_out
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)  :: u_vel
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)  :: v_vel
    real(dp), intent(in) :: dt_step
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), optional, intent(out) :: flux_x_out
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), optional, intent(out) :: flux_y_out

    real(dp), dimension(0:nx, 1:ny) :: Flux_x
    real(dp), dimension(1:nx, 0:ny) :: Flux_y
    real(dp) :: u_f, v_f, cfl_loc, r_ratio, num_diff, denom_diff
    real(dp) :: C_face, inv_dx, inv_dy
    integer(ip) :: i, j

    inv_dx = one / dx_grid
    inv_dy = one / dy_grid

    ! Flux horizontaux sur faces u : i = 0 .. nx
    do j = 1, ny
      do i = 0, nx
        u_f = u_vel(i, j)
        if (u_f >= zero) then
          num_diff   = C_in(i, j) - C_in(i-1, j)
          denom_diff = C_in(i+1, j) - C_in(i, j)
          if (abs(denom_diff) > 1.0e-12_dp) then
            r_ratio = num_diff / denom_diff
          else
            r_ratio = zero
          end if
          cfl_loc = abs(u_f) * dt_step * inv_dx
          C_face = C_in(i, j) + half * (one - cfl_loc) * van_leer_limiter(r_ratio) * denom_diff
        else
          num_diff   = C_in(i+2, j) - C_in(i+1, j)
          denom_diff = C_in(i+1, j) - C_in(i, j)
          if (abs(denom_diff) > 1.0e-12_dp) then
            r_ratio = num_diff / denom_diff
          else
            r_ratio = zero
          end if
          cfl_loc = abs(u_f) * dt_step * inv_dx
          C_face = C_in(i+1, j) - half * (one - cfl_loc) * van_leer_limiter(r_ratio) * denom_diff
        end if
        Flux_x(i, j) = u_f * max(zero, min(one, C_face))
      end do
    end do

    ! Flux verticaux sur faces v : j = 0 .. ny
    do j = 0, ny
      do i = 1, nx
        v_f = v_vel(i, j)
        if (v_f >= zero) then
          num_diff   = C_in(i, j) - C_in(i, j-1)
          denom_diff = C_in(i, j+1) - C_in(i, j)
          if (abs(denom_diff) > 1.0e-12_dp) then
            r_ratio = num_diff / denom_diff
          else
            r_ratio = zero
          end if
          cfl_loc = abs(v_f) * dt_step * inv_dy
          C_face = C_in(i, j) + half * (one - cfl_loc) * van_leer_limiter(r_ratio) * denom_diff
        else
          num_diff   = C_in(i, j+2) - C_in(i, j+1)
          denom_diff = C_in(i, j+1) - C_in(i, j)
          if (abs(denom_diff) > 1.0e-12_dp) then
            r_ratio = num_diff / denom_diff
          else
            r_ratio = zero
          end if
          cfl_loc = abs(v_f) * dt_step * inv_dy
          C_face = C_in(i, j+1) - half * (one - cfl_loc) * van_leer_limiter(r_ratio) * denom_diff
        end if
        Flux_y(i, j) = v_f * max(zero, min(one, C_face))
      end do
    end do

    ! Bilan conservatif aux centres
    C_out = C_in
    do j = 1, ny
      do i = 1, nx
        C_out(i, j) = C_in(i, j) - dt_step * ( &
          (Flux_x(i, j) - Flux_x(i-1, j)) * inv_dx + &
          (Flux_y(i, j) - Flux_y(i, j-1)) * inv_dy )
        ! Bornage strict physique C in [0, 1]
        C_out(i, j) = max(zero, min(one, C_out(i, j)))
      end do
    end do

    if (present(flux_x_out)) then
      flux_x_out = zero
      flux_x_out(0:nx, 1:ny) = Flux_x(0:nx, 1:ny)
      call mpi_exchange_halos_u(flux_x_out, nx, ny, nh)
    end if
    if (present(flux_y_out)) then
      flux_y_out = zero
      flux_y_out(1:nx, 0:ny) = Flux_y(1:nx, 0:ny)
      call mpi_exchange_halos_v(flux_y_out, nx, ny, nh)
    end if

    call mpi_exchange_halos_scalar(C_out, nx, ny, nh)
  end subroutine vof_advection_tvd

  ! ----------------------------------------------------------------------------
  ! 2. Filtre Helmholtz régularisant SVOF : C^s - delta_s^2 * Lapl(C^s) = C
  ! ----------------------------------------------------------------------------
  subroutine svof_helmholtz_filter(C_raw, C_smooth, delta_filter, nx, ny, nh, max_iter, tol)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: C_raw
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: C_smooth
    real(dp), intent(in) :: delta_filter
    integer(ip), optional, intent(in) :: max_iter
    real(dp), optional, intent(in)    :: tol

    integer(ip) :: it, it_max, i, j, ierr
    real(dp) :: eps_tol, coeff_x, coeff_y, diag_c, res_loc, res_glob
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh) :: C_old

    it_max  = 100
    eps_tol = 1.0e-6_dp
    if (present(max_iter)) it_max  = max_iter
    if (present(tol))      eps_tol = tol

    coeff_x = (delta_filter / dx_grid)**2
    coeff_y = (delta_filter / dy_grid)**2
    diag_c  = one + two * (coeff_x + coeff_y)

    C_smooth = C_raw

    do it = 1, it_max
      call mpi_exchange_halos_scalar(C_smooth, nx, ny, nh)
      call apply_neumann_bc(C_smooth, nx, ny, nh)
      C_old = C_smooth

      do j = 1, ny
        do i = 1, nx
          C_smooth(i, j) = ( C_raw(i, j) + &
            coeff_x * (C_old(i+1, j) + C_old(i-1, j)) + &
            coeff_y * (C_old(i, j+1) + C_old(i, j-1)) ) / diag_c
        end do
      end do

      ! Évaluation de la convergence
      res_loc = zero
      do j = 1, ny
        do i = 1, nx
          res_loc = max(res_loc, abs(C_smooth(i, j) - C_old(i, j)))
        end do
      end do
      call MPI_Allreduce(res_loc, res_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_cart, ierr)
      if (res_glob < eps_tol) exit
    end do

    call mpi_exchange_halos_scalar(C_smooth, nx, ny, nh)
    call apply_neumann_bc(C_smooth, nx, ny, nh)
  end subroutine svof_helmholtz_filter

  ! ----------------------------------------------------------------------------
  ! Routine utilitaire interne : Neumann homogène (dC/dn = 0) aux parois
  ! ----------------------------------------------------------------------------
  subroutine apply_neumann_bc(phi, nx, ny, nh)
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
  end subroutine apply_neumann_bc

end module mod_vof
