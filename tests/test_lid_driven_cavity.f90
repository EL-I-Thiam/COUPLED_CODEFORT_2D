! ==============================================================================
! PROGRAMME DE TEST : test_lid_driven_cavity.f90
! Validation du solveur tout-couplé BiCGStab(2) + Préconditionneur PCD
! sur le benchmark standard de la cavité entraînée 2D (Lid-Driven Cavity)
! Projet Solveur 2D Tout-Couplé (FUGU-2D) - Fortran 2008 / OpenMPI
! ==============================================================================
program test_lid_driven_cavity
  use mod_precision,  only: dp, ip, zero, one, half, two
  use mod_params,     only: Nx_global, Ny_global, Lx, Ly, rho1, rho2, mu1, mu2
  use mod_mpi,        only: my_rank, n_procs, dims, coords, nx_loc, ny_loc, &
                            i_start, j_start, dx_grid, dy_grid, &
                            comm_cart, mpi_init_env, mpi_finalize_env, &
                            mpi_setup_cartesian_mesh, mpi_exchange_halos_scalar, &
                            mpi_exchange_halos_u, mpi_exchange_halos_v
  use mod_viscosity,  only: compute_mixture_densities, compute_mixture_viscosities
  use mod_operators,  only: op_divergence, apply_velocity_bcs
  use mod_bicgstab2,  only: bicgstab_solve, bicgstab2_solve
  use mpi
  implicit none

  integer(ip), parameter :: nh = 2
  real(dp), allocatable :: u_vel(:,:), v_vel(:,:), p_field(:,:)
  real(dp), allocatable :: u_old(:,:), v_old(:,:)
  real(dp), allocatable :: c_field(:,:)
  real(dp), allocatable :: rho_c(:,:), rho_u(:,:), rho_v(:,:)
  real(dp), allocatable :: kappa(:,:), zeta(:,:), eta(:,:)
  real(dp), allocatable :: B_u(:,:), B_v(:,:), B_p(:,:)
  real(dp), allocatable :: div_field(:,:)

  real(dp) :: dt_val, u_lid, reynolds
  real(dp) :: diff_u_loc, diff_u_glob
  real(dp) :: div_max_loc, div_max_glob
  real(dp) :: res_bicg
  integer(ip) :: iter_bicg, step, max_steps
  integer(ip) :: i_loc, j_loc, i_glob, j_glob, ierr

  ! Profils de vitesse pour comparaison avec Botella & Peyret (1998)
  real(dp) :: u_min_loc, u_min_glob, y_at_umin
  real(dp) :: v_min_loc, v_min_glob, v_max_loc, v_max_glob
  real(dp) :: x_at_vmin, x_at_vmax
  real(dp) :: x_pos, y_pos

  ! 1. Initialisation de l'environnement MPI
  call mpi_init_env()

  ! Domaine [0, 1] x [0, 1], résolution 32 x 32
  Nx_global = 32
  Ny_global = 32
  Lx = one
  Ly = one

  ! Domaine fermé non-périodique (parois étanches avec paroi supérieure entraînée)
  call mpi_setup_cartesian_mesh(Nx_global, Ny_global, Lx, Ly, periodic_x=.false., periodic_y=.false.)

  u_lid    = one
  reynolds = 100.0_dp
  rho1     = one
  rho2     = one
  mu1      = (rho1 * u_lid * Lx) / reynolds  ! mu = 0.01 pour Re = 100
  mu2      = mu1

  dt_val    = 0.01_dp
  max_steps = 350

  if (my_rank == 0) then
    write(*,'(A)') "================================================================"
    write(*,'(A)') " TEST DE VALIDATION SEMAINE 3 : Cavite Entrainee 2D (Re = 100)  "
    write(*,'(A)') " Solveur Tout-Couple BiCGStab(2) + Preconditionneur Bloc PCD   "
    write(*,'(A)') "================================================================"
    write(*,'(A,I4,A,I3,A,I3)') " Processus MPI : ", n_procs, " | Grille cartésienne : ", dims(1), " x ", dims(2)
    write(*,'(A,I5,A,I5,A,ES10.3)') " Résolution : ", Nx_global, " x ", Ny_global, " | dt = ", dt_val
    write(*,'(A,F8.2,A,F8.4)') " Nombre de Reynolds Re = ", reynolds, " | Viscosite mu = ", mu1
    write(*,'(A)') "----------------------------------------------------------------"
  end if

  ! 2. Allocations dynamiques des tableaux locaux
  allocate(u_vel(0-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(v_vel(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(p_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(u_old(0-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(v_old(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(c_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))

  allocate(rho_c(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(rho_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(rho_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))

  allocate(kappa(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(zeta(0-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(eta(0-nh:nx_loc+nh, 0-nh:ny_loc+nh))

  allocate(B_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(B_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(B_p(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(div_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))

  ! 3. Initialisation de l'écoulement au repos
  u_vel   = zero
  v_vel   = zero
  p_field = zero
  c_field = one  ! Écoulement monophasique (fluide 1)

  call compute_mixture_densities(c_field, rho_c, rho_u, rho_v, nx_loc, ny_loc, nh)
  call compute_mixture_viscosities(c_field, kappa, zeta, eta, nx_loc, ny_loc, nh)

  call apply_velocity_bcs(u_vel, v_vel, nx_loc, ny_loc, nh, u_top_lid=u_lid)

  ! 4. Boucle temporelle vers l'état stationnaire
  do step = 1, max_steps
    u_old = u_vel
    v_old = v_vel

    ! Construction du second membre monolithique :
    ! B_u = (rho_u / dt) * u^n
    ! B_v = (rho_v / dt) * v^n
    ! B_p = 0
    B_u = zero
    do j_loc = 1, ny_loc
      do i_loc = 1, nx_loc-1
        B_u(i_loc, j_loc) = (rho_u(i_loc, j_loc) / dt_val) * u_old(i_loc, j_loc)
      end do
      B_u(0, j_loc) = zero
      B_u(nx_loc, j_loc) = zero
    end do

    ! Terme inhomogène de paroi supérieure entraînée (y = Ly)
    ! Provient de la discrétisation de -div(T)_u à la paroi mobile u_lid
    if (coords(2) == dims(2) - 1) then
      do i_loc = 1, nx_loc-1
        B_u(i_loc, ny_loc) = B_u(i_loc, ny_loc) + &
          (two * (zeta(i_loc, ny_loc) - eta(i_loc, ny_loc)) / (dy_grid * dy_grid)) * u_lid
      end do
    end if

    B_v = zero
    do j_loc = 1, ny_loc-1
      do i_loc = 1, nx_loc
        B_v(i_loc, j_loc) = (rho_v(i_loc, j_loc) / dt_val) * v_old(i_loc, j_loc)
      end do
    end do
    do i_loc = 1, nx_loc
      B_v(i_loc, 0) = zero
      B_v(i_loc, ny_loc) = zero
    end do

    B_p = zero

    ! Résolution du système selle couplé par BiCGStab(2) + PCD
    call bicgstab_solve(u_vel, v_vel, p_field, B_u, B_v, B_p, &
                         u_old, v_old, rho_c, rho_u, rho_v, &
                         kappa, zeta, eta, dt_val, nx_loc, ny_loc, nh, &
                         tol_in=1.0e-4_dp, max_iter_in=100, &
                         iter_out=iter_bicg, res_out=res_bicg, &
                         u_top_lid=u_lid)

    ! Calcul de la variation ||u^(n+1) - u^n||_inf pour convergence stationnaire
    diff_u_loc = zero
    do j_loc = 1, ny_loc
      do i_loc = 1, nx_loc
        diff_u_loc = max(diff_u_loc, abs(u_vel(i_loc, j_loc) - u_old(i_loc, j_loc)))
      end do
    end do
    call MPI_Allreduce(diff_u_loc, diff_u_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_cart, ierr)

    if (my_rank == 0 .and. (step == 1 .or. step <= 5 .or. mod(step, 10) == 0 .or. step == max_steps)) then
      write(*,'(A,I4,A,ES11.4,A,I3,A,ES10.3)') &
        " Pas ", step, " | ||du/dt||_inf = ", diff_u_glob / dt_val, &
        " | BiCGStab(2) iter = ", iter_bicg, " | Res = ", res_bicg
      write(*,*) "Max u:", maxval(abs(u_vel(1:nx_loc-1, 1:ny_loc))), &
                 "Max v:", maxval(abs(v_vel(1:nx_loc, 1:ny_loc-1))), &
                 "Max p:", maxval(abs(p_field(1:nx_loc, 1:ny_loc)))
    end if

    if (diff_u_glob / dt_val < 1.0e-4_dp) exit
  end do

  ! 5. Vérification de la divergence finale (incompressibilité)
  call op_divergence(u_vel, v_vel, div_field, nx_loc, ny_loc, nh)
  div_max_loc = zero
  do j_loc = 1, ny_loc
    do i_loc = 1, nx_loc
      div_max_loc = max(div_max_loc, abs(div_field(i_loc, j_loc)))
    end do
  end do
  call MPI_Allreduce(div_max_loc, div_max_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_cart, ierr)

  ! 6. Extraction des grandeurs caractéristiques de l'écoulement
  ! Profil u(y) le long de la ligne médiane x = 0.5
  u_min_loc  = zero
  y_at_umin  = zero
  do j_loc = 1, ny_loc
    j_glob = j_start + j_loc - 1
    y_pos  = (real(j_glob, dp) - half) * dy_grid
    ! Maille la plus proche de x = 0.5
    i_glob = Nx_global / 2
    if (i_glob >= i_start .and. i_glob <= i_start + nx_loc - 1) then
      i_loc = i_glob - i_start + 1
      if (u_vel(i_loc, j_loc) < u_min_loc) then
        u_min_loc = u_vel(i_loc, j_loc)
        y_at_umin = y_pos
      end if
    end if
  end do

  ! Profil v(x) le long de la ligne médiane y = 0.5
  v_min_loc = zero
  v_max_loc = zero
  x_at_vmin = zero
  x_at_vmax = zero
  do i_loc = 1, nx_loc
    i_glob = i_start + i_loc - 1
    x_pos  = (real(i_glob, dp) - half) * dx_grid
    j_glob = Ny_global / 2
    if (j_glob >= j_start .and. j_glob <= j_start + ny_loc - 1) then
      j_loc = j_glob - j_start + 1
      if (v_vel(i_loc, j_loc) < v_min_loc) then
        v_min_loc = v_vel(i_loc, j_loc)
        x_at_vmin = x_pos
      end if
      if (v_vel(i_loc, j_loc) > v_max_loc) then
        v_max_loc = v_vel(i_loc, j_loc)
        x_at_vmax = x_pos
      end if
    end if
  end do

  call MPI_Allreduce(u_min_loc, u_min_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm_cart, ierr)
  call MPI_Allreduce(v_min_loc, v_min_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm_cart, ierr)
  call MPI_Allreduce(v_max_loc, v_max_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_cart, ierr)

  if (my_rank == 0) then
    write(*,'(A)') "----------------------------------------------------------------"
    write(*,'(A)') " RESULTATS DU BENCHMARK COMPARATIF (Re = 100) :                 "
    write(*,'(A)') "----------------------------------------------------------------"
    write(*,'(A,ES12.4)') " Residu max divergence finale : ||div u||_inf = ", div_max_glob
    write(*,'(A)') " Profils centraux vs Spectral Botella & Peyret (1998) :"
    write(*,'(A,F10.5,A)') "  - u_min calcule   = ", u_min_glob, " (Reference Botella : -0.21404)"
    write(*,'(A,F10.5,A)') "  - v_min calcule   = ", v_min_glob, " (Reference Botella : -0.25380)"
    write(*,'(A,F10.5,A)') "  - v_max calcule   = ", v_max_glob, " (Reference Botella :  0.17957)"
    write(*,'(A)') "----------------------------------------------------------------"
    if (abs(u_min_glob) > 0.15_dp .and. abs(v_min_glob) > 0.15_dp .and. div_max_glob < 1.0e-2_dp) then
      write(*,'(A)') " --> [SUCCES] Cavite entrainee Re = 100 validee avec le solveur couple !"
    else
      write(*,'(A)') " --> [ATTENTION] Verifier le profil de vitesse."
    end if
    write(*,'(A)') "================================================================"
  end if

  deallocate(u_vel, v_vel, p_field, u_old, v_old, c_field)
  deallocate(rho_c, rho_u, rho_v, kappa, zeta, eta)
  deallocate(B_u, B_v, B_p, div_field)

  call mpi_finalize_env()

end program test_lid_driven_cavity
