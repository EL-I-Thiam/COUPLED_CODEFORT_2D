! ==============================================================================
! PROGRAMME DE TEST : test_static_droplet.f90
! Validation du couplage VOF, lissage SVOF et force capillaire CSF
! sur le cas standard de la goutte statique 2D (Loi de Laplace)
! Projet Solveur 2D Tout-Couplé (FUGU-2D) - Fortran 2008 / OpenMPI
! ==============================================================================
program test_static_droplet
  use mod_precision,       only: dp, ip, zero, one, half
  use mod_params,          only: Nx_global, Ny_global, Lx, Ly, sigma
  use mod_mpi,             only: my_rank, n_procs, dims, nx_loc, ny_loc, &
                                 i_start, j_start, dx_grid, dy_grid, comm_cart, &
                                 mpi_init_env, mpi_finalize_env, mpi_setup_cartesian_mesh, &
                                 mpi_exchange_halos_scalar, mpi_exchange_halos_u, &
                                 mpi_exchange_halos_v
  use mod_vof,             only: svof_helmholtz_filter
  use mod_surface_tension, only: compute_interface_curvature, compute_csf_force
  use mpi
  implicit none

  integer(ip), parameter :: nh = 2
  real(dp), allocatable :: C_raw(:,:), C_smooth(:,:), kappa_curv(:,:)
  real(dp), allocatable :: f_sigma_u(:,:), f_sigma_v(:,:)

  real(dp) :: x_c, y_c, r_dist, R_drop, delta_filter
  real(dp) :: kappa_exact, kappa_sum_loc, kappa_sum_glob
  real(dp) :: delta_p_exact, kappa_mean, err_kappa
  real(dp) :: f_max_u_loc, f_max_v_loc, f_max_glob
  integer(ip) :: i_loc, j_loc, i_glob, j_glob, ierr
  integer(ip) :: n_int_loc, n_int_glob

  call mpi_init_env()

  Nx_global = 64
  Ny_global = 64
  Lx = one
  Ly = one
  sigma = 1.0_dp
  R_drop = 0.25_dp

  call mpi_setup_cartesian_mesh(Nx_global, Ny_global, Lx, Ly, periodic_x=.false., periodic_y=.false.)

  if (my_rank == 0) then
    write(*,'(A)') "================================================================"
    write(*,'(A)') " TEST DE VALIDATION SEMAINE 4 : Goutte Statique & Tension CSF   "
    write(*,'(A)') "================================================================"
    write(*,'(A,I4,A,I3,A,I3)') " Processus MPI : ", n_procs, " | Grille cartésienne : ", dims(1), " x ", dims(2)
    write(*,'(A,I5,A,I5,A,ES10.3)') " Résolution : ", Nx_global, " x ", Ny_global, " | dx = ", dx_grid
    write(*,'(A,F6.3,A,F6.3)') " Rayon R0 = ", R_drop, " m | Tension de surface sigma = ", sigma
  end if

  allocate(C_raw(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(C_smooth(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(kappa_curv(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(f_sigma_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(f_sigma_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))

  ! 1. Initialisation géométrique de la goutte au centre (0.5, 0.5)
  do j_loc = 1-nh, ny_loc+nh
    j_glob = j_start + j_loc - 1
    y_c = (real(j_glob, dp) - half) * dy_grid
    do i_loc = 1-nh, nx_loc+nh
      i_glob = i_start + i_loc - 1
      x_c = (real(i_glob, dp) - half) * dx_grid
      r_dist = sqrt((x_c - 0.5_dp)**2 + (y_c - 0.5_dp)**2)
      if (r_dist <= R_drop) then
        C_raw(i_loc, j_loc) = one
      else
        C_raw(i_loc, j_loc) = zero
      end if
    end do
  end do

  call mpi_exchange_halos_scalar(C_raw, nx_loc, ny_loc, nh)

  ! 2. Filtrage SVOF régularisant par équation de Helmholtz (delta = 1.5 * dx)
  delta_filter = 1.5_dp * dx_grid
  call svof_helmholtz_filter(C_raw, C_smooth, delta_filter, nx_loc, ny_loc, nh, max_iter=150, tol=1.0e-7_dp)

  ! 3. Calcul de la courbure kappa_Gamma = - div(n)
  call compute_interface_curvature(C_smooth, kappa_curv, nx_loc, ny_loc, nh)

  ! 4. Calcul de la force volumique capillaire CSF sur les faces décalées
  call compute_csf_force(C_raw, kappa_curv, f_sigma_u, f_sigma_v, nx_loc, ny_loc, nh)

  ! 5. Évaluation statistique de la courbure à l'interface (zone de transition 0.05 < C < 0.95)
  kappa_exact = one / R_drop  ! 4.0 m^-1 en 2D
  delta_p_exact = sigma / R_drop  ! 4.0 Pa

  kappa_sum_loc = zero
  n_int_loc = 0
  do j_loc = 1, ny_loc
    j_glob = j_start + j_loc - 1
    y_c = (real(j_glob, dp) - half) * dy_grid
    do i_loc = 1, nx_loc
      i_glob = i_start + i_loc - 1
      x_c = (real(i_glob, dp) - half) * dx_grid
      r_dist = sqrt((x_c - 0.5_dp)**2 + (y_c - 0.5_dp)**2)
      ! Sélection des cellules traversées par l'interface (distance à R0 < 1.5 dx)
      if (abs(r_dist - R_drop) <= 1.5_dp * dx_grid) then
        kappa_sum_loc = kappa_sum_loc + kappa_curv(i_loc, j_loc)
        n_int_loc = n_int_loc + 1
      end if
    end do
  end do

  call MPI_Allreduce(kappa_sum_loc, kappa_sum_glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)
  call MPI_Allreduce(n_int_loc, n_int_glob, 1, MPI_INTEGER, MPI_SUM, comm_cart, ierr)

  if (n_int_glob > 0) then
    kappa_mean = kappa_sum_glob / real(n_int_glob, dp)
  else
    kappa_mean = zero
  end if
  err_kappa = abs(kappa_mean - kappa_exact) / kappa_exact * 100.0_dp

  ! 6. Évaluation de l'intensité de la force CSF
  f_max_u_loc = zero
  f_max_v_loc = zero
  do j_loc = 1, ny_loc
    do i_loc = 1, nx_loc
      f_max_u_loc = max(f_max_u_loc, abs(f_sigma_u(i_loc, j_loc)))
      f_max_v_loc = max(f_max_v_loc, abs(f_sigma_v(i_loc, j_loc)))
    end do
  end do
  call MPI_Allreduce(max(f_max_u_loc, f_max_v_loc), f_max_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_cart, ierr)

  if (my_rank == 0) then
    write(*,'(A)') "----------------------------------------------------------------"
    write(*,'(A)') " RÉSULTATS DU BENCHMARK CSF (Loi de Laplace 2D) :              "
    write(*,'(A)') "----------------------------------------------------------------"
    write(*,'(A,F10.4,A)') " Courbure théorique kappa_exact = ", kappa_exact, " m^-1"
    write(*,'(A,F10.4,A)') " Courbure calculée (SVOF) kappa = ", kappa_mean,  " m^-1"
    write(*,'(A,F8.3,A)')  " Erreur relative sur la courbure = ", err_kappa, " %"
    write(*,'(A,ES11.4)')  " Intensité max de la force CSF  = ", f_max_glob
    write(*,'(A,F10.4,A)') " Saut de pression de Laplace     = ", delta_p_exact, " Pa"
    write(*,'(A)') "----------------------------------------------------------------"
    if (err_kappa < 5.0_dp .and. f_max_glob > 10.0_dp) then
      write(*,'(A)') " --> [SUCCÈS] Courbure SVOF et force CSF conformes à la loi de Laplace !"
    else
      write(*,'(A)') " --> [ATTENTION] Écart sur la courbure ou la force capillaire."
    end if
    write(*,'(A)') "================================================================"
  end if

  deallocate(C_raw, C_smooth, kappa_curv, f_sigma_u, f_sigma_v)
  call mpi_finalize_env()

end program test_static_droplet
