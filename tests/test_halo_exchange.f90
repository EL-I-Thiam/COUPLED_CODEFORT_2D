program test_halo_exchange
  use mod_precision, only: dp, ip, zero, one, pi
  use mod_params,    only: Nx_global, Ny_global, Lx, Ly
  use mod_mpi,       only: my_rank, n_procs, dims, nx_loc, ny_loc, &
                           i_start, j_start, dx_grid, dy_grid, &
                           comm_cart, mpi_init_env, mpi_finalize_env, &
                           mpi_setup_cartesian_mesh, mpi_exchange_halos_scalar
  use mpi
  implicit none

  integer(ip), parameter :: nh = 2
  real(dp), allocatable :: phi(:,:)
  real(dp) :: x_glob, y_glob, phi_exact, err_local, err_global
  integer(ip) :: i_loc, j_loc, i_glob, j_glob, ierr

  ! 1. Initialisation MPI et maillage
  call mpi_init_env()

  Nx_global = 64
  Ny_global = 64
  Lx = one
  Ly = one

  ! Configuration périodique en X et Y pour tester toutes les frontières et raccords
  call mpi_setup_cartesian_mesh(Nx_global, Ny_global, Lx, Ly, periodic_x=.true., periodic_y=.true.)

  if (my_rank == 0) then
    write(*,'(A)') "================================================================"
    write(*,'(A)') " TEST DE VALIDATION SEMAINE 1 : Topologie MPI & Echange de halos "
    write(*,'(A)') "================================================================"
    write(*,'(A,I4,A,I3,A,I3)') " Nombre de processus : ", n_procs, &
                               " | Grille de calcul MPI : ", dims(1), " x ", dims(2)
    write(*,'(A,I5,A,I5)') " Taille globale : ", Nx_global, " x ", Ny_global
  end if

  ! 2. Allocation du champ local avec halos
  allocate(phi(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  phi = zero

  ! 3. Initialisation de l'intérieur du sous-domaine avec une fonction analytique lisse
  do j_loc = 1, ny_loc
    j_glob = j_start + j_loc - 1
    y_glob = (real(j_glob, dp) - 0.5_dp) * dy_grid
    do i_loc = 1, nx_loc
      i_glob = i_start + i_loc - 1
      x_glob = (real(i_glob, dp) - 0.5_dp) * dx_grid
      phi(i_loc, j_loc) = sin(2.0_dp * pi * x_glob / Lx) * cos(2.0_dp * pi * y_glob / Ly)
    end do
  end do

  ! 4. Échange des couches de cellules fantômes
  call mpi_exchange_halos_scalar(phi, nx_loc, ny_loc, nh)

  ! 5. Vérification de l'exactitude des valeurs reçues dans les halos
  err_local = zero

  ! Vérification des halos Est et Ouest (y compris les coins)
  do j_loc = 1-nh, ny_loc+nh
    j_glob = j_start + j_loc - 1
    ! Gestion de la périodicité globale en Y
    if (j_glob < 1) j_glob = j_glob + Ny_global
    if (j_glob > Ny_global) j_glob = j_glob - Ny_global
    y_glob = (real(j_glob, dp) - 0.5_dp) * dy_grid

    ! Halo Ouest
    do i_loc = 1-nh, 0
      i_glob = i_start + i_loc - 1
      if (i_glob < 1) i_glob = i_glob + Nx_global
      x_glob = (real(i_glob, dp) - 0.5_dp) * dx_grid
      phi_exact = sin(2.0_dp * pi * x_glob / Lx) * cos(2.0_dp * pi * y_glob / Ly)
      err_local = max(err_local, abs(phi(i_loc, j_loc) - phi_exact))
    end do

    ! Halo Est
    do i_loc = nx_loc+1, nx_loc+nh
      i_glob = i_start + i_loc - 1
      if (i_glob > Nx_global) i_glob = i_glob - Nx_global
      x_glob = (real(i_glob, dp) - 0.5_dp) * dx_grid
      phi_exact = sin(2.0_dp * pi * x_glob / Lx) * cos(2.0_dp * pi * y_glob / Ly)
      err_local = max(err_local, abs(phi(i_loc, j_loc) - phi_exact))
    end do
  end do

  ! Vérification des halos Sud et Nord (intérieurs en X)
  do i_loc = 1, nx_loc
    i_glob = i_start + i_loc - 1
    x_glob = (real(i_glob, dp) - 0.5_dp) * dx_grid

    ! Halo Sud
    do j_loc = 1-nh, 0
      j_glob = j_start + j_loc - 1
      if (j_glob < 1) j_glob = j_glob + Ny_global
      y_glob = (real(j_glob, dp) - 0.5_dp) * dy_grid
      phi_exact = sin(2.0_dp * pi * x_glob / Lx) * cos(2.0_dp * pi * y_glob / Ly)
      err_local = max(err_local, abs(phi(i_loc, j_loc) - phi_exact))
    end do

    ! Halo Nord
    do j_loc = ny_loc+1, ny_loc+nh
      j_glob = j_start + j_loc - 1
      if (j_glob > Ny_global) j_glob = j_glob - Ny_global
      y_glob = (real(j_glob, dp) - 0.5_dp) * dy_grid
      phi_exact = sin(2.0_dp * pi * x_glob / Lx) * cos(2.0_dp * pi * y_glob / Ly)
      err_local = max(err_local, abs(phi(i_loc, j_loc) - phi_exact))
    end do
  end do

  ! 6. Réduction de l'erreur maximale globale
  call MPI_Allreduce(err_local, err_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_cart, ierr)

  if (my_rank == 0) then
    write(*,'(A,ES13.5)') " Erreur maximale dans les cellules fantomes : ", err_global
    if (err_global < 1.0e-13_dp) then
      write(*,'(A)') " --> [SUCCES] Test valide a la precision machine (L_inf < 1e-13) !"
    else
      write(*,'(A)') " --> [ECHEC] Incoherence detectee dans le transfert de halos."
    end if
    write(*,'(A)') "================================================================"
  end if

  deallocate(phi)
  call mpi_finalize_env()

end program test_halo_exchange
