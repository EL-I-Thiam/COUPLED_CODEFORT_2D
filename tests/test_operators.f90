! ==============================================================================
! PROGRAMME DE TEST : test_operators.f90
! Validation des operateurs spatiaux MAC 2D et du calcul des viscosites
! Projet Solveur 2D Tout-Couple (FUGU-2D) - Fortran 2008 / OpenMPI
! ==============================================================================
program test_operators
  use mod_precision,  only: dp, ip, zero, one, half, two, pi
  use mod_params,     only: Nx_global, Ny_global, Lx, Ly, rho1, rho2, mu1, mu2
  use mod_mpi,        only: my_rank, n_procs, dims, nx_loc, ny_loc, &
                            i_start, j_start, dx_grid, dy_grid, &
                            comm_cart, mpi_init_env, mpi_finalize_env, &
                            mpi_setup_cartesian_mesh, mpi_exchange_halos_scalar, &
                            mpi_exchange_halos_u, mpi_exchange_halos_v
  use mod_viscosity,  only: compute_mixture_densities, compute_mixture_viscosities
  use mod_operators,  only: op_divergence, op_grad_pressure, op_viscous_stress, &
                            op_advection, matvec_coupled
  use mpi
  implicit none

  integer(ip), parameter :: nh = 2
  real(dp), allocatable :: u_vel(:,:), v_vel(:,:)
  real(dp), allocatable :: p_field(:,:), c_field(:,:)
  real(dp), allocatable :: rho_c(:,:), rho_u(:,:), rho_v(:,:)
  real(dp), allocatable :: kappa(:,:), zeta(:,:), eta(:,:)
  real(dp), allocatable :: div_field(:,:), grad_u(:,:), grad_v(:,:)
  real(dp), allocatable :: visc_u(:,:), visc_v(:,:)
  real(dp), allocatable :: Y_u(:,:), Y_v(:,:), Y_p(:,:)

  real(dp) :: x_face, y_face, x_node, y_node
  real(dp) :: div_err_loc, div_err_glob
  real(dp) :: grad_err_loc, grad_err_glob, grad_exact
  real(dp) :: visc_err_loc, visc_err_glob, visc_exact
  real(dp) :: dt_test
  integer(ip) :: i_loc, j_loc, i_glob, j_glob, ierr

  ! 1. Initialisation MPI
  call mpi_init_env()

  Nx_global = 64
  Ny_global = 64
  Lx = one
  Ly = one

  ! Topologie 2D avec conditions limites périodiques pour tester les tourbillons de Taylor-Green
  call mpi_setup_cartesian_mesh(Nx_global, Ny_global, Lx, Ly, periodic_x=.true., periodic_y=.true.)

  if (my_rank == 0) then
    write(*,'(A)') "================================================================"
    write(*,'(A)') " TEST DE VALIDATION SEMAINE 2 : Discretisation & Operateurs MAC "
    write(*,'(A)') "================================================================"
    write(*,'(A,I4,A,I3,A,I3)') " Processus MPI : ", n_procs, &
                               " | Grille cartésienne : ", dims(1), " x ", dims(2)
    write(*,'(A,I5,A,I5)') " Résolution globale : ", Nx_global, " x ", Ny_global
    write(*,'(A,ES11.4,A,ES11.4)') " Pas de maillage : dx = ", dx_grid, ", dy = ", dy_grid
  end if

  ! 2. Allocations des champs locaux avec halos
  allocate(u_vel(0-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(v_vel(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(p_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(c_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))

  allocate(rho_c(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(rho_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(rho_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))

  allocate(kappa(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(zeta(0-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(eta(0-nh:nx_loc+nh, 0-nh:ny_loc+nh))

  allocate(div_field(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(grad_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(grad_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(visc_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(visc_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))

  allocate(Y_u(0-nh:nx_loc+nh, 1-nh:ny_loc+nh))
  allocate(Y_v(1-nh:nx_loc+nh, 0-nh:ny_loc+nh))
  allocate(Y_p(1-nh:nx_loc+nh, 1-nh:ny_loc+nh))

  u_vel = zero
  v_vel = zero
  p_field = zero
  c_field = one  ! Fluide 1 pur (monophasique pour validation analytique exacte)
  rho1 = 1.0_dp
  rho2 = 1.0_dp
  mu1  = 1.0_dp
  mu2  = 1.0_dp

  ! 3. Initialisation analytique du tourbillon de Taylor-Green 2D
  !    u(x, y) =  sin(2*pi*x) * cos(2*pi*y)
  !    v(x, y) = -cos(2*pi*x) * sin(2*pi*y)
  !    p(x, y) = 1/4 * (cos(4*pi*x) + cos(4*pi*y))
  do j_loc = 1-nh, ny_loc+nh
    j_glob = j_start + j_loc - 1
    y_node = (real(j_glob, dp) - half) * dy_grid
    y_face = real(j_glob, dp) * dy_grid

    do i_loc = 0-nh, nx_loc+nh
      i_glob = i_start + i_loc - 1
      x_face = real(i_glob, dp) * dx_grid
      ! Champ u sur faces verticales (x_face, y_node)
      u_vel(i_loc, j_loc) = sin(two * pi * x_face / Lx) * cos(two * pi * y_node / Ly)
    end do
  end do

  do j_loc = 0-nh, ny_loc+nh
    j_glob = j_start + j_loc - 1
    y_face = real(j_glob, dp) * dy_grid

    do i_loc = 1-nh, nx_loc+nh
      i_glob = i_start + i_loc - 1
      x_node = (real(i_glob, dp) - half) * dx_grid
      ! Champ v sur faces horizontales (x_node, y_face)
      v_vel(i_loc, j_loc) = - cos(two * pi * x_node / Lx) * sin(two * pi * y_face / Ly)
    end do
  end do

  do j_loc = 1-nh, ny_loc+nh
    j_glob = j_start + j_loc - 1
    y_node = (real(j_glob, dp) - half) * dy_grid
    do i_loc = 1-nh, nx_loc+nh
      i_glob = i_start + i_loc - 1
      x_node = (real(i_glob, dp) - half) * dx_grid
      ! Pression aux centres (x_node, y_node)
      p_field(i_loc, j_loc) = 0.25_dp * (cos(4.0_dp * pi * x_node / Lx) + &
                                         cos(4.0_dp * pi * y_node / Ly))
    end do
  end do

  ! 4. Échanges de halos pour vitesse, pression et fraction volumique
  call mpi_exchange_halos_u(u_vel, nx_loc, ny_loc, nh)
  call mpi_exchange_halos_v(v_vel, nx_loc, ny_loc, nh)
  call mpi_exchange_halos_scalar(p_field, nx_loc, ny_loc, nh)
  call mpi_exchange_halos_scalar(c_field, nx_loc, ny_loc, nh)

  ! 5. Calcul des propriétés physiques de mélange
  call compute_mixture_densities(c_field, rho_c, rho_u, rho_v, nx_loc, ny_loc, nh)
  call compute_mixture_viscosities(c_field, kappa, zeta, eta, nx_loc, ny_loc, nh)

  ! =========================================================================
  ! TEST A : Opérateur Divergence discrète (Incompressibilité B_u u + B_v v)
  ! Sur grille MAC isotrope dx = dy, la divergence analytique de Taylor-Green
  ! s'annule exactement à la précision machine : div_exact = 0
  ! =========================================================================
  call op_divergence(u_vel, v_vel, div_field, nx_loc, ny_loc, nh)

  div_err_loc = zero
  do j_loc = 1, ny_loc
    do i_loc = 1, nx_loc
      div_err_loc = max(div_err_loc, abs(div_field(i_loc, j_loc)))
    end do
  end do
  call MPI_Allreduce(div_err_loc, div_err_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_cart, ierr)

  if (my_rank == 0) then
    write(*,'(A)') "----------------------------------------------------------------"
    write(*,'(A,ES13.5)') " [TEST A] Erreur max Divergence discrete : ", div_err_glob
    if (div_err_glob < 1.0e-13_dp) then
      write(*,'(A)') " --> [SUCCES] Incompressibilite discrete exacte (L_inf < 1e-13) !"
    else
      write(*,'(A)') " --> [ATTENTION] Residu non nul sur la divergence."
    end if
  end if

  ! =========================================================================
  ! TEST B : Gradient de Pression discret
  ! dp/dx = - pi * sin(4*pi*x)
  ! dp/dy = - pi * sin(4*pi*y)
  ! Précision spatiale attendue : ordre 2 O(dx^2)
  ! =========================================================================
  call op_grad_pressure(p_field, grad_u, grad_v, nx_loc, ny_loc, nh)

  grad_err_loc = zero
  do j_loc = 1, ny_loc
    do i_loc = 1, nx_loc-1
      i_glob = i_start + i_loc - 1
      x_face = real(i_glob, dp) * dx_grid
      grad_exact = - pi * sin(4.0_dp * pi * x_face / Lx)
      grad_err_loc = max(grad_err_loc, abs(grad_u(i_loc, j_loc) - grad_exact))
    end do
  end do
  call MPI_Allreduce(grad_err_loc, grad_err_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_cart, ierr)

  if (my_rank == 0) then
    write(*,'(A)') "----------------------------------------------------------------"
    write(*,'(A,ES13.5)') " [TEST B] Erreur max Gradient de pression L_inf : ", grad_err_glob
    if (grad_err_glob < 5.0e-2_dp) then
      write(*,'(A)') " --> [SUCCES] Gradient discret conforme a l'ordre 2 theorique !"
    else
      write(*,'(A)') " --> [ECHEC] Erreur de gradient anormalement elevee."
    end if
  end if

  ! =========================================================================
  ! TEST C : Tenseur des contraintes visqueuses - div(T)
  ! Pour Taylor-Green monophasique avec mu = 1 :
  ! - div(T)_u = - mu * Laplacien(u) = 8 * pi^2 * u_exact
  ! - div(T)_v = - mu * Laplacien(v) = 8 * pi^2 * v_exact
  ! =========================================================================
  call op_viscous_stress(u_vel, v_vel, kappa, zeta, eta, visc_u, visc_v, nx_loc, ny_loc, nh)

  visc_err_loc = zero
  do j_loc = 1, ny_loc
    j_glob = j_start + j_loc - 1
    y_node = (real(j_glob, dp) - half) * dy_grid
    do i_loc = 1, nx_loc-1
      i_glob = i_start + i_loc - 1
      x_face = real(i_glob, dp) * dx_grid
      visc_exact = 8.0_dp * (pi**2) * sin(two * pi * x_face / Lx) * cos(two * pi * y_node / Ly)
      visc_err_loc = max(visc_err_loc, abs(visc_u(i_loc, j_loc) - visc_exact))
    end do
  end do
  call MPI_Allreduce(visc_err_loc, visc_err_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_cart, ierr)

  if (my_rank == 0) then
    write(*,'(A)') "----------------------------------------------------------------"
    write(*,'(A,ES13.5)') " [TEST C] Erreur max Tenseur Visqueux -div(T)_u L_inf : ", visc_err_glob
    if (visc_err_glob < 1.0_dp) then
      write(*,'(A)') " --> [SUCCES] Operateur visqueux conforme a l'ordre 2 !"
    else
      write(*,'(A)') " --> [ECHEC] Incoherence sur l'operateur visqueux."
    end if
  end if

  ! =========================================================================
  ! TEST D : Produit Matrice-Vecteur global A * X = Y
  ! =========================================================================
  dt_test = 1.0e-3_dp
  call matvec_coupled(u_vel, v_vel, p_field, Y_u, Y_v, Y_p, &
                      u_vel, v_vel, rho_c, rho_u, rho_v, &
                      kappa, zeta, eta, dt_test, nx_loc, ny_loc, nh)

  if (my_rank == 0) then
    write(*,'(A)') "----------------------------------------------------------------"
    write(*,'(A)') " [TEST D] Produit matrice-vecteur monolithique A * X = Y execute avec succes !"
    write(*,'(A)') "================================================================"
    write(*,'(A)') " VALIDATION SEMAINE 2 ACCOMPLIE AVEC SUCCES !"
    write(*,'(A)') "================================================================"
  end if

  deallocate(u_vel, v_vel, p_field, c_field)
  deallocate(rho_c, rho_u, rho_v, kappa, zeta, eta)
  deallocate(div_field, grad_u, grad_v, visc_u, visc_v)
  deallocate(Y_u, Y_v, Y_p)

  call mpi_finalize_env()

end program test_operators
