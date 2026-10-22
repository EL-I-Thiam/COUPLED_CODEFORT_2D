module mod_mpi
  use mod_precision, only: dp, ip, zero, one
  use mpi
  implicit none
  public

  ! =========================================================================
  ! Variables globales MPI et Topologie Cartésienne 2D
  ! =========================================================================
  integer(ip) :: my_rank, n_procs
  integer(ip) :: comm_cart = MPI_COMM_NULL
  integer(ip) :: dims(2)
  integer(ip) :: coords(2)
  logical     :: periods(2) = [.false., .false.]

  ! Voisins cartésiens (Ouest, Est, Sud, Nord)
  integer(ip) :: nbr_west, nbr_east
  integer(ip) :: nbr_south, nbr_north

  ! Tailles locales et bornes globales du sous-domaine
  integer(ip) :: nx_loc, ny_loc
  integer(ip) :: i_start, i_end
  integer(ip) :: j_start, j_end

  ! Métrique locale
  real(dp) :: dx_grid, dy_grid

contains

  ! -------------------------------------------------------------------------
  ! Initialisation MPI
  ! -------------------------------------------------------------------------
  subroutine mpi_init_env()
    integer(ip) :: ierr
    logical :: is_init

    call MPI_Initialized(is_init, ierr)
    if (.not. is_init) then
      call MPI_Init(ierr)
    end if
    call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, n_procs, ierr)
  end subroutine mpi_init_env

  ! -------------------------------------------------------------------------
  ! Finalisation MPI
  ! -------------------------------------------------------------------------
  subroutine mpi_finalize_env()
    integer(ip) :: ierr
    if (comm_cart /= MPI_COMM_NULL .and. comm_cart /= MPI_COMM_WORLD) then
      call MPI_Comm_free(comm_cart, ierr)
      comm_cart = MPI_COMM_NULL
    end if
    call MPI_Finalize(ierr)
  end subroutine mpi_finalize_env

  ! -------------------------------------------------------------------------
  ! Configuration de la décomposition de domaine cartésienne 2D
  ! -------------------------------------------------------------------------
  subroutine mpi_setup_cartesian_mesh(Nx_global, Ny_global, Lx, Ly, periodic_x, periodic_y)
    integer(ip), intent(in) :: Nx_global, Ny_global
    real(dp),    intent(in) :: Lx, Ly
    logical, optional, intent(in) :: periodic_x, periodic_y
    integer(ip) :: ierr
    integer(ip) :: rem_x, rem_y

    if (present(periodic_x)) periods(1) = periodic_x
    if (present(periodic_y)) periods(2) = periodic_y

    ! Calcul automatique de la grille Px x Py
    dims = 0
    call MPI_Dims_create(n_procs, 2, dims, ierr)

    ! Création du communicateur cartésien
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, .true., comm_cart, ierr)
    call MPI_Comm_rank(comm_cart, my_rank, ierr)
    call MPI_Cart_coords(comm_cart, my_rank, 2, coords, ierr)

    ! Identification des voisins
    call MPI_Cart_shift(comm_cart, 0, 1, nbr_west, nbr_east, ierr)
    call MPI_Cart_shift(comm_cart, 1, 1, nbr_south, nbr_north, ierr)

    ! Découpage régulier équilibré
    nx_loc = Nx_global / dims(1)
    rem_x  = mod(Nx_global, dims(1))
    if (coords(1) < rem_x) then
      nx_loc = nx_loc + 1
      i_start = coords(1) * nx_loc + 1
    else
      i_start = rem_x * (nx_loc + 1) + (coords(1) - rem_x) * nx_loc + 1
    end if
    i_end = i_start + nx_loc - 1

    ny_loc = Ny_global / dims(2)
    rem_y  = mod(Ny_global, dims(2))
    if (coords(2) < rem_y) then
      ny_loc = ny_loc + 1
      j_start = coords(2) * ny_loc + 1
    else
      j_start = rem_y * (ny_loc + 1) + (coords(2) - rem_y) * ny_loc + 1
    end if
    j_end = j_start + ny_loc - 1

    dx_grid = Lx / real(Nx_global, dp)
    dy_grid = Ly / real(Ny_global, dp)

  end subroutine mpi_setup_cartesian_mesh

  ! -------------------------------------------------------------------------
  ! Échange de cellules fantômes pour scalaires centrés (p, C, rho)
  ! Domaine tableau : (1-nh : nx+nh, 1-nh : ny+nh)
  ! -------------------------------------------------------------------------
  subroutine mpi_exchange_halos_scalar(phi, nx, ny, nh)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: phi
    integer(ip) :: ierr
    integer(ip) :: status_arr(MPI_STATUS_SIZE)
    real(dp), allocatable :: s_buf_x(:,:), r_buf_x(:,:)
    real(dp), allocatable :: s_buf_y(:,:), r_buf_y(:,:)

    allocate(s_buf_x(nh, 1-nh:ny+nh), r_buf_x(nh, 1-nh:ny+nh))
    allocate(s_buf_y(1-nh:nx+nh, nh), r_buf_y(1-nh:nx+nh, nh))

    ! 1. Direction X : Est -> Ouest
    s_buf_x = phi(nx-nh+1:nx, :)
    call MPI_Sendrecv(s_buf_x, nh * (ny + 2*nh), MPI_DOUBLE_PRECISION, nbr_east, 101, &
                      r_buf_x, nh * (ny + 2*nh), MPI_DOUBLE_PRECISION, nbr_west, 101, &
                      comm_cart, status_arr, ierr)
    if (nbr_west /= MPI_PROC_NULL) phi(1-nh:0, :) = r_buf_x

    ! 2. Direction X : Ouest -> Est
    s_buf_x = phi(1:nh, :)
    call MPI_Sendrecv(s_buf_x, nh * (ny + 2*nh), MPI_DOUBLE_PRECISION, nbr_west, 102, &
                      r_buf_x, nh * (ny + 2*nh), MPI_DOUBLE_PRECISION, nbr_east, 102, &
                      comm_cart, status_arr, ierr)
    if (nbr_east /= MPI_PROC_NULL) phi(nx+1:nx+nh, :) = r_buf_x

    ! 3. Direction Y : Nord -> Sud
    s_buf_y = phi(:, ny-nh+1:ny)
    call MPI_Sendrecv(s_buf_y, (nx + 2*nh) * nh, MPI_DOUBLE_PRECISION, nbr_north, 201, &
                      r_buf_y, (nx + 2*nh) * nh, MPI_DOUBLE_PRECISION, nbr_south, 201, &
                      comm_cart, status_arr, ierr)
    if (nbr_south /= MPI_PROC_NULL) phi(:, 1-nh:0) = r_buf_y

    ! 4. Direction Y : Sud -> Nord
    s_buf_y = phi(:, 1:nh)
    call MPI_Sendrecv(s_buf_y, (nx + 2*nh) * nh, MPI_DOUBLE_PRECISION, nbr_south, 202, &
                      r_buf_y, (nx + 2*nh) * nh, MPI_DOUBLE_PRECISION, nbr_north, 202, &
                      comm_cart, status_arr, ierr)
    if (nbr_north /= MPI_PROC_NULL) phi(:, ny+1:ny+nh) = r_buf_y

    deallocate(s_buf_x, r_buf_x, s_buf_y, r_buf_y)
  end subroutine mpi_exchange_halos_scalar

  ! -------------------------------------------------------------------------
  ! Échange de cellules fantômes pour vitesse u (faces verticales)
  ! Taille tableau : (0-nh : nx+nh, 1-nh : ny+nh)
  ! -------------------------------------------------------------------------
  subroutine mpi_exchange_halos_u(u_vel, nx, ny, nh)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: u_vel
    integer(ip) :: ierr
    integer(ip) :: status_arr(MPI_STATUS_SIZE)
    real(dp), allocatable :: s_buf_x(:,:), r_buf_x(:,:)
    real(dp), allocatable :: s_buf_y(:,:), r_buf_y(:,:)

    allocate(s_buf_x(nh, 1-nh:ny+nh), r_buf_x(nh, 1-nh:ny+nh))
    allocate(s_buf_y(0-nh:nx+nh, nh), r_buf_y(0-nh:nx+nh, nh))

    ! Envoi vers l'Est / Réception de l'Ouest
    s_buf_x = u_vel(nx-nh:nx-1, :)
    call MPI_Sendrecv(s_buf_x, nh * (ny + 2*nh), MPI_DOUBLE_PRECISION, nbr_east, 301, &
                      r_buf_x, nh * (ny + 2*nh), MPI_DOUBLE_PRECISION, nbr_west, 301, &
                      comm_cart, status_arr, ierr)
    if (nbr_west /= MPI_PROC_NULL) u_vel(-nh:-1, :) = r_buf_x

    ! Envoi vers l'Ouest / Réception de l'Est
    s_buf_x = u_vel(1:nh, :)
    call MPI_Sendrecv(s_buf_x, nh * (ny + 2*nh), MPI_DOUBLE_PRECISION, nbr_west, 302, &
                      r_buf_x, nh * (ny + 2*nh), MPI_DOUBLE_PRECISION, nbr_east, 302, &
                      comm_cart, status_arr, ierr)
    if (nbr_east /= MPI_PROC_NULL) u_vel(nx+1:nx+nh, :) = r_buf_x

    ! Envoi vers le Nord / Réception du Sud
    s_buf_y = u_vel(:, ny-nh+1:ny)
    call MPI_Sendrecv(s_buf_y, (nx + 1 + 2*nh) * nh, MPI_DOUBLE_PRECISION, nbr_north, 401, &
                      r_buf_y, (nx + 1 + 2*nh) * nh, MPI_DOUBLE_PRECISION, nbr_south, 401, &
                      comm_cart, status_arr, ierr)
    if (nbr_south /= MPI_PROC_NULL) u_vel(:, 1-nh:0) = r_buf_y

    ! Envoi vers le Sud / Réception du Nord
    s_buf_y = u_vel(:, 1:nh)
    call MPI_Sendrecv(s_buf_y, (nx + 1 + 2*nh) * nh, MPI_DOUBLE_PRECISION, nbr_south, 402, &
                      r_buf_y, (nx + 1 + 2*nh) * nh, MPI_DOUBLE_PRECISION, nbr_north, 402, &
                      comm_cart, status_arr, ierr)
    if (nbr_north /= MPI_PROC_NULL) u_vel(:, ny+1:ny+nh) = r_buf_y

    deallocate(s_buf_x, r_buf_x, s_buf_y, r_buf_y)
  end subroutine mpi_exchange_halos_u

  ! -------------------------------------------------------------------------
  ! Échange de cellules fantômes pour vitesse v (faces horizontales)
  ! Taille tableau : (1-nh : nx+nh, 0-nh : ny+nh)
  ! -------------------------------------------------------------------------
  subroutine mpi_exchange_halos_v(v_vel, nx, ny, nh)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(inout) :: v_vel
    integer(ip) :: ierr
    integer(ip) :: status_arr(MPI_STATUS_SIZE)
    real(dp), allocatable :: s_buf_x(:,:), r_buf_x(:,:)
    real(dp), allocatable :: s_buf_y(:,:), r_buf_y(:,:)

    allocate(s_buf_x(nh, 0-nh:ny+nh), r_buf_x(nh, 0-nh:ny+nh))
    allocate(s_buf_y(1-nh:nx+nh, nh), r_buf_y(1-nh:nx+nh, nh))

    ! Envoi vers l'Est / Réception de l'Ouest
    s_buf_x = v_vel(nx-nh+1:nx, :)
    call MPI_Sendrecv(s_buf_x, nh * (ny + 1 + 2*nh), MPI_DOUBLE_PRECISION, nbr_east, 501, &
                      r_buf_x, nh * (ny + 1 + 2*nh), MPI_DOUBLE_PRECISION, nbr_west, 501, &
                      comm_cart, status_arr, ierr)
    if (nbr_west /= MPI_PROC_NULL) v_vel(1-nh:0, :) = r_buf_x

    ! Envoi vers l'Ouest / Réception de l'Est
    s_buf_x = v_vel(1:nh, :)
    call MPI_Sendrecv(s_buf_x, nh * (ny + 1 + 2*nh), MPI_DOUBLE_PRECISION, nbr_west, 502, &
                      r_buf_x, nh * (ny + 1 + 2*nh), MPI_DOUBLE_PRECISION, nbr_east, 502, &
                      comm_cart, status_arr, ierr)
    if (nbr_east /= MPI_PROC_NULL) v_vel(nx+1:nx+nh, :) = r_buf_x

    ! Envoi vers le Nord / Réception du Sud
    s_buf_y = v_vel(:, ny-nh:ny-1)
    call MPI_Sendrecv(s_buf_y, (nx + 2*nh) * nh, MPI_DOUBLE_PRECISION, nbr_north, 601, &
                      r_buf_y, (nx + 2*nh) * nh, MPI_DOUBLE_PRECISION, nbr_south, 601, &
                      comm_cart, status_arr, ierr)
    if (nbr_south /= MPI_PROC_NULL) v_vel(:, -nh:-1) = r_buf_y

    ! Envoi vers le Sud / Réception du Nord
    s_buf_y = v_vel(:, 1:nh)
    call MPI_Sendrecv(s_buf_y, (nx + 2*nh) * nh, MPI_DOUBLE_PRECISION, nbr_south, 602, &
                      r_buf_y, (nx + 2*nh) * nh, MPI_DOUBLE_PRECISION, nbr_north, 602, &
                      comm_cart, status_arr, ierr)
    if (nbr_north /= MPI_PROC_NULL) v_vel(:, ny+1:ny+nh) = r_buf_y

    deallocate(s_buf_x, r_buf_x, s_buf_y, r_buf_y)
  end subroutine mpi_exchange_halos_v

  ! -------------------------------------------------------------------------
  ! Réduction globale pour le produit scalaire distribué (BiCGStab2)
  ! -------------------------------------------------------------------------
  function mpi_global_dot_product(vec_loc_a, vec_loc_b, n_elem) result(dot_glob)
    integer(ip), intent(in) :: n_elem
    real(dp), dimension(n_elem), intent(in) :: vec_loc_a, vec_loc_b
    real(dp) :: dot_loc, dot_glob
    integer(ip) :: ierr

    dot_loc = dot_product(vec_loc_a, vec_loc_b)
    call MPI_Allreduce(dot_loc, dot_glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_cart, ierr)
  end function mpi_global_dot_product

  ! -------------------------------------------------------------------------
  ! Réduction globale pour le pas de temps critique (CFL global)
  ! -------------------------------------------------------------------------
  function mpi_global_min_dt(dt_local) result(dt_glob)
    real(dp), intent(in) :: dt_local
    real(dp) :: dt_glob
    integer(ip) :: ierr

    call MPI_Allreduce(dt_local, dt_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm_cart, ierr)
  end function mpi_global_min_dt

end module mod_mpi
