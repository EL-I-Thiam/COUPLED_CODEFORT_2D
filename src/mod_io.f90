! ==============================================================================
! MODULE : mod_io
! Solveur 2D Tout-Couplé (FUGU-2D) - Fortran 2008 / OpenMPI
!
! Entrées / Sorties parallèles :
!   - Exportation au format VTK Legacy (.vtk STRUCTURED_POINTS) par rassemblement
!     sur le processus maître (rang 0) pour visualisation sous ParaView.
! ==============================================================================
module mod_io
  use mod_precision, only: dp, ip, zero, half
  use mod_params,    only: Nx_global, Ny_global
  use mod_mpi,       only: my_rank, n_procs, nx_loc, ny_loc, i_start, i_end, &
                           j_start, j_end, dx_grid, dy_grid, comm_cart
  use mpi
  implicit none
  private

  public :: write_vtk_fields
  public :: write_checkpoint
  public :: read_checkpoint

contains

  ! ----------------------------------------------------------------------------
  ! Rassemble les champs locaux sur le rang 0 et écrit un fichier VTK Legacy
  ! Contient : Pression (p), Fraction volumique (C), Masse volumique (rho),
  !            Courbure (kappa), et Vitesse centrée (u, v, 0).
  ! ----------------------------------------------------------------------------
  subroutine write_vtk_fields(filename, p_field, c_field, rho_c, kappa_curv, &
                              u_vel, v_vel, nx, ny, nh)
    character(len=*), intent(in) :: filename
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in) :: p_field, c_field, rho_c, kappa_curv
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in) :: u_vel
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in) :: v_vel

    real(dp), allocatable :: uc_loc(:,:), vc_loc(:,:)
    real(dp), allocatable :: p_glob(:,:), c_glob(:,:), rho_glob(:,:), kap_glob(:,:)
    real(dp), allocatable :: u_glob(:,:), v_glob(:,:)
    real(dp), allocatable :: buf_send(:), buf_recv(:)
    integer(ip) :: bounds(4), r_bounds(4)
    integer(ip) :: i, j, idx, src, n_loc_pts, n_recv_pts, ierr
    integer(ip) :: status_arr(MPI_STATUS_SIZE)
    integer(ip) :: unit_vtk

    allocate(uc_loc(1:nx, 1:ny), vc_loc(1:nx, 1:ny))
    do j = 1, ny
      do i = 1, nx
        uc_loc(i, j) = half * (u_vel(i-1, j) + u_vel(i, j))
        vc_loc(i, j) = half * (v_vel(i, j-1) + v_vel(i, j))
      end do
    end do

    n_loc_pts = nx * ny
    allocate(buf_send(6 * n_loc_pts))
    idx = 0
    do j = 1, ny
      do i = 1, nx
        idx = idx + 1
        buf_send(idx)                 = p_field(i, j)
        buf_send(idx +   n_loc_pts)   = c_field(i, j)
        buf_send(idx + 2*n_loc_pts)   = rho_c(i, j)
        buf_send(idx + 3*n_loc_pts)   = kappa_curv(i, j)
        buf_send(idx + 4*n_loc_pts)   = uc_loc(i, j)
        buf_send(idx + 5*n_loc_pts)   = vc_loc(i, j)
      end do
    end do

    bounds = [i_start, i_end, j_start, j_end]

    if (my_rank == 0) then
      allocate(p_glob(Nx_global, Ny_global),   c_glob(Nx_global, Ny_global))
      allocate(rho_glob(Nx_global, Ny_global), kap_glob(Nx_global, Ny_global))
      allocate(u_glob(Nx_global, Ny_global),   v_glob(Nx_global, Ny_global))

      ! Copie du bloc local du rang 0
      p_glob(i_start:i_end, j_start:j_end)   = p_field(1:nx, 1:ny)
      c_glob(i_start:i_end, j_start:j_end)   = c_field(1:nx, 1:ny)
      rho_glob(i_start:i_end, j_start:j_end) = rho_c(1:nx, 1:ny)
      kap_glob(i_start:i_end, j_start:j_end) = kappa_curv(1:nx, 1:ny)
      u_glob(i_start:i_end, j_start:j_end)   = uc_loc(1:nx, 1:ny)
      v_glob(i_start:i_end, j_start:j_end)   = vc_loc(1:nx, 1:ny)

      ! Réception des autres rangs
      do src = 1, n_procs - 1
        call MPI_Recv(r_bounds, 4, MPI_INTEGER, src, 901, comm_cart, status_arr, ierr)
        n_recv_pts = (r_bounds(2) - r_bounds(1) + 1) * (r_bounds(4) - r_bounds(3) + 1)
        allocate(buf_recv(6 * n_recv_pts))
        call MPI_Recv(buf_recv, 6 * n_recv_pts, MPI_DOUBLE_PRECISION, src, 902, comm_cart, status_arr, ierr)

        idx = 0
        do j = r_bounds(3), r_bounds(4)
          do i = r_bounds(1), r_bounds(2)
            idx = idx + 1
            p_glob(i, j)   = buf_recv(idx)
            c_glob(i, j)   = buf_recv(idx +   n_recv_pts)
            rho_glob(i, j) = buf_recv(idx + 2*n_recv_pts)
            kap_glob(i, j) = buf_recv(idx + 3*n_recv_pts)
            u_glob(i, j)   = buf_recv(idx + 4*n_recv_pts)
            v_glob(i, j)   = buf_recv(idx + 5*n_recv_pts)
          end do
        end do
        deallocate(buf_recv)
      end do

      ! Écriture du fichier VTK Legacy
      open(newunit=unit_vtk, file=trim(filename), status='replace', action='write')
      write(unit_vtk, '(A)') '# vtk DataFile Version 3.0'
      write(unit_vtk, '(A)') 'FUGU-2D Coupled Two-Phase Solver Output'
      write(unit_vtk, '(A)') 'ASCII'
      write(unit_vtk, '(A)') 'DATASET STRUCTURED_POINTS'
      write(unit_vtk, '(A,3I8)') 'DIMENSIONS ', Nx_global, Ny_global, 1
      write(unit_vtk, '(A,3ES15.7)') 'ORIGIN ', half*dx_grid, half*dy_grid, zero
      write(unit_vtk, '(A,3ES15.7)') 'SPACING ', dx_grid, dy_grid, dx_grid
      write(unit_vtk, '(A,I10)') 'POINT_DATA ', Nx_global * Ny_global

      write(unit_vtk, '(A)') 'SCALARS VOF_C double 1'
      write(unit_vtk, '(A)') 'LOOKUP_TABLE default'
      do j = 1, Ny_global
        do i = 1, Nx_global
          write(unit_vtk, '(ES15.7)') c_glob(i, j)
        end do
      end do

      write(unit_vtk, '(A)') 'SCALARS Pressure double 1'
      write(unit_vtk, '(A)') 'LOOKUP_TABLE default'
      do j = 1, Ny_global
        do i = 1, Nx_global
          write(unit_vtk, '(ES15.7)') p_glob(i, j)
        end do
      end do

      write(unit_vtk, '(A)') 'SCALARS Density double 1'
      write(unit_vtk, '(A)') 'LOOKUP_TABLE default'
      do j = 1, Ny_global
        do i = 1, Nx_global
          write(unit_vtk, '(ES15.7)') rho_glob(i, j)
        end do
      end do

      write(unit_vtk, '(A)') 'SCALARS Curvature double 1'
      write(unit_vtk, '(A)') 'LOOKUP_TABLE default'
      do j = 1, Ny_global
        do i = 1, Nx_global
          write(unit_vtk, '(ES15.7)') kap_glob(i, j)
        end do
      end do

      write(unit_vtk, '(A)') 'VECTORS Velocity double'
      do j = 1, Ny_global
        do i = 1, Nx_global
          write(unit_vtk, '(3ES15.7)') u_glob(i, j), v_glob(i, j), zero
        end do
      end do

      close(unit_vtk)
      deallocate(p_glob, c_glob, rho_glob, kap_glob, u_glob, v_glob)
    else
      call MPI_Send(bounds, 4, MPI_INTEGER, 0, 901, comm_cart, ierr)
      call MPI_Send(buf_send, 6 * n_loc_pts, MPI_DOUBLE_PRECISION, 0, 902, comm_cart, ierr)
    end if

    deallocate(uc_loc, vc_loc, buf_send)
  end subroutine write_vtk_fields

  ! ----------------------------------------------------------------------------
  ! Sauvegarde binaire par rang MPI (Checkpoint / Reprise exacte bit-a-bit)
  ! ----------------------------------------------------------------------------
  subroutine write_checkpoint(prefix, step, t_curr, u_vel, v_vel, u_prev, v_prev, &
                              p_field, c_field, nx, ny, nh)
    character(len=*), intent(in) :: prefix
    integer(ip), intent(in) :: step, nx, ny, nh
    real(dp), intent(in) :: t_curr
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in) :: u_vel, u_prev
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in) :: v_vel, v_prev
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in) :: p_field, c_field

    character(len=256) :: fname
    integer(ip) :: unit_chk

    write(fname, '(A,A,I4.4,A)') trim(prefix), '_r', my_rank, '.bin'
    open(newunit=unit_chk, file=trim(fname), form='unformatted', access='stream', &
         status='replace', action='write')
    write(unit_chk) step, t_curr, nx, ny, nh
    write(unit_chk) u_vel, v_vel, u_prev, v_prev, p_field, c_field
    close(unit_chk)
  end subroutine write_checkpoint

  ! ----------------------------------------------------------------------------
  ! Lecture binaire par rang MPI pour reprise de calcul
  ! ----------------------------------------------------------------------------
  subroutine read_checkpoint(prefix, step, t_curr, u_vel, v_vel, u_prev, v_prev, &
                             p_field, c_field, nx, ny, nh, ok)
    character(len=*), intent(in) :: prefix
    integer(ip), intent(out) :: step
    real(dp), intent(out) :: t_curr
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(out) :: u_vel, u_prev
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(out) :: v_vel, v_prev
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(out) :: p_field, c_field
    logical, intent(out) :: ok

    character(len=256) :: fname
    integer(ip) :: unit_chk, nx_f, ny_f, nh_f
    logical :: exists_f

    write(fname, '(A,A,I4.4,A)') trim(prefix), '_r', my_rank, '.bin'
    inquire(file=trim(fname), exist=exists_f)
    if (.not. exists_f) then
      ok = .false.
      return
    end if

    open(newunit=unit_chk, file=trim(fname), form='unformatted', access='stream', &
         status='old', action='read')
    read(unit_chk) step, t_curr, nx_f, ny_f, nh_f
    if (nx_f /= nx .or. ny_f /= ny .or. nh_f /= nh) then
      close(unit_chk)
      ok = .false.
      return
    end if
    read(unit_chk) u_vel, v_vel, u_prev, v_prev, p_field, c_field
    close(unit_chk)
    ok = .true.
  end subroutine read_checkpoint

end module mod_io
