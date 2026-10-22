! ==============================================================================
! MODULE : mod_operators
! Solveur 2D Tout-Couple (FUGU-2D) - Fortran 2008 / OpenMPI
!
! Discretisation des operateurs spatiaux sur maillage cartésien décale (MAC) :
!   1. Divergence de la vitesse aux centres des mailles (B_u u + B_v v)
!   2. Gradient de pression discret sur les faces decalees (B_pu^T p, B_pv^T p)
!   3. Tenseur des contraintes visqueuses avec decomposition elongation/cisaillement/rotation
!   4. Convection semi-implicite centree d'ordre 2 avec champ advectant u_adv
!   5. Application des conditions aux limites (parois adherentes / entrainees / periodiques)
!   6. Produit matrice-vecteur monolithique global A * X = Y pour le systeme point-selle
! ==============================================================================
module mod_operators
  use mod_precision, only: dp, ip, zero, one, half, two
  use mod_mpi,       only: dx_grid, dy_grid, coords, dims, periods, &
                           mpi_exchange_halos_scalar, mpi_exchange_halos_u, &
                           mpi_exchange_halos_v
  implicit none
  private

  public :: op_divergence
  public :: op_grad_pressure
  public :: op_viscous_stress
  public :: op_advection
  public :: apply_velocity_bcs
  public :: matvec_coupled

contains

  ! ----------------------------------------------------------------------------
  ! 1. Operateur Divergence discret au centre des mailles
  !    div(u, v)_{i,j} = (u_{i,j} - u_{i-1,j}) / dx + (v_{i,j} - v_{i,j-1}) / dy
  ! ----------------------------------------------------------------------------
  subroutine op_divergence(u_vel, v_vel, div_field, nx, ny, nh)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)  :: u_vel
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)  :: v_vel
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(out) :: div_field

    integer(ip) :: i, j
    real(dp) :: inv_dx, inv_dy

    inv_dx = one / dx_grid
    inv_dy = one / dy_grid

    div_field = zero
    do j = 1, ny
      do i = 1, nx
        div_field(i, j) = (u_vel(i, j) - u_vel(i-1, j)) * inv_dx + &
                          (v_vel(i, j) - v_vel(i, j-1)) * inv_dy
      end do
    end do
  end subroutine op_divergence

  ! ----------------------------------------------------------------------------
  ! 2. Operateur Gradient de pression discret
  !    (grad p)_u sur faces verticales : (p_{i+1,j} - p_{i,j}) / dx
  !    (grad p)_v sur faces horizontales : (p_{i,j+1} - p_{i,j}) / dy
  ! ----------------------------------------------------------------------------
  subroutine op_grad_pressure(p_field, grad_u, grad_v, nx, ny, nh)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)  :: p_field
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(out) :: grad_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(out) :: grad_v

    integer(ip) :: i, j
    real(dp) :: inv_dx, inv_dy

    inv_dx = one / dx_grid
    inv_dy = one / dy_grid

    grad_u = zero
    grad_v = zero

    ! Faces u : i = 0 .. nx (faces verticales internes et frontieres locales)
    do j = 1, ny
      do i = 0, nx
        grad_u(i, j) = (p_field(i+1, j) - p_field(i, j)) * inv_dx
      end do
    end do

    ! Faces v : j = 0 .. ny (faces horizontales internes et frontieres locales)
    do j = 0, ny
      do i = 1, nx
        grad_v(i, j) = (p_field(i, j+1) - p_field(i, j)) * inv_dy
      end do
    end do
  end subroutine op_grad_pressure

  ! ----------------------------------------------------------------------------
  ! 3. Operateur Tenseur des Contraintes Visqueuses : - div(T)
  !    Discretisation par volumes finis avec viscosites decomposees :
  !      - T_xx aux centres = kappa * (u_{i,j} - u_{i-1,j}) / dx
  !      - T_yy aux centres = kappa * (v_{i,j} - v_{i,j-1}) / dy
  !      - T_xy aux sommets = (zeta - eta) * (u_{i,j+1} - u_{i,j})/dy + eta * (v_{i+1,j} - v_{i,j})/dx
  !      - T_yx aux sommets = eta * (u_{i,j+1} - u_{i,j})/dy + (zeta - eta) * (v_{i+1,j} - v_{i,j})/dx
  !    (Avec zeta = 2*mu, eta = mu => zeta - eta = eta = mu)
  ! ----------------------------------------------------------------------------
  subroutine op_viscous_stress(u_vel, v_vel, kappa, zeta, eta, visc_u, visc_v, nx, ny, nh)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)  :: u_vel
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)  :: v_vel
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)  :: kappa
    real(dp), dimension(0-nh:nx+nh, 0-nh:ny+nh), intent(in)  :: zeta
    real(dp), dimension(0-nh:nx+nh, 0-nh:ny+nh), intent(in)  :: eta
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(out) :: visc_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(out) :: visc_v

    integer(ip) :: i, j
    real(dp) :: inv_dx, inv_dy
    real(dp), dimension(0:nx+1, 1:ny)   :: T_xx
    real(dp), dimension(1:nx,   0:ny+1) :: T_yy
    real(dp), dimension(0:nx,   0:ny)   :: T_xy, T_yx
    real(dp) :: du_dy, dv_dx

    inv_dx = one / dx_grid
    inv_dy = one / dy_grid

    ! a) Contraintes d'elongation T_xx et T_yy aux centres des cellules
    do j = 1, ny
      do i = 0, nx+1
        T_xx(i, j) = kappa(i, j) * (u_vel(i, j) - u_vel(i-1, j)) * inv_dx
      end do
    end do

    do j = 0, ny+1
      do i = 1, nx
        T_yy(i, j) = kappa(i, j) * (v_vel(i, j) - v_vel(i, j-1)) * inv_dy
      end do
    end do

    ! b) Contraintes de cisaillement et rotation T_xy et T_yx aux sommets
    do j = 0, ny
      do i = 0, nx
        du_dy = (u_vel(i, j+1) - u_vel(i, j)) * inv_dy
        dv_dx = (v_vel(i+1, j) - v_vel(i, j)) * inv_dx
        T_xy(i, j) = (zeta(i, j) - eta(i, j)) * du_dy + eta(i, j) * dv_dx
        T_yx(i, j) = eta(i, j) * du_dy + (zeta(i, j) - eta(i, j)) * dv_dx
      end do
    end do

    ! c) Divergence pour la composante u (sur faces u : i = 0 .. nx)
    !    - div(T)_u = - ( (T_xx(i+1, j) - T_xx(i, j))/dx + (T_xy(i, j) - T_xy(i, j-1))/dy )
    visc_u = zero
    do j = 1, ny
      do i = 0, nx
        visc_u(i, j) = - ( (T_xx(i+1, j) - T_xx(i, j)) * inv_dx + &
                           (T_xy(i, j)   - T_xy(i, j-1)) * inv_dy )
      end do
    end do

    ! d) Divergence pour la composante v (sur faces v : j = 0 .. ny)
    !    - div(T)_v = - ( (T_yx(i, j) - T_yx(i-1, j))/dx + (T_yy(i, j+1) - T_yy(i, j))/dy )
    visc_v = zero
    do j = 0, ny
      do i = 1, nx
        visc_v(i, j) = - ( (T_yx(i, j)   - T_yx(i-1, j)) * inv_dx + &
                           (T_yy(i, j+1) - T_yy(i, j))   * inv_dy )
      end do
    end do

  end subroutine op_viscous_stress

  ! ----------------------------------------------------------------------------
  ! 4. Operateur Convectif discret : div(rho * u_adv (x) u)
  !    Schema centre d'ordre 2 sur maillage MAC
  ! ----------------------------------------------------------------------------
  subroutine op_advection(u_adv, v_adv, rho_c, rho_u, rho_v, u_vel, v_vel, &
                          adv_u, adv_v, nx, ny, nh)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)  :: u_adv, u_vel
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)  :: v_adv, v_vel
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)  :: rho_c
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)  :: rho_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)  :: rho_v
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(out) :: adv_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(out) :: adv_v

    integer(ip) :: i, j
    real(dp) :: inv_dx, inv_dy
    real(dp) :: mdot_E, mdot_W, mdot_N, mdot_S
    real(dp) :: u_E, u_W, u_N, u_S
    real(dp) :: v_E, v_W, v_N, v_S

    inv_dx = one / dx_grid
    inv_dy = one / dy_grid

    adv_u = zero
    adv_v = zero

    ! Advection de la composante u (sur volume de contrôle centré sur la face u)
    do j = 1, ny
      do i = 1, nx-1
        ! Flux Est au centre de la cellule (i+1, j)
        mdot_E = rho_c(i+1, j) * half * (u_adv(i, j) + u_adv(i+1, j))
        u_E    = half * (u_vel(i, j) + u_vel(i+1, j))

        ! Flux Ouest au centre de la cellule (i, j)
        mdot_W = rho_c(i, j) * half * (u_adv(i-1, j) + u_adv(i, j))
        u_W    = half * (u_vel(i-1, j) + u_vel(i, j))

        ! Flux Nord au sommet (i, j)
        mdot_N = half * (rho_v(i, j) + rho_v(i+1, j)) * half * (v_adv(i, j) + v_adv(i+1, j))
        u_N    = half * (u_vel(i, j) + u_vel(i, j+1))

        ! Flux Sud au sommet (i, j-1)
        mdot_S = half * (rho_v(i, j-1) + rho_v(i+1, j-1)) * half * (v_adv(i, j-1) + v_adv(i+1, j-1))
        u_S    = half * (u_vel(i, j-1) + u_vel(i, j))

        adv_u(i, j) = (mdot_E * u_E - mdot_W * u_W) * inv_dx + &
                      (mdot_N * u_N - mdot_S * u_S) * inv_dy
      end do
    end do

    ! Advection de la composante v (sur volume de contrôle centré sur la face v)
    do j = 1, ny-1
      do i = 1, nx
        ! Flux Est au sommet (i, j)
        mdot_E = half * (rho_u(i, j) + rho_u(i, j+1)) * half * (u_adv(i, j) + u_adv(i, j+1))
        v_E    = half * (v_vel(i, j) + v_vel(i+1, j))

        ! Flux Ouest au sommet (i-1, j)
        mdot_W = half * (rho_u(i-1, j) + rho_u(i-1, j+1)) * half * (u_adv(i-1, j) + u_adv(i-1, j+1))
        v_W    = half * (v_vel(i-1, j) + v_vel(i, j))

        ! Flux Nord au centre de la cellule (i, j+1)
        mdot_N = rho_c(i, j+1) * half * (v_adv(i, j) + v_adv(i, j+1))
        v_N    = half * (v_vel(i, j) + v_vel(i, j+1))

        ! Flux Sud au centre de la cellule (i, j)
        mdot_S = rho_c(i, j) * half * (v_adv(i, j-1) + v_adv(i, j))
        v_S    = half * (v_vel(i, j-1) + v_vel(i, j))

        adv_v(i, j) = (mdot_E * v_E - mdot_W * v_W) * inv_dx + &
                      (mdot_N * v_N - mdot_S * v_S) * inv_dy
      end do
    end do

  end subroutine op_advection

  ! ----------------------------------------------------------------------------
  ! 5. Application des Conditions aux Limites Physiques (No-slip, Cavité, Périodique)
  ! ----------------------------------------------------------------------------
  subroutine apply_velocity_bcs(u_vel, v_vel, nx, ny, nh, u_top_lid)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: u_vel
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(inout) :: v_vel
    real(dp), optional, intent(in) :: u_top_lid

    integer(ip) :: i, j
    real(dp) :: u_lid

    u_lid = zero
    if (present(u_top_lid)) u_lid = u_top_lid

    ! 1. Si domaine non-périodique en X : parois Ouest et Est
    if (.not. periods(1)) then
      ! Paroi Ouest (x = 0) sur le premier processus de la ligne
      if (coords(1) == 0) then
        do j = 1-nh, ny+nh
          u_vel(0, j) = zero
          ! Extrapolation antissymetrique pour vitesse tangentielle v
          do i = 1, nh
            u_vel(-i, j) = - u_vel(i, j)
          end do
        end do
        do j = 0-nh, ny+nh
          do i = 1, nh
            v_vel(1-i, j) = - v_vel(i, j)
          end do
        end do
      end if

      ! Paroi Est (x = Lx) sur le dernier processus de la ligne
      if (coords(1) == dims(1) - 1) then
        do j = 1-nh, ny+nh
          u_vel(nx, j) = zero
          do i = 1, nh
            u_vel(nx+i, j) = - u_vel(nx-i, j)
          end do
        end do
        do j = 0-nh, ny+nh
          do i = 1, nh
            v_vel(nx+i, j) = - v_vel(nx-i+1, j)
          end do
        end do
      end if
    end if

    ! 2. Si domaine non-périodique en Y : parois Sud et Nord
    if (.not. periods(2)) then
      ! Paroi Sud (y = 0) sur le premier processus de la colonne
      if (coords(2) == 0) then
        do i = 1-nh, nx+nh
          v_vel(i, 0) = zero
          do j = 1, nh
            v_vel(i, -j) = - v_vel(i, j)
          end do
        end do
        do i = 0-nh, nx+nh
          do j = 1, nh
            u_vel(i, 1-j) = - u_vel(i, j)
          end do
        end do
      end if

      ! Paroi Nord (y = Ly) sur le dernier processus de la colonne (avec paroi entraînée éventuelle)
      if (coords(2) == dims(2) - 1) then
        do i = 1-nh, nx+nh
          v_vel(i, ny) = zero
          do j = 1, nh
            v_vel(i, ny+j) = - v_vel(i, ny-j)
          end do
        end do
        do i = 0-nh, nx+nh
          do j = 1, nh
            u_vel(i, ny+j) = two * u_lid - u_vel(i, ny-j+1)
          end do
        end do
      end if
    end if

  end subroutine apply_velocity_bcs

  ! ----------------------------------------------------------------------------
  ! 6. Produit Matrice-Vecteur global pour le systeme selle complet A * X = Y
  !    [ F_uu   F_uv   B_pu^T ] [ X_u ]   [ Y_u ]
  !    [ F_vu   F_vv   B_pv^T ] [ X_v ] = [ Y_v ]
  !    [ B_u    B_v      0    ] [ X_p ]   [ Y_p ]
  ! ----------------------------------------------------------------------------
  subroutine matvec_coupled(X_u, X_v, X_p, Y_u, Y_v, Y_p, &
                            u_adv, v_adv, rho_c, rho_u, rho_v, &
                            kappa, zeta, eta, dt_step, nx, ny, nh, include_advection)
    integer(ip), intent(in) :: nx, ny, nh
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: X_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(inout) :: X_v
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(inout) :: X_p
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(out)   :: Y_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(out)   :: Y_v
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(out)   :: Y_p
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: u_adv
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)    :: v_adv
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: rho_c
    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: rho_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh), intent(in)    :: rho_v
    real(dp), dimension(1-nh:nx+nh, 1-nh:ny+nh), intent(in)    :: kappa
    real(dp), dimension(0-nh:nx+nh, 0-nh:ny+nh), intent(in)    :: zeta, eta
    real(dp), intent(in) :: dt_step
    logical, optional, intent(in) :: include_advection

    real(dp), dimension(0-nh:nx+nh, 1-nh:ny+nh) :: visc_u, adv_u, grad_p_u
    real(dp), dimension(1-nh:nx+nh, 0-nh:ny+nh) :: visc_v, adv_v, grad_p_v
    logical :: with_adv
    integer(ip) :: i, j

    with_adv = .true.
    if (present(include_advection)) with_adv = include_advection

    ! a) Échanges de halos pour assurer la cohérence des cellules fantômes
    call mpi_exchange_halos_u(X_u, nx, ny, nh)
    call mpi_exchange_halos_v(X_v, nx, ny, nh)
    call mpi_exchange_halos_scalar(X_p, nx, ny, nh)

    ! b) Calcul du tenseur visqueux - div(T)
    call op_viscous_stress(X_u, X_v, kappa, zeta, eta, visc_u, visc_v, nx, ny, nh)

    ! c) Calcul du gradient de pression B_p^T * X_p
    call op_grad_pressure(X_p, grad_p_u, grad_p_v, nx, ny, nh)

    ! d) Calcul convectif si active
    if (with_adv) then
      call op_advection(u_adv, v_adv, rho_c, rho_u, rho_v, X_u, X_v, adv_u, adv_v, nx, ny, nh)
    else
      adv_u = zero
      adv_v = zero
    end if

    ! e) Assemblage du bloc de quantite de mouvement en u
    !    Y_u = (rho_u / dt) * X_u + Adv_u - div(T)_u + Grad_x(X_p)
    Y_u = zero
    do j = 1, ny
      do i = 1, nx-1
        Y_u(i, j) = (rho_u(i, j) / dt_step) * X_u(i, j) + adv_u(i, j) + &
                    visc_u(i, j) + grad_p_u(i, j)
      end do
      if (periods(1)) then
        Y_u(0, j)  = (rho_u(0, j) / dt_step) * X_u(0, j) + adv_u(0, j) + visc_u(0, j) + grad_p_u(0, j)
        Y_u(nx, j) = (rho_u(nx, j) / dt_step) * X_u(nx, j) + adv_u(nx, j) + visc_u(nx, j) + grad_p_u(nx, j)
      else
        if (coords(1) == 0) then
          Y_u(0, j) = X_u(0, j)
        else
          Y_u(0, j) = (rho_u(0, j) / dt_step) * X_u(0, j) + adv_u(0, j) + visc_u(0, j) + grad_p_u(0, j)
        end if
        if (coords(1) == dims(1) - 1) then
          Y_u(nx, j) = X_u(nx, j)
        else
          Y_u(nx, j) = (rho_u(nx, j) / dt_step) * X_u(nx, j) + adv_u(nx, j) + visc_u(nx, j) + grad_p_u(nx, j)
        end if
      end if
    end do

    ! f) Assemblage du bloc de quantite de mouvement en v
    !    Y_v = (rho_v / dt) * X_v + Adv_v - div(T)_v + Grad_y(X_p)
    Y_v = zero
    do j = 1, ny-1
      do i = 1, nx
        Y_v(i, j) = (rho_v(i, j) / dt_step) * X_v(i, j) + adv_v(i, j) + &
                    visc_v(i, j) + grad_p_v(i, j)
      end do
    end do
    do i = 1, nx
      if (periods(2)) then
        Y_v(i, 0)  = (rho_v(i, 0) / dt_step) * X_v(i, 0) + adv_v(i, 0) + visc_v(i, 0) + grad_p_v(i, 0)
        Y_v(i, ny) = (rho_v(i, ny) / dt_step) * X_v(i, ny) + adv_v(i, ny) + visc_v(i, ny) + grad_p_v(i, ny)
      else
        if (coords(2) == 0) then
          Y_v(i, 0) = X_v(i, 0)
        else
          Y_v(i, 0) = (rho_v(i, 0) / dt_step) * X_v(i, 0) + adv_v(i, 0) + visc_v(i, 0) + grad_p_v(i, 0)
        end if
        if (coords(2) == dims(2) - 1) then
          Y_v(i, ny) = X_v(i, ny)
        else
          Y_v(i, ny) = (rho_v(i, ny) / dt_step) * X_v(i, ny) + adv_v(i, ny) + visc_v(i, ny) + grad_p_v(i, ny)
        end if
      end if
    end do

    ! g) Assemblage du bloc d'incompressibilite (Divergence B_u X_u + B_v X_v)
    call op_divergence(X_u, X_v, Y_p, nx, ny, nh)

  end subroutine matvec_coupled

end module mod_operators
