module mod_params
  use mod_precision, only: dp, ip, zero, one, half
  implicit none
  public

  ! =========================================================================
  ! 1. Paramètres géométriques du domaine physique [0, Lx] x [0, Ly]
  ! =========================================================================
  real(dp) :: Lx = one
  real(dp) :: Ly = one

  ! Nombre global de cellules scalaires (pression, couleur C)
  integer(ip) :: Nx_global = 64
  integer(ip) :: Ny_global = 64

  ! =========================================================================
  ! 2. Propriétés physiques des deux fluides (Modèle One-Fluid)
  ! =========================================================================
  real(dp) :: rho1 = 1000.0_dp   ! kg/m^3 (Fluide 1 - porteur / lourd)
  real(dp) :: rho2 = 1.0_dp      ! kg/m^3 (Fluide 2 - dispersé / léger)
  real(dp) :: mu1  = 1.0_dp      ! Pa.s
  real(dp) :: mu2  = 0.01_dp     ! Pa.s
  real(dp) :: sigma = zero       ! Tension de surface (N/m)
  real(dp) :: gx    = zero       ! Accélération gravité x (m/s^2)
  real(dp) :: gy    = -9.81_dp   ! Accélération gravité y (m/s^2)

  ! =========================================================================
  ! 3. Contrôle temporel et stabilité
  ! =========================================================================
  real(dp) :: dt = 1.0e-3_dp
  real(dp) :: t_final = one
  real(dp) :: cfl_target = half
  integer(ip) :: max_time_steps = 10000

  ! Schéma temporel : 1 = Euler ordre 1, 2 = Gear / BDF2 ordre 2
  integer(ip) :: time_order = 1

  ! =========================================================================
  ! 4. Paramètres du solveur Krylov BiCGStab(2) et Préconditionneur PCD
  ! =========================================================================
  real(dp)    :: tol_bicgstab2     = 1.0e-6_dp
  integer(ip) :: max_iter_bicgstab2 = 500
  real(dp)    :: tol_subsolver     = 1.0e-7_dp
  integer(ip) :: max_iter_subsolver = 100

  ! =========================================================================
  ! 5. Paramètres de lissage SVOF (Helmholtz) et Pénalisation frontière
  ! =========================================================================
  real(dp) :: svof_diff_coeff = half   ! Paramètre Li pour D = Li * Delta_h^2
  real(dp) :: alpha_penalty   = 1.0e12_dp ! Facteur pénalisation de Darcy

end module mod_params
