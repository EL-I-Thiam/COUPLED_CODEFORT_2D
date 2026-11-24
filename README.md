# Solveur 2D Tout-Couplé Volumes Finis (One-Fluid + VOF) en Fortran & OpenMPI

Ce projet implémente un solveur numérique bidimensionnel haute performance en volumes finis sur maillage cartésien décalé (MAC), dédié à la simulation directe des écoulements diphasiques incompressibles à forts contrastes de masse volumique et de viscosité.

La méthodologie est directement issue des travaux du laboratoire MSME (Université Gustave Eiffel) :
* **Thèse de doctorat (2022) :** *Développement d'un solveur tout-couplé parallèle 3D pour la simulation des écoulements diphasiques incompressibles à forts rapports de viscosités et de masses volumiques* — Mohamed El Ouafa (Dir. S. Vincent, V. Le Chenadec).
* **Article (Computers & Fluids, 2023) :** *Fully-coupled parallel solver for the simulation of two-phase incompressible flows* — S. El Ouafa, S. Vincent, V. Le Chenadec, B. Trouette.

---

## 1. Caractéristiques numériques du solveur

1. **Modèle diphasique :**
   - Modèle *One-Fluid* (MUF) de Kataoka.
   - Suivi d'interface par Volume of Fluid (VOF) conservateur (Weymouth & Yue, 2010).
   - Lissage régularisé par équation de Helmholtz (SVOF) pour le calcul sans bruit de la normale $\mathbf{n}$, de la courbure $\kappa$ et de la force capillaire CSF.
   - Moyenne harmonique des viscosités aux sommets des mailles pour une convergence d'ordre 2.

2. **Couplage vitesse-pression monolithique (Tout-couplé 2D) :**
   - Résolution simultanée de $(\mathbf{u}, p)$ sans fractionnement d'opérateurs ni couches limites artificielles de pression.
   - Discrétisation sur grille décalée MAC à 11 diagonales non nulles en 2D.
   - Conditions aux limites imposées par méthode de pénalisation (Tableau 3.2 de la thèse).

3. **Algèbre linéaire & Parallélisme HPC :**
   - Solveur de Krylov : **BiCGStab(2)** préconditionné à gauche.
   - Préconditionneur par blocs $\mathcal{P}$ : bloc Gauss-Seidel d'ordre 0 pour la vitesse et PCD (*Pressure Convection Diffusion*) pour le complément de Schur de pression.
   - Parallélisation massive par décomposition de domaine 2D avec **OpenMPI**.

---

## 2. Compilation et exécution

### Prérequis
* **Compilateur Fortran :** `gfortran` (version $\ge 11$) ou `ifort`/`ifx`.
* **Bibliothèque MPI :** `OpenMPI` ou `MPICH`.
* **Outils :** `make`.

### Installation des outils

* **Sur macOS (via Homebrew) :**
  ```bash
  brew install gcc open-mpi cmake
  ```

* **Sur Linux (Ubuntu/Debian) :**
  ```bash
  sudo apt update && sudo apt install -y gfortran libopenmpi-dev openmpi-bin cmake make git
  ```

* **Sur Linux (Rocky/CentOS/Fedora) :**
  ```bash
  sudo dnf install -y gcc-gfortran openmpi openmpi-devel cmake make git
  ```

### Commandes Make et exécution
* Compiler le projet (`fugu2d` et les 8 tests unitaires/benchmarks) :
  ```bash
  make
  ```
* Exécuter le solveur configurable (`input.nml`) sur 4 cœurs MPI :
  ```bash
  mpirun -np 4 ./bin/fugu2d
  ```
* Lancer l'ensemble des 8 tests de validation sur 4 processus :
  ```bash
  make tests
  ```
* Lancer un test spécifique (`test_mpi`, `test_operators`, `test_cavity`, `test_droplet`, `test_laplace`, `test_bubble`, `test_poiseuille`, `test_week6`) :
  ```bash
  make test_week6
  ```
* Nettoyer les fichiers objets et modules :
  ```bash
  make clean
  ```

---

## 3. Feuille de route de développement hebdomadaire

* **Semaine 1 (Terminée) :** Socle MPI 2D, topologie cartésienne `MPI_Cart_create` et routines d'échange de halos pour scalaires, vitesses décalées et sommets.
* **Semaine 2 (Terminée) :** Opérateurs discrets de quantité de mouvement, tenseur visqueux à viscosités variables, divergence et gradient de pression.
* **Semaine 3 (Terminée) :** Solveur Krylov BiCGStab / BiCGStab(2) distribué et préconditionneur par blocs PCD (Validation : Cavité entraînée de Botella & Peyret).
* **Semaine 4 (Terminée) :** Transport VOF-TVD (Superbee), lissage SVOF Helmholtz et force capillaire CSF (Validation : Goutte statique / loi de Laplace).
* **Semaine 5 (Terminée) :** Couplage diphasique complet NS + VOF + CSF et sorties VTK (`.vtk`) pour ParaView (Validation : Goutte statique couplée et bulle ascendante de Hysing Cas 1).
* **Semaine 6 (Terminée) :** Schéma temporel BDF2 (`time_order = 2`), transport consistant masse-quantité de mouvement (AMP, Hysing Cas 2 à $\rho_1/\rho_2 = 10^3$), pénalisation de Darcy ($\alpha = 10^{12}$), Poiseuille diphasique analytique, pilote `fugu2d` (`input.nml`) et reprise binaire MPI (*checkpoint/restart*).
