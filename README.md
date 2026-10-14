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

### Commandes Make
* Compiler le projet et les tests unitaires :
  ```bash
  make
  ```
* Lancer le test de communication MPI cartésienne (sur 4 processus) :
  ```bash
  make test_mpi
  ```
* Nettoyer les fichiers objets et modules :
  ```bash
  make clean
  ```

---

## 3. Feuille de route de développement hebdomadaire

* **Semaine 1 (En cours) :** Socle MPI 2D, topologie cartésienne `MPI_Cart_create` et routines d'échange de halos pour scalaires, vitesses décalées et sommets.
* **Semaine 2 :** Opérateurs discrets de quantité de mouvement, tenseur visqueux à viscosités variables, divergence et gradient de pression.
* **Semaine 3 :** Solveur Krylov BiCGStab(2) distribué et préconditionneur par blocs PCD (Validation : Cavité entraînée de Botella & Peyret).
* **Semaine 4 :** Transport VOF de Weymouth & Yue, lissage SVOF Helmholtz et force capillaire CSF (Validation : Goutte statique / loi de Laplace).
* **Semaine 5 :** Couplage diphasique complet et schéma de transport consistant AMP (Validation : Poiseuille diphasique analytique).
* **Semaine 6 :** Cas-tests raides (Bulle de Hysing, chute de cylindre dense à $\rho_1/\rho_2 = 10^5, \mu_1/\mu_2 = 10^{10}$) et sorties VTK (.vtr) pour ParaView.
