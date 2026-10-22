# ==============================================================================
# Makefile - Solveur 2D Tout-Couplé Fortran & OpenMPI (FUGU-2D)
# Compatible macOS (Apple Silicon / Homebrew) et Linux (Ubuntu / Debian)
# ==============================================================================

# Répertoires
SRC_DIR  := src
TEST_DIR := tests
OBJ_DIR  := build
MOD_DIR  := mod
BIN_DIR  := bin

# Détection automatique de mpifort et mpirun (gère /opt/homebrew/bin sur Mac et /usr/bin sur Ubuntu)
override FC     := $(shell command -v mpifort 2>/dev/null || echo /opt/homebrew/bin/mpifort)
override MPIRUN := $(shell command -v mpirun 2>/dev/null || echo /opt/homebrew/bin/mpirun)
override FFLAGS := -O3 -Wall -Wextra -fcheck=all -g -J$(MOD_DIR) -I$(MOD_DIR)

# Modules principaux
SRCS     := $(SRC_DIR)/mod_precision.f90 \
            $(SRC_DIR)/mod_params.f90 \
            $(SRC_DIR)/mod_mpi.f90 \
            $(SRC_DIR)/mod_viscosity.f90 \
            $(SRC_DIR)/mod_operators.f90 \
            $(SRC_DIR)/mod_subsolvers.f90 \
            $(SRC_DIR)/mod_pcd_precond.f90 \
            $(SRC_DIR)/mod_bicgstab2.f90

OBJS     := $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRCS))

# Cibles principales
all: directories $(BIN_DIR)/test_halo_exchange $(BIN_DIR)/test_operators $(BIN_DIR)/test_lid_driven_cavity

directories:
	@mkdir -p $(OBJ_DIR) $(MOD_DIR) $(BIN_DIR)

# Compilation des modules (gestion rigoureuse des dépendances)
$(OBJ_DIR)/mod_precision.o: $(SRC_DIR)/mod_precision.f90 | directories
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/mod_params.o: $(SRC_DIR)/mod_params.f90 $(OBJ_DIR)/mod_precision.o | directories
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/mod_mpi.o: $(SRC_DIR)/mod_mpi.f90 $(OBJ_DIR)/mod_precision.o $(OBJ_DIR)/mod_params.o | directories
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/mod_viscosity.o: $(SRC_DIR)/mod_viscosity.f90 $(OBJ_DIR)/mod_precision.o $(OBJ_DIR)/mod_params.o $(OBJ_DIR)/mod_mpi.o | directories
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/mod_operators.o: $(SRC_DIR)/mod_operators.f90 $(OBJ_DIR)/mod_precision.o $(OBJ_DIR)/mod_params.o $(OBJ_DIR)/mod_mpi.o $(OBJ_DIR)/mod_viscosity.o | directories
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/mod_subsolvers.o: $(SRC_DIR)/mod_subsolvers.f90 $(OBJ_DIR)/mod_precision.o $(OBJ_DIR)/mod_params.o $(OBJ_DIR)/mod_mpi.o | directories
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/mod_pcd_precond.o: $(SRC_DIR)/mod_pcd_precond.f90 $(OBJ_DIR)/mod_precision.o $(OBJ_DIR)/mod_params.o $(OBJ_DIR)/mod_mpi.o $(OBJ_DIR)/mod_subsolvers.o | directories
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/mod_bicgstab2.o: $(SRC_DIR)/mod_bicgstab2.f90 $(OBJ_DIR)/mod_precision.o $(OBJ_DIR)/mod_params.o $(OBJ_DIR)/mod_mpi.o $(OBJ_DIR)/mod_operators.o $(OBJ_DIR)/mod_pcd_precond.o | directories
	$(FC) $(FFLAGS) -c $< -o $@

# Exécutables de tests
$(BIN_DIR)/test_halo_exchange: $(TEST_DIR)/test_halo_exchange.f90 $(OBJ_DIR)/mod_precision.o $(OBJ_DIR)/mod_params.o $(OBJ_DIR)/mod_mpi.o | directories
	$(FC) $(FFLAGS) $< $(OBJ_DIR)/mod_precision.o $(OBJ_DIR)/mod_params.o $(OBJ_DIR)/mod_mpi.o -o $@

$(BIN_DIR)/test_operators: $(TEST_DIR)/test_operators.f90 $(OBJ_DIR)/mod_precision.o $(OBJ_DIR)/mod_params.o $(OBJ_DIR)/mod_mpi.o $(OBJ_DIR)/mod_viscosity.o $(OBJ_DIR)/mod_operators.o | directories
	$(FC) $(FFLAGS) $< $(OBJ_DIR)/mod_precision.o $(OBJ_DIR)/mod_params.o $(OBJ_DIR)/mod_mpi.o $(OBJ_DIR)/mod_viscosity.o $(OBJ_DIR)/mod_operators.o -o $@

$(BIN_DIR)/test_lid_driven_cavity: $(TEST_DIR)/test_lid_driven_cavity.f90 $(OBJS) | directories
	$(FC) $(FFLAGS) $< $(OBJS) -o $@

# Exécution des tests unitaires
test_mpi: $(BIN_DIR)/test_halo_exchange
	@echo "Lancement du test d'échange de halos MPI sur 4 cœurs..."
	$(MPIRUN) -np 4 ./$(BIN_DIR)/test_halo_exchange

test_operators: $(BIN_DIR)/test_operators
	@echo "Lancement du test de validation des opérateurs spatiaux sur 4 cœurs..."
	$(MPIRUN) -np 4 ./$(BIN_DIR)/test_operators

test_cavity: $(BIN_DIR)/test_lid_driven_cavity
	@echo "Lancement du benchmark cavité entraînée (solveur tout-couplé) sur 4 cœurs..."
	$(MPIRUN) -np 4 ./$(BIN_DIR)/test_lid_driven_cavity

tests: test_mpi test_operators test_cavity

clean:
	rm -rf $(OBJ_DIR) $(MOD_DIR) $(BIN_DIR)

.PHONY: all directories test_mpi test_operators test_cavity tests clean
