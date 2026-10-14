# ==============================================================================
# Makefile - Solveur 2D Tout-Couplé Fortran & OpenMPI (FUGU-2D)
# ==============================================================================

# Répertoires
SRC_DIR  := src
TEST_DIR := tests
OBJ_DIR  := build
MOD_DIR  := mod
BIN_DIR  := bin

# Forcer l'utilisation de mpifort (ignore le f77 par défaut de GNU Make et l'environnement)
override FC     := mpifort
override FFLAGS := -O3 -Wall -Wextra -fcheck=all -g -J$(MOD_DIR) -I$(MOD_DIR)

# Modules principaux
SRCS     := $(SRC_DIR)/mod_precision.f90 \
            $(SRC_DIR)/mod_params.f90 \
            $(SRC_DIR)/mod_mpi.f90

OBJS     := $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRCS))

# Cibles principales
all: directories $(BIN_DIR)/test_halo_exchange

directories:
	@mkdir -p $(OBJ_DIR) $(MOD_DIR) $(BIN_DIR)

# Compilation des modules
$(OBJ_DIR)/mod_precision.o: $(SRC_DIR)/mod_precision.f90 | directories
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/mod_params.o: $(SRC_DIR)/mod_params.f90 $(OBJ_DIR)/mod_precision.o | directories
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/mod_mpi.o: $(SRC_DIR)/mod_mpi.f90 $(OBJ_DIR)/mod_precision.o $(OBJ_DIR)/mod_params.o | directories
	$(FC) $(FFLAGS) -c $< -o $@

# Exécutable de test
$(BIN_DIR)/test_halo_exchange: $(TEST_DIR)/test_halo_exchange.f90 $(OBJS) | directories
	$(FC) $(FFLAGS) $< $(OBJS) -o $@

# Exécution du test unitaire
test_mpi: $(BIN_DIR)/test_halo_exchange
	@echo "Lancement du test d'échange de halos MPI sur 4 cœurs..."
	mpirun -np 4 ./$(BIN_DIR)/test_halo_exchange

clean:
	rm -rf $(OBJ_DIR) $(MOD_DIR) $(BIN_DIR)

.PHONY: all directories test_mpi clean
