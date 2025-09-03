# Compiler
CC = gcc

# Include and library paths
INCLUDE = -Iinclude
LIBPATH = -Llib

# Libraries to link
LIBS = -lglfw3 -lopengl32 -lgdi32 -lglu32

# Output executable
OUT_SEQ = bin/main_seq.exe
OUT_OMP = bin/main_omp.exe

# Flags
CFLAGS      = -Ofast -march=native -ffast-math
CFLAGS_OMP  = -Ofast -march=native -ffast-math -fopenmp -fopenmp-simd

# Source files (explícitos para no mezclar mains) ----
MAIN_SEQ   = src/main.c
MAIN_PAR   = src/main_parallel.c
SMOKE_SEQ  = src/smoke.c
SMOKE_OMP  = src/smoke_parallel.c
STATE  = src/state.c
GLAD = src/glad.c

SRC_SEQ = $(MAIN_SEQ) $(STATE) $(SMOKE_SEQ)
SRC_OMP = $(GLAD) $(MAIN_PAR) $(STATE) $(SMOKE_OMP)

# Utils
BIN_DIR = bin

.PHONY: all seq run-seq omp run-omp clean

all: seq omp

# Asegura bin/
$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

# --------- Secuencial ----------
$(OUT_SEQ): $(SRC_SEQ) | $(BIN_DIR)
	@echo "Compiling sequential version..."
	$(CC) $(CFLAGS) $(SRC_SEQ) $(INCLUDE) $(LIBPATH) $(LIBS) -o $(OUT_SEQ)

seq: $(OUT_SEQ)

run-seq: $(OUT_SEQ)
	@echo "Running sequential..."
	$(OUT_SEQ)

# --------- Paralelo (OpenMP) ----------
$(OUT_OMP): $(SRC_OMP) | $(BIN_DIR)
	@echo "Compiling parallel version with OpenMP..."
	$(CC) $(CFLAGS_OMP) $(SRC_OMP) $(INCLUDE) $(LIBPATH) $(LIBS) -o $(OUT_OMP)

omp: $(OUT_OMP)

run-omp: $(OUT_OMP)
	@echo "Running parallel..."
	$(OUT_OMP)

# --------- Limpieza ----------
clean:
	@echo "Cleaning..."
	@rm -f $(OUT_SEQ) $(OUT_OMP)