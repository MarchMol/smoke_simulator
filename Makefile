# Compiler
CC = gcc

# Include and library paths
INCLUDE = -Iinclude
LIBPATH = -Llib

# Libraries to link
LIBS = -lglfw3 -lopengl32 -lgdi32

# Output executable
OUT = bin/main.exe

# Source files
SRC = main.c

# Default target
all: compile

# Compile the program
compile:
	@echo "Compiling..."
	$(CC) $(SRC) $(INCLUDE) $(LIBPATH) $(LIBS) -o $(OUT)
	@echo "Done!"

# Run the program
run: compile
	@echo "Running..."
	$(OUT)