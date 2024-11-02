# Makefile for compiling the Peeling C program with architecture-specific optimizations

# ===============================
# Variables
# ===============================

# Detect the operating system and architecture
OS := $(shell uname)

# Set the compiler based on the OS
ifeq ($(OS), Darwin)
    CC := clang
else
    CC := gcc
endif

# Compiler flags
CFLAGS := -O3 -march=native -mtune=native -Wall -Wextra -std=c11

# Library flags
LIBS := -lm

# Handle OpenMP on macOS
ifeq ($(OS), Darwin)
    # Detect Homebrew installation prefix for libomp
    HOMEBREW_PREFIX := $(shell brew --prefix libomp 2>/dev/null)

    # Default to /usr/local if libomp is not found via brew
    ifeq ($(HOMEBREW_PREFIX),)
        HOMEBREW_PREFIX := /usr/local
    endif

    # Add OpenMP include and library paths
    CFLAGS += -I$(HOMEBREW_PREFIX)/include
    LDFLAGS += -L$(HOMEBREW_PREFIX)/lib -lomp

    # Optional: Specify the path to the OpenMP library if installed in a non-standard location
    # Uncomment and modify the following lines if necessary
    # CFLAGS += -I/path/to/libomp/include
    # LDFLAGS += -L/path/to/libomp/lib -lomp
else
	# Linker flags
	LDFLAGS := -fopenmp
endif


# Source and target
SRC := peeling_omp.c
TARGET := Peeling_omp

# ===============================
# Targets
# ===============================

.PHONY: all clean

# Default target
all: $(TARGET)

# Link the executable
$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS) $(LIBS)

# Clean up build artifacts
clean:
	@echo "Cleaning up..."
	@rm -f $(TARGET)
