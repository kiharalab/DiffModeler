# Detect the operating system and architecture
OS := $(shell uname)

# Set the compiler based on the OS
ifeq ($(OS), Darwin)
    CC := clang++
else
    CC := g++
endif

# Compiler flags
CFLAGS := -O3 -march=native -mtune=native -Wall -Wextra -std=c++0x

# Handle OpenMP on macOS
ifeq ($(OS),Darwin)
    CFLAGS  = -std=c++0x -Wall -Wextra -O3
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

SRC1    = src/training.cpp
SRC2    = src/scoring.cpp
SRC3    = src/scoring_omp.cpp
EXE1    = training
EXE2    = scoring
EXE3    = scoring_omp

all: $(EXE1) $(EXE2) $(EXE3)

$(EXE1): $(SRC1)
	$(CC) $(CFLAGS) -o $(EXE1) $(SRC1)

$(EXE2): $(SRC2)
	$(CC) $(CFLAGS) -o $(EXE2) $(SRC2)

$(EXE3): $(SRC3)
	$(CC) $(CFLAGS) -o $(EXE3) $(SRC3) $(LDFLAGS)
