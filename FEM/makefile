# Project Name
TARGET = main

# Compiler
CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -fopenmp

# CUDA Support (Uncomment if using CUDA)
# NVCC = nvcc
# CUDAFLAGS = -arch=sm_60 -O2

# Directories
SRC_DIR = src
INCLUDE_DIR = $(PWD)/include
BUILD_DIR = build

# Find all C++ source files
SRCS := $(wildcard $(SRC_DIR)/**/*.cpp $(SRC_DIR)/*.cpp main.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS))
OBJS := $(patsubst main.cpp, $(BUILD_DIR)/main.o, $(OBJS))

# Include Directories
INCLUDES = -I$(INCLUDE_DIR) -I/usr/include/openblas

# Libraries (BLAS, Eigen, etc.)
LIBS = -L/usr/lib64 -lopenblas -lcblas -llapack


# Default rule (build the project)
all: $(TARGET)

# Compile the project
$(TARGET): $(OBJS)
	@echo "Building $(TARGET)..."
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET) $(LIBS)

# Build object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Compile the main file
$(BUILD_DIR)/main.o: main.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c main.cpp -o $(BUILD_DIR)/main.o

# Run the program
run: $(TARGET)
	./$(TARGET)

# Debug mode
debug: CXXFLAGS += -g
debug: clean all

# Clean build files
clean:
	rm -rf $(BUILD_DIR) $(TARGET)
	@echo "Cleaned build files."

# Phony targets
.PHONY: all clean debug run