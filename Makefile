CXX = g++
CXXFLAGS = -std=c++17 -O3
INCLUDES = -Iinclude -I.
LDFLAGS = -lstdc++fs

# Directories
SRC_DIR = src
OBJ_DIR = build/obj
BIN_DIR = bin

# Source files
CORE_SRC = $(wildcard $(SRC_DIR)/core/*.cpp)
CROSSBAR_SRC = $(wildcard $(SRC_DIR)/crossbar_model/*.cpp)
MEMRISTOR_SRC = $(wildcard $(SRC_DIR)/memristor_model/*.cpp)
SRC_FILES = $(CORE_SRC) $(CROSSBAR_SRC) $(MEMRISTOR_SRC)

# Object files with path
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

# Target executable
TARGET = $(BIN_DIR)/xbar_simulator

all: create_dirs $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS)

# Create necessary directories for object files
create_dirs:
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(OBJ_DIR)/core
	@mkdir -p $(OBJ_DIR)/crossbar_model
	@mkdir -p $(OBJ_DIR)/memristor_model

# Pattern rule for object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR)
	rm -f $(TARGET)

.PHONY: all clean create_dirs
