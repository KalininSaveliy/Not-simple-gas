# OS commands (for win use cmd)
ifeq ($(OS),Windows_NT)
	EXE = .exe
	RM = del /Q
	RMDIR = rmdir /S /Q
	SRCS = $(shell cmd /C "for /R $(SRC_DIR) %%f in (*.cpp) do @echo %%f")
	MKDIR = if not exist "$(subst /,\,$(dir $@))" mkdir "$(subst /,\,$(dir $@))"
	PATHSEP = \\
else
	EXE =
	RM = rm -f
	RMDIR = rm -rf
	FIND_CPP = $(shell find $(SRC_DIR) -type f -name "*.cpp")
	MKDIR = mkdir -p $(dir $@)
	PATHSEP = /
endif


# Python
PYTHON = .\.venv\Scripts\python.exe

# Compiler and Flags
CXX = g++
CXXFLAGS = -fdiagnostics-color=always\
		   -Wall\
		   -g\
		   -fopenmp\
		   -O2\
		   -Wextra\
		   -std=c++23

# Pathes
SRC_DIR = source
BUILD_DIR = build
TARGET = $(BUILD_DIR)/main$(EXE)

SRCS = $(wildcard $(SRC_DIR)/*/*.cpp) $(wildcard $(SRC_DIR)/*/*/*.cpp)
OBJS = $(SRCS:$(SRC_DIR)%.cpp=$(BUILD_DIR)%.o)

INC_DIRS = $(wildcard $(SRC_DIR)/*/*[!.h][!.cpp])
INC_FLAGS = $(addprefix -I , $(INC_DIRS))

# Config
CONFIG_FOLDER = config/
CONFIG_MAIN = $(CONFIG_FOLDER)/main.json
CONFIG_CPP_MAIN  = $(CONFIG_FOLDER)main_Cpp.txt
CONFIG_CPP_DEBUG = $(CONFIG_FOLDER)debug_Cpp.txt


all: $(TARGET) run

gif:
	$(PYTHON) .\make_graph\gif_flat.py


$(CONFIG_CPP_MAIN): $(CONFIG_MAIN)
	$(PYTHON) .\config\cfg.py

run: $(TARGET)
	$(TARGET) $(CONFIG_CPP_MAIN)
	$(MAKE) gif

rundebug: $(TARGET)
	$(TARGET) $(CONFIG_CPP_DEBUG)

info: $(OBJS) $(SRCS)
	@echo $(OBJS)
	@echo $(SRCS)
	@echo $(INC_DIRS)

build: $(TARGET) $(CONFIG_CPP_MAIN) ;

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(INC_FLAGS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@$(MKDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INC_FLAGS)

clean:
	-$(RMDIR) $(BUILD_DIR)

cleandata:
	$(RM) data\main\*
	$(RM) data/debug/*
