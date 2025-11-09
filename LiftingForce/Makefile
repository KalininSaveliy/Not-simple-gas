# RM
ifeq ($(OS),Windows_NT)
	RM = del /Q
	RMDIR = rmdir /S /Q
else
	RM = rm -f
	RMDIR = rm -rf
endif

# Python
PYTHON = .\.venv\Scripts\python.exe

# compiler and flags
CXX = g++
CXXFLAGS = -fdiagnostics-color=always\
		   -Wall\
		   -g\
		   -fopenmp\
		   -std=c++23

TARGET = main.exe

CONFIG_FOLDER = config/
CONFIG_MAIN = $(CONFIG_FOLDER)/main.json
CONFIG_CPP_MAIN  = $(CONFIG_FOLDER)main_Cpp.txt
CONFIG_CPP_DEBUG = $(CONFIG_FOLDER)debug_Cpp.txt

SRCS = main.cpp
OBJS = $(SRCS:.cpp=.o)


all: $(TARGET) run

gif:
	$(PYTHON) .\make_graph\gif_flat.py

build: $(TARGET) $(CONFIG_CPP_MAIN) ;

$(CONFIG_CPP_MAIN): $(CONFIG_MAIN)
	$(PYTHON) .\config\cfg.py

run: $(TARGET)
	$(TARGET) $(CONFIG_CPP_MAIN)
	$(MAKE) gif

rundebug: $(TARGET)
	$(TARGET) $(CONFIG_CPP_DEBUG)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# .PHONY: clean
clean:
	$(RM) $(OBJS)
	$(RM) $(TARGET)

cleandata:
	$(RM) data\main\*
	$(RM) data/debug/*