# Python
PYTHON = .\.venv\Scripts\python.exe

# compiler and flags
CXX = g++
CXXFLAGS = -fdiagnostics-color=always\
		   -Wall\
		   -g\
		   -std=c++23

TARGET = main.exe

CONFIG_FOLDER = config/
CONFIG_MAIN  = $(CONFIG_FOLDER)main_Cpp.txt
CONFIG_DEBUG = $(CONFIG_FOLDER)debug_Cpp.txt

SRCS = main.cpp
OBJS = $(SRCS:.cpp=.o)


all: $(TARGET) run

gif:
	PYTHON .\make_graph\gif_flat.py

build: $(TARGET)
	PYTHON .\config\cfg.py

run: $(TARGET)
	$(TARGET) $(CONFIG_MAIN)

rundebug: $(TARGET)
	$(TARGET) $(CONFIG_DEBUG)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# .PHONY: clean
clean:
	del $(OBJS)
	del $(TARGET)

cleandata:
	del .\data\main\*
	del .\data\debug\*

# varr:
# 	gnuplot -c make_graph/gif_flat.gp @.\config\tmp.txt 
