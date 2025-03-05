# compiler and flags
CXX = g++
CXXFLAGS = -fdiagnostics-color=always\
		   -Wall\
		   -g

TARGET = main.exe
CONFIG_MAIN  = main
CONFIG_DEBUG = debug

SRCS = main.cpp
OBJS = $(SRCS:.cpp=.o)


all: $(TARGET) run

build: $(TARGET)

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
