CC = g++  # Компилятор
CFLAGS = -c  # флаги компиляции
# LDFLAGS=smth  флаги линковки, компановки и связывания
SOURCES = .\source\KorobovGrid.cpp main1.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = main1.exe


all: $(EXECUTABLE)

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

main1.o: .\source\GasRelaxation.h .\source\Iterator3dim.h
KorobovGrid.o: .\source\KorobovGrid.h
# main.exe: main.o KorobovGrid.o
# 	g++ main.o KorobovGrid.o -o main.exe

# main.o: main.cpp
# 	g++ -c main.cpp -o main.o

# KorobovGrid.o: KorobovGrid.cpp
# 	g++ -c KorobovGrid.cpp -o KorobovGrid.o

clean:
	del $(OBJECTS) $(EXECUTABLE)
