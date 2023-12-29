all: relax.exe

relax.exe: main.o KorobovGrid.o
	g++ main.o KorobovGrid.o -o relax.exe

KorobovGrid.o: .\source\KorobovGrid.cpp
	g++ -c .\source\KorobovGrid.cpp -o KorobovGrid.o

main.o: relaxation.cpp
	g++ -c relaxation.cpp -o main.o

clean:
	del *.o relax.exe
