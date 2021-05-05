
all: main

build: LEGACY.o

main: main.o build
	g++ $(LDFLAGS) -o main.exe  main.o LEGACY.o -lntl -lgmp -lm

main.o:main.cpp

LEGACY.o: LEGACY.cpp LEGACY.h

clean:
	rm *.o
	rm *.txt
	rm *.exe