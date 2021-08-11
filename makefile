
all: main 

build: LEGACY.o

main: main_trash.o build
	g++ $(LDFLAGS) -o main_trash.exe main_trash.o LEGACY.o -lntl -lgmp -lm			
	
main.o:main_trash.cpp 

LEGACY.o: LEGACY.cpp LEGACY.h

clean:
	rm *.exe
	rm *.o
	rm *.txt

