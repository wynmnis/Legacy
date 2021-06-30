
all: main 

build: LEGACY.o

main: main.o main_AE.o main_CF.o main_PFA_Rader.o build
	g++ $(LDFLAGS) -o main.exe  main.o  LEGACY.o -lntl -lgmp -lm
	g++ $(LDFLAGS) -o main_AE.exe main_AE.o LEGACY.o -lntl -lgmp -lm
	g++ $(LDFLAGS) -o main_CF.exe main_CF.o LEGACY.o -lntl -lgmp -lm
	g++ $(LDFLAGS) -o main_PFA_Rader.exe main_PFA_Rader.o LEGACY.o -lntl -lgmp -lm	
main.o:main.cpp main_AE.cpp main_CF.cpp main_PFA_Rader.cpp

LEGACY.o: LEGACY.cpp LEGACY.h

clean:
	rm *.exe
	rm *.o
	rm *.txt

