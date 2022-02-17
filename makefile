
all: main 

build: LEGACY.o

main: main_trash.o main.o main_reconfig.o main_winograd.o main_reconfig_0214.o build
	g++ $(LDFLAGS) -o main_trash.exe main_trash.o LEGACY.o -lntl -lgmp -lm			
	g++ $(LDFLAGS) -o main.exe main.o LEGACY.o -lntl -lgmp -lm	
	g++ $(LDFLAGS) -o main_reconfig.exe main_reconfig.o LEGACY.o -lntl -lgmp -lm	
	g++ $(LDFLAGS) -o main_winograd.exe main_winograd.o LEGACY.o -lntl -lgmp -lm	
	g++ $(LDFLAGS) -o main_reconfig_0214.exe main_reconfig_0214.o LEGACY.o -lntl -lgmp -lm		
main.o:main_trash.cpp main.cpp 

LEGACY.o: LEGACY.cpp LEGACY.h

clean:
	rm *.exe
	rm *.o
	rm *.txt

