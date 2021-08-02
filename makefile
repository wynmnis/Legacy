
all: main 

build: LEGACY.o

main: main.o main_AE.o main_CF.o main_PFA_Rader.o main_trash.o main_21845p.o main_21845_ZZ.o main_test_config_rader.o build
	g++ $(LDFLAGS) -o main.exe  main.o  LEGACY.o -lntl -lgmp -lm
	g++ $(LDFLAGS) -o main_AE.exe main_AE.o LEGACY.o -lntl -lgmp -lm
	g++ $(LDFLAGS) -o main_CF.exe main_CF.o LEGACY.o -lntl -lgmp -lm	
	g++ $(LDFLAGS) -o main_PFA_Rader.exe main_PFA_Rader.o LEGACY.o -lntl -lgmp -lm	
	g++ $(LDFLAGS) -o main_trash.exe main_trash.o LEGACY.o -lntl -lgmp -lm		
	g++ $(LDFLAGS) -o main_21845p.exe main_21845p.o LEGACY.o -lntl -lgmp -lm	
	g++ $(LDFLAGS) -o main_21845_ZZ.exe main_21845_ZZ.o LEGACY.o -lntl -lgmp -lm	
	g++ $(LDFLAGS) -o main_test_config_rader.exe main_test_config_rader.o LEGACY.o -lntl -lgmp -lm	
	
main.o:main.cpp main_AE.cpp main_CF.cpp main_PFA_Rader.cpp main_trash.cpp main_21845p.cpp main_21845_ZZ.cpp main_test_config_rader.cpp

LEGACY.o: LEGACY.cpp LEGACY.h

clean:
	rm *.exe
	rm *.o
	rm *.txt

