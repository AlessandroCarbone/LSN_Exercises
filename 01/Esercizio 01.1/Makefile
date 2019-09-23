CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o
	g++ random.o main.o -o main.exe  
main.o : main.cpp
	g++ -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
