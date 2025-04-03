a.out: BLTC.o directsum.o quicksort.o main.o
	@echo "Linking outputs ..."
	nvcc -arch=native BLTC.o directsum.o quicksort.o main.o

BLTC.o: BLTC.cu
	@echo "Building BLTC ..."
	nvcc -arch=native --default-stream per-thread -dc BLTC.cu

directsum.o: directsum.cu
	@echo "Building directsum ..."
	nvcc -arch=native -rdc=true -dc directsum.cu -lcudadevrt

quicksort.o: quicksort.cpp
	@echo "Building quicksort ..."
	g++ -c quicksort.cpp

main.o: main.cpp
	@echo "Building main ..."
	g++ -c main.cpp

clean:
	rm $(wildcard *.o)
	rm a.out
