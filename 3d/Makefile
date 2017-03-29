.SUFFIXES: .cxx .o

CXX = g++ -g -Wall -Wfatal-errors -O3 -fopenmp

.cxx.o  :
	$(CXX) -c $? -o $@

all:
	@make kernel
	@make fmm

kernel: kernel.o
	$(CXX) $? -o $@
	./kernel 20
	./kernel 40
	./kernel 60

fmm: fmm.o
	$(CXX) $? -o $@
	./fmm

clean:
	$(RM) ./*.o ./kernel ./fmm
