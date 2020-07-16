
CXX=g++
LINK=ld
CU=nvcc

FLAGS=-std=c++17 -Wall -I./src/ -c -Ofast

SRC=./src/mandelbrot.cpp ./src/newton.cpp
CUSRC=./src/cu_newton.cu ./src/cu_mandelbrot.cu
COMMONSRC=./src/newton_common.cpp

OBJ=./obj/mandelbrot.o ./obj/newton.o
CUOBJ=./obj/cu_newton.o ./obj/cu_mandelbrot.o
COMMONOBJ=./obj/newton_common.o

INC=./src/fractal.hpp
CUINC=./src/cu_fractal.hpp
COMMONINC=./src/common.hpp

TARGET=./bin/fractal

all: $(TARGET)

$(TARGET):	$(OBJ) $(CUOBJ) $(COMMONOBJ)
	$(LINK) $(OBJ) $(CUOBJ) $(COMMONOBJ) 

obj/mandelbrot.o: ./src/mandelbrot.cpp $(INC) $(COMMONINC)
	$(CXX) $(FLAGS) ./src/mandelbrot.cpp

obj/newton.o: ./src/newton.cpp $(INC) $(COMMONINC)
	$(CXX) $(FLAGS) ./src/newton.cpp

obj/cu_mandelbrot.o: ./src/cu_mandelbrot.cu $(CUINC) $(COMMONINC)
	$(CU) -c -o ./obj/cu_mandelbrot.o -I./src/ ./src/cu_mandelbrot.cu

obj/cu_newton.o: ./src/cu_newton.cu $(CUINC) $(COMMONINC)
	$(CU) -c -o ./obj/cu_newton.o -I./src/ ./src/cu_newton.cu

obj/newton_common.o: ./src/cu_newton.cu $(COMMONINC)
	$(CXX) $(FLAGS) ./src/newton_common.cpp
