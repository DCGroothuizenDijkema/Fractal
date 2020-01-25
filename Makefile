
# Windows nmake makefile

CXX=cl
LINK=link

FLAGS=/EHsc /std:c++17 /I/lib/boost/ /I./src/ /c /Fo:./obj/

SRC=./src/newton.cpp ./src/mandelbrot.cpp
OBJ=./obj/newton.obj ./obj/mandelbrot.obj

INC=./src/fractal.hpp

TARGET=fractal.dll

all: dir $(TARGET)
test: dir $(TEST)

dir: 
	-@ if NOT EXIST "./bin/" mkdir "./bin/"
	-@ if NOT EXIST "./obj/" mkdir "./obj"

$(TARGET):	$(OBJ)
	$(LINK) /DLL /OUT:./bin/$(TARGET) $(OBJ)

obj/newton.obj: ./src/newton.cpp $(INC)
	$(CXX) $(FLAGS) ./src/newton.cpp

obj/mandelbrot.obj: ./src/newton.cpp $(INC)
	$(CXX) $(FLAGS) ./src/mandelbrot.cpp
