
# Windows nmake makefile

CXX=cl
LINK=link

FLAGS=/EHsc /std:c++17 /I./src/ /c /Fo:./obj/ /O2
TESTFLAGS=/EHsc /std:c++17 /I/lib/boost/ /I./src/ /I./test/ /c /Fo:./obj/

SRC=./src/newton.cpp ./src/mandelbrot.cpp
OBJ=./obj/newton.obj ./obj/mandelbrot.obj

INC=./src/fractal.hpp
TESTINC=./test/test-newton.hpp

TESTSRC=./test/test.cpp ./src/newton.cpp
TESTOBJ=./obj/test.obj ./obj/newton.obj

TARGET=fractal.dll
TEST=test.exe

all: dir $(TARGET)
test: dir $(TEST)

dir: 
	-@ if NOT EXIST "./bin/" mkdir "./bin/"
	-@ if NOT EXIST "./obj/" mkdir "./obj"

$(TARGET):	$(OBJ)
	$(LINK) /DLL /OUT:./bin/$(TARGET) $(OBJ)

$(TEST):	$(TESTOBJ)
	$(LINK) /OUT:./bin/$(TEST) $(TESTOBJ)

obj/newton.obj: ./src/newton.cpp $(INC)
	$(CXX) $(FLAGS) ./src/newton.cpp

obj/mandelbrot.obj: ./src/mandelbrot.cpp $(INC)
	$(CXX) $(FLAGS) ./src/mandelbrot.cpp

obj/test.obj: ./test/test.cpp $(INC) $(TESTINC)
  $(CXX) $(TESTFLAGS) ./test/test.cpp
