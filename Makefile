
# Windows nmake makefile

CXX=cl
LINK=link

FLAGS=/EHsc /std:c++17 /I./src/ /c /Fo:./obj/ /O2
TESTFLAGS=/EHsc /std:c++17 /I/lib/boost/ /I./src/ /I./test/ /c /Fo:./obj/

SRC=./src/linalg.cpp ./src/mandelbrot.cpp ./src/newton.cpp
OBJ=./obj/linalg.obj ./obj/mandelbrot.obj ./obj/newton.obj

INC=./src/fractal.hpp
TESTINC=./test/test-linalg.hpp ./test/test-newton.hpp

TESTSRC=./test/test.cpp /src/linalg.cpp ./src/newton.cpp
TESTOBJ=./obj/test.obj ./obj/linalg.obj ./obj/newton.obj

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

obj/linalg.obj: ./src/linalg.cpp $(INC)
	$(CXX) $(FLAGS) ./src/linalg.cpp

obj/mandelbrot.obj: ./src/mandelbrot.cpp $(INC)
	$(CXX) $(FLAGS) ./src/mandelbrot.cpp

obj/newton.obj: ./src/newton.cpp $(INC)
	$(CXX) $(FLAGS) ./src/newton.cpp

obj/test.obj: ./test/test.cpp $(INC) $(TESTINC)
  $(CXX) $(TESTFLAGS) ./test/test.cpp
