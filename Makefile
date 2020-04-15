
# Windows nmake makefile

CXX=cl
LINK=link
CU=nvcc

FLAGS=/EHsc /std:c++17 /I./src/ /c /Fo:./obj/ /O2
TESTFLAGS=/EHsc /std:c++17 /I/lib/boost/ /I./src/ /I./test/ /c /Fo:./obj/

SRC=./src/mandelbrot.cpp ./src/newton.cpp
OBJ=./obj/mandelbrot.obj ./obj/newton.obj
CUOBJ=./obj/cumandelbrot.obj ./obj/cunewton.obj

INC=./src/fractal.hpp
TESTINC=./test/test-newton.hpp

TESTSRC=./test/test.cpp ./src/newton.cpp
TESTOBJ=./obj/test.obj ./obj/newton.obj

TARGET=./bin/fractal.dll
CUTARGET=./bin/cufractal.dll
TEST=./bin/test.exe

all: dir $(TARGET)
test: dir $(TEST)
cuda: dir $(CUTARGET)

dir: 
	-@ if NOT EXIST "./bin/" mkdir "./bin/"
	-@ if NOT EXIST "./obj/" mkdir "./obj"

$(TARGET):	$(OBJ)
	$(LINK) /DLL /OUT:$(TARGET) $(OBJ)

$(TEST):	$(TESTOBJ)
	$(LINK) /OUT:$(TEST) $(TESTOBJ)

$(CUTARGET):	$(CUOBJ)
	$(CU) -o $(CUTARGET) --shared $(CUOBJ)

obj/mandelbrot.obj: ./src/mandelbrot.cpp $(INC)
	$(CXX) $(FLAGS) ./src/mandelbrot.cpp

obj/newton.obj: ./src/newton.cpp $(INC)
	$(CXX) $(FLAGS) ./src/newton.cpp

obj/cumandelbrot.obj: ./src/mandelbrot.cpp $(INC)
	$(CU) -c -o ./obj/cumandelbrot.obj -I./src/ ./src/mandelbrot.cpp

obj/cunewton.obj: ./src/newton.cpp $(INC)
	$(CU) -c -o ./obj/cunewton.obj -I./src/ ./src/newton.cpp

obj/test.obj: ./test/test.cpp $(INC) $(TESTINC)
  $(CXX) $(TESTFLAGS) ./test/test.cpp
