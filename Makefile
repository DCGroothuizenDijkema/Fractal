
# Windows nmake makefile

CXX=cl
LINK=link
CU=nvcc

FLAGS=/EHsc /std:c++17 /I./src/ /c /Fo:./obj/ /O2
TESTFLAGS=/EHsc /std:c++17 /I/lib/boost/ /I./src/ /I./test/ /c /Fo:./obj/

SRC=./src/mandelbrot.cpp ./src/newton.cpp
CUSRC=./src/cu_newton.cu ./src/cu_mandelbrot.cu
COMMONSRC=./src/newton_common.cpp

OBJ=./obj/mandelbrot.obj ./obj/newton.obj
CUOBJ=./obj/cu_newton.obj ./obj/cu_mandelbrot.obj
COMMONOBJ=./obj/newton_common.obj

INC=./src/fractal.hpp
CUINC=./src/cu_fractal.hpp
COMMONINC=./src/common.hpp
TESTINC=./test/test-newton.hpp

TESTSRC=./test/test.cpp ./src/newton.cpp
TESTOBJ=./obj/test.obj ./obj/newton.obj

TARGET=./bin/fractal.dll
CUTARGET=./bin/cufractal.dll
TEST=./bin/test.exe

dll: dir $(TARGET)
test: dir $(TEST)
cuda: dir $(CUTARGET)
all: dir $(TARGET) $(CUTARGET) $(TEST)
reset: clean all

clean:
	-@ if EXIST "./bin/" del /F /Q /S "./bin/" > NUL
	-@ if EXIST "./bin/" rmdir /Q /S "./bin/"
	-@ if EXIST "./obj/" del /F /Q /S "./obj/" > NUL
	-@ if EXIST "./obj/" rmdir /Q /S "./obj/"

dir: 
	-@ if NOT EXIST "./bin/" mkdir "./bin/"
	-@ if NOT EXIST "./obj/" mkdir "./obj"

$(TARGET):	$(OBJ) $(COMMONOBJ)
	$(LINK) /DLL /OUT:$(TARGET) $(OBJ) $(COMMONOBJ)

$(TEST):	$(TESTOBJ)
	$(LINK) /OUT:$(TEST) $(TESTOBJ)

$(CUTARGET):	$(CUOBJ) $(COMMONOBJ)
	$(CU) -o $(CUTARGET) --shared $(CUOBJ) $(COMMONOBJ)

obj/mandelbrot.obj: ./src/mandelbrot.cpp $(INC) $(COMMONINC)
	$(CXX) $(FLAGS) ./src/mandelbrot.cpp

obj/newton.obj: ./src/newton.cpp $(INC) $(COMMONINC)
	$(CXX) $(FLAGS) ./src/newton.cpp

obj/cu_mandelbrot.obj: ./src/cu_mandelbrot.cu $(CUINC) $(COMMONINC)
  $(CU) -c -o ./obj/cu_mandelbrot.obj -I./src/ ./src/cu_mandelbrot.cu

obj/cu_newton.obj: ./src/cu_newton.cu $(CUINC) $(COMMONINC)
	$(CU) -c -o ./obj/cu_newton.obj -I./src/ ./src/cu_newton.cu

obj/newton_common.obj: ./src/cu_newton.cu $(COMMONINC)
	$(CXX) $(FLAGS) ./src/newton_common.cpp

obj/test.obj: ./test/test.cpp $(INC) $(TESTINC)
  $(CXX) $(TESTFLAGS) ./test/test.cpp
