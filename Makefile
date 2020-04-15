
# Windows nmake makefile

CXX=cl
LINK=link
CU=nvcc

FLAGS=/EHsc /std:c++17 /I./src/ /c /Fo:./obj/ /O2
TESTFLAGS=/EHsc /std:c++17 /I/lib/boost/ /I./src/ /I./test/ /c /Fo:./obj/

SRC=./src/mandelbrot.cpp ./src/newton.cpp
CUSRC=./src/cu_newton.cu
OBJ=./obj/mandelbrot.obj ./obj/newton.obj
# CUOBJ=./obj/cu_mandelbrot.obj ./obj/cu_newton.obj
CUOBJ=./obj/cu_newton.obj

INC=./src/fractal.hpp
CUINC=./src/cu_fractal.hpp
TESTINC=./test/test-newton.hpp

TESTSRC=./test/test.cpp ./src/newton.cpp
TESTOBJ=./obj/test.obj ./obj/newton.obj

TARGET=./bin/fractal.dll
CUTARGET=./bin/cufractal.dll
TEST=./bin/test.exe

dll: dir $(TARGET)
test: dir $(TEST)
cuda: dir $(CUTARGET)

clean:
	-@ del /F /Q /S "./bin/" > NUL
	-@ rmdir /Q /S "./bin/"
	-@ del /F /Q /S "./obj/" > NUL
	-@ rmdir /Q /S "./obj/"

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

# obj/cu_mandelbrot.obj: ./src/mandelbrot.cpp $(INC)
	# $(CU) -c -o ./obj/cumandelbrot.obj -I./src/ ./src/mandelbrot.cpp

obj/cu_newton.obj: ./src/cu_newton.cu $(CUINC)
	$(CU) -c -o ./obj/cu_newton.obj -I./src/ ./src/cu_newton.cu

obj/test.obj: ./test/test.cpp $(INC) $(TESTINC)
  $(CXX) $(TESTFLAGS) ./test/test.cpp
