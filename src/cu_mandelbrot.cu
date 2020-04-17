
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// cu_mandelbrot.cu                                                                                                                      //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - April, 2020                                                                                               //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// CUDA accelerated implementation for producing fractals from the Mandelbrot Set


#include <cu_fractal.hpp>

__device__ int iterate(cuDoubleComplex x, const cuDoubleComplex &c, const int max_itr);

__global__ void compute_mandelbrot(int * const d_iterations, const int max_itr, const int xresolution, const int yresolution
  , const double startx, const double starty, const double deltax, const double deltay, const int total, bool verbose);

int __declspec(dllexport) sample_mandelbrot(int * const h_iterations, const int max_itr, const int xresolution, const int yresolution
  , const double startx, const double endx, const double starty, const double endy, const bool verbose);
