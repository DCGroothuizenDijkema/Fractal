
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// cu_mandelbrot.cu                                                                                                                      //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - April, 2020                                                                                               //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// CUDA accelerated implementation for producing fractals from the Mandelbrot Set


#include <cu_fractal.hpp>

__device__ int iterate(cuDoubleComplex x, const cuDoubleComplex &c, const int max_itr)
{
  for (int itr=0;itr<max_itr;++itr)
  {
    if (cuCabs(x)>2) { return itr; }
    x=cuCadd(cuCmul(x,x),c);
  }
  return NPP_MAX_32S;
}

__global__ void compute_mandelbrot(int * const d_iterations, const int max_itr, const int xresolution, const int yresolution
  , const double startx, const double starty, const double deltax, const double deltay, const int total, bool verbose);

int __declspec(dllexport) sample_mandelbrot(int * const h_iterations, const int max_itr, const int xresolution, const int yresolution
  , const double startx, const double endx, const double starty, const double endy, const bool verbose);
