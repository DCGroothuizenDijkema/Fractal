
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// cu_fractal.hpp                                                                                                                        //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - April, 2020                                                                                               //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// Header file for the CUDA accelerated implementation of the Fractal project


#pragma once

#ifndef CUFRACTAL_H__
#define CUFRACTAL_H__

#include <cuComplex.h>
#include <cuda_runtime.h>
#include <driver_types.h>
#include <math_constants.h>
#include <npp.h>

#include <thrust/pair.h>
#include <thrust/tuple.h>

#include <common.hpp>

#define CUDA_REQUIRE_SUCCESS(ans) { cuda_check((ans),__FILE__,__func__,#ans,__LINE__); }
inline void cuda_check(const cudaError_t code, const char * const file, const char * const func, const char * const call, const int line)
{
  if (code!=cudaSuccess)
  {
    std::cerr << file << ":" << line << ": CUDA ERROR (" << code << "): " << cudaGetErrorName(code) << ": " << cudaGetErrorString(code) 
      << '\n' << "  " << func << "()" << '\n'
      << "  {" << '\n'
      << "    " << call << '\n'
      << "  }" << std::endl;
    exit(code);
  }
}

__device__ int iterate(cuDoubleComplex x, const cuDoubleComplex &c, const int max_itr);
__global__ void compute_mandelbrot(int * const d_iterations, const int max_itr, const int xresolution, const int yresolution
  , const double startx, const double starty, const double deltax, const double deltay);
__global__ void compute_julia(int * const d_iterations, const double re, const double im, const int max_itr
  , const int xresolution, const int yresolution, const double startx, const double starty, const double deltax, const double deltay);

int sample_mandelbrot(int * const h_iterations, const int max_itr, const int xresolution, const int yresolution
  , const double startx, const double endx, const double starty, const double endy, const bool verbose);
int sample_julia(int * const h_itr, const double re, const double im, const int max_itr
  , const int xresolution, const int yresolution, const double startx, const double endx, const double starty, const double endy
  , const bool verbose);

__device__ thrust::pair<cuDoubleComplex,cuDoubleComplex> polynomial_and_deriv(const cuDoubleComplex &x, const double * const coeffs
  , const int degree);
__device__ cuDoubleComplex newton_root(const double * const coeffs, int * const itr_taken, cuDoubleComplex x, const int degree
  , const int max_itr, const double tol);
__global__ void compute_newton(double * const d_re, double * const d_im, int * const d_itr, const double * const d_coeffs, const int max_itr
  , const int degree, const int xresolution, const int yresolution, const double startx, const double starty, const double deltax
  , const double deltay);

int sample_newton(double * const h_re, double * const h_im, int * const h_itr, const double * const h_coeffs
  , const int max_itr, const int degree, const int xresolution, const int yresolution, const double startx, const double endx
  , const double starty, const double endy, const bool verbose);

#endif // !CUFRACTAL_H__
