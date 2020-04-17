
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

#define CUDA_REQUIRE_SUCCESS(ans) { cuda_check((ans),__FILE__,__LINE__); }
inline void cuda_check(const cudaError_t code, const char *file, const int line)
{
  if (code!=cudaSuccess)
  {
    std::cerr << "CUDA ERROR: " << cudaGetErrorString(code) << file << line << std::endl;
    exit(code);
  }
}

__device__ thrust::pair<cuDoubleComplex,cuDoubleComplex> polynomial_and_deriv(const cuDoubleComplex &x, const double * const coeffs
  , const int degree);
__device__ cuDoubleComplex newton_root(const double * const coeffs, int * const itr_taken, cuDoubleComplex x, const int degree
  , const int max_itr, const double tol);
__global__ void compute_newton(double * const d_re, double * const d_im, int * const d_itr, const double * const d_coeffs, const int max_itr
  , const int degree, const int xresolution, const int yresolution, const double startx, const double starty, const double deltax
  , const double deltay);

int __declspec(dllexport) sample_newton(double * const h_re, double * const h_im, int * const h_itr, const double * const h_coeffs
  , const int max_itr, const int degree, const int xresolution, const int yresolution, const double startx, const double endx
  , const double starty, const double endy, const bool verbose);

#endif // !CUFRACTAL_H__
