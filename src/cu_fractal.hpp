
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
#include <math_constants.h>
#include <cuda_runtime.h>
#include <driver_types.h>
#include <npp.h>
#include <thrust/pair.h>
#include <thrust/tuple.h>

#include <common.hpp>

#define CUDA_ASSERT_SUCCESS(ans) { _cuda_assert_success((ans),__FILE__,__LINE__); }
inline void _cuda_assert_success(cudaError_t code, const char *file, int line)
{
  if (code!=cudaSuccess)
  {
    std::cerr << "CUDA ERROR: " << cudaGetErrorString(code) << file << line << std::endl;;
    exit(code);
  }
}

__device__  thrust::pair<cuDoubleComplex,cuDoubleComplex> polynomial_and_deriv(const cuDoubleComplex &x, const double * const coeffs
  , const int degree);
__device__  cuDoubleComplex newton_root(const double * const coeffs, int * const itr_taken, cuDoubleComplex x, const int degree
  , const int max_itr, const double tol);
__global__ void compute_newton(double *d_re, double *d_im, int *d_itr, double * const coeffs, const int max_itr, const int degree
  , const int xresolution, const int yresolution, const double startx, const double starty, const double deltax, const double deltay);

int __declspec(dllexport) sample_newton(double *h_re, double *h_im, int *h_itr, double *coeffs, const int max_itr, const int degree
  , const int xresolution, const int yresolution, const double startx, const double endx, const double starty, const double endy);
void __declspec(dllexport) assign_roots(int * const * const index, const double * const * const re, const double * const * const im
  , const double * const roots_re, const double * const roots_im, const int degree, const int xresolution, const int yresolution);

#endif // !CUFRACTAL_H__
