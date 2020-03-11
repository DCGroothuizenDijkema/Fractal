
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// fractal.hpp                                                                                                                           //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - January, 2020                                                                                             //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// Header file for Fractal project


#pragma once

#ifndef FRACTAL_H__
#define FRACTAL_H__

#include <cmath>

#include <algorithm>
#include <chrono>
#include <complex>
#include <iostream>
#include <limits>
#include <thread>
#include <type_traits>
#include <vector>

template <typename OutputIt>
inline OutputIt iteration_limits(const int num_threads, const int yresolution, OutputIt out)
{
  const int span=yresolution/num_threads;

  for (int itr=0;itr<num_threads;++itr) { *out++=span*itr; }
  *out++=yresolution;

  return out;
}

template <typename InputIt, typename OutputIt>
OutputIt zip(InputIt first1, InputIt last1, InputIt first2, InputIt last2, OutputIt out)
{
  using scalar=std::decay_t<decltype(*first1)>;
  while (first1!=last1&&first2!=last2)
  {
    *out++=std::complex<double>(static_cast<scalar>(*first1++),static_cast<scalar>(*first2++));
  }

  return out;
}

template <typename It, typename Compare>
std::size_t argmin(It first, It last, Compare comp)
{
  return std::min_element(first,last,comp)-first;
}


int iterate(std::complex<double> x, const std::complex<double> &c, const int max_itr);
void compute_mandelbrot_range(int **iterations, const int max_itr, const int xresolution, const int start_itr, const int end_itr
  , const double startx, const double starty, const double deltax, const double deltay, const int total, bool verbose);
void __declspec(dllexport) sample_mandelbrot(int **iterations, const int max_itr, const int num_threads, const int xresolution
  , const int yresolution, int * const limit, const double startx, const double endx, const double starty, const double endy
  , const bool verbose);

std::pair<std::complex<double>,std::complex<double>> polynomial_and_deriv(const std::complex<double> &x, double const * const coeffs
  , const int degree);
std::complex<double> newton_root( double const * const coeffs, int * const itr_taken, std::complex<double> x, const int degree
  , const int max_itr, const double tol);
void compute_newton_range(double **re, double **im, int **iterations, double * const coeffs, const int max_itr, const int degree
  , const int xresolution, const int start_itr, const int end_itr, const double startx, const double starty, const double deltax
  , const double deltay, const int total, bool verbose);
void __declspec(dllexport) sample_newton(double **re, double **im, int **iterations, double *coeffs, const int max_itr
  , const int num_threads, const int degree, const int xresolution, const int yresolution, int * const limit, const double startx
  , const double endx, const double starty, const double endy, const bool verbose);
void __declspec(dllexport) assign_roots(int * const * const index, const double * const * const re, const double * const * const im
  , const double * const roots_re, const double * const roots_im, const int degree, const int xresolution, const int yresolution);

#endif // FRACTAL_H__
