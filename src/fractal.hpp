
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

#include <chrono>
#include <complex>
#include <iostream>
#include <limits>
#include <thread>
#include <tuple>
#include <vector>

template <typename T>
struct eigenpair
{
  eigenpair(const size_t size);
  eigenpair(const T value, const size_t size);
  eigenpair(const eigenpair &other);
  eigenpair(eigenpair &&other);

  ~eigenpair();

  // assignment
  eigenpair &operator=(eigenpair other);
  // getters and setters
  [[nodiscard]] T &operator[](size_t idx);
  [[nodiscard]] const T &operator[](const size_t idx);

  // swap
  friend void swap(eigenpair& first, eigenpair& second) noexcept;
private:
  T value;
  T* vector;
  size_t size;
};

inline std::vector<int> iteration_limits(const int num_threads, const int yresolution)
{
  const int span=yresolution/num_threads;
  std::vector<int> increments;

  for (int itr=0;itr<num_threads;++itr) { increments.push_back(span*itr); }
  increments.push_back(yresolution);

  return increments;
}

template <typename T>
T dot(T const * const vector_one, T const * const vector_two, const size_t size);
template <typename T>
void dot(T const * const * const mat, T const * const vector, T * const out, const size_t size);

int iterate(std::complex<double> x, const std::complex<double> &c, const int max_itr);
void compute_mandelbrot_range(int **iterations, const int max_itr, const int xresolution, const int start_itr, const int end_itr
  , const double startx, const double starty, const double deltax, const double deltay, const int total, bool verbose);
void __declspec(dllexport) sample_mandelbrot(int **iterations, const int max_itr, const int num_threads, const int xresolution
  , const int yresolution, int * const limit, const double startx, const double endx, const double starty, const double endy
  , const bool verbose);

std::pair<std::complex<double>,std::complex<double>> polynomial_and_deriv(const std::complex<double> &x, double * const coeffs
  , const int degree);
std::complex<double> newton_root(double * const coeffs, int * const itr_taken, const std::complex<double> x, const int degree
  , const int max_itr, const double tol);
void compute_newton_range(double **re, double **im, int **iterations, double * coeffs, const int max_itr, const int degree
  , const int xresolution, const int start_itr, const int end_itr, const double startx, const double starty, const double deltax
  , const double deltay, const int total, bool verbose);
void __declspec(dllexport) sample_newton(double **re, double **im, int **iterations, double * coeffs, const int max_itr
  , const int num_threads, const int degree, const int xresolution, const int yresolution, int * const limit, const double startx
  , const double endx, const double starty, const double endy, const bool verbose);

void initialise_companion_matrix(double * const * const mat, const int degree);
void assign_companion_matrix(double * const * const mat, double const * const coeffs, const int degree);
void __declspec(dllexport) roots(double const * const coeffs, const int degree);

#endif // FRACTAL_H__
