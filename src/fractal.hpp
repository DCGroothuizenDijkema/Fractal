
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
#include <tuple>

int iterate(std::complex<double> x, const std::complex<double> &c, const int max_itr);
void __declspec(dllexport) sample_mandelbrot(int **iterations, const int max_itr, const int xresolution, const int yresolution,
  int * const limit, const double startx, const double endx, const double starty, const double endy, const bool verbose);

std::pair<std::complex<double>,std::complex<double>> polynomial_and_deriv(const std::complex<double> &x, double * const coeffs
  , const int degree);

#endif // FRACTAL_H__
