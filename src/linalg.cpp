
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// linalg.cpp                                                                                                                            //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - February, 2020                                                                                            //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// Linear algebra routines for fractals


#include <fractal.hpp>

void initialise_companion_matrix(std::complex<double> * const * const mat, const int degree)
{
  if (degree<2) { throw std::invalid_argument("`degree` must be greater than or equal to 2"); }
  
  for (int itr=0;itr<degree;++itr) { std::fill(*(mat+itr),*(mat+itr)+degree,0.); }
  for (int itr=1,jtr=0;itr<degree,jtr<degree-1;++itr,++jtr) { *(*(mat+itr)+jtr)=1.; }
}

void assign_companion_matrix(std::complex<double> * const * const mat, double const * const coeffs, const int degree)
{
  for (int itr=0;itr<degree;++itr) { *(*(mat+itr)+degree-1)=*(coeffs+itr); }
}

eigenpair<std::complex<double>> power_iteration(std::complex<double> const * const * const mat, std::complex<double> const * const vec
  , const size_t size, const double tol)
{
}

void deflate(std::complex<double> * const * const mat, const eigenpair<std::complex<double>> &pair)
{
}

void __declspec(dllexport) roots(double const * const coeffs, double * const roots_re, double * const roots_im, const int degree)
{
}
