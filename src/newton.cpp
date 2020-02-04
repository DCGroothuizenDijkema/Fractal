
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// newton.cpp                                                                                                                            //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - January, 2020                                                                                             //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// Producing Newton's fractals


#include <fractal.hpp>

std::pair<std::complex<double>,std::complex<double>> polynomial_and_deriv(const std::complex<double> &x, double * const coeffs
  , const int degree)
{
  std::complex<double> p=*(coeffs+degree),p_prime=0;

  for (int itr=degree-1;itr>=0;--itr)
  {
    p_prime=x*p_prime+p;
    p=x*p+*(coeffs+itr);
  }

  return std::make_pair(p,p_prime);
}
