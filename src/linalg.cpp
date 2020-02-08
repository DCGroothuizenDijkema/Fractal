
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// linalg.cpp                                                                                                                            //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - February, 2020                                                                                            //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// Linear algebra routines for fractals


#include <fractal.hpp>

void assign_companion_matrix(std::complex<double> * const * const mat, double const * const coeffs, const int degree)
{
  for (int itr=0;itr<degree;++itr) { *(*(mat+itr)+degree-1)=*(coeffs+itr); }
}

eigenpair<std::complex<double>> power_iteration(std::complex<double> const * const * const mat, std::complex<double> const * const vec
  , const size_t size, const double tol, const int max_itr)
{
  double delta;
  int itr=0;
  eigenpair<std::complex<double>> pair(vec,std::sqrt(dot(vec,vec,size)),size);

  pair()=pair.norm();
  pair.normalise();

  do
  {
    std::complex<double> prev_value=std::complex<double>(pair);
    // update guess
    std::complex<double> *product=new std::complex<double>[size];
    dot(mat,*pair,product,size);
    for (int jtr=0;jtr<size;++jtr) { pair[jtr]=*(product+jtr); }
    // normalise and update eigenvalue
    pair()=pair.norm();
    pair.normalise();
    // compute convergence
    delta=(std::abs(pair())-std::norm(prev_value))/(1+std::norm(pair()));
    if (itr==max_itr) { throw convergence_error(); }
  } while (fabs(delta)>=tol);

  return pair;
}

void deflate(std::complex<double> * const * const mat, const eigenpair<std::complex<double>> &pair)
{
}

void __declspec(dllexport) roots(double const * const coeffs, double * const roots_re, double * const roots_im, const int degree)
{
}
