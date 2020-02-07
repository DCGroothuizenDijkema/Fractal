
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// linalg.cpp                                                                                                                            //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - February, 2020                                                                                            //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// Linear algebra routines for fractals


#include <fractal.hpp>

double **initialise_companion_matrix(const int degree)
{
  double **mat=new double*[degree];
  for (int itr=0;itr<degree;++itr)
  {
    *(mat+itr)=new double[degree];
    std::fill(*(mat+itr),*(mat+itr)+degree,0.);
  }
  for (int itr=1;itr<degree;++itr)
  {
    for (int jtr=0;jtr<degree-1;++jtr) { *(*(mat+itr)+jtr)=1.; }
  }
}

void assign_companion_matrix(double * const * const mat, double const * const coeffs, const int degree)
{
  for (int itr=0;itr<degree;++itr) { *(*(mat+itr)+degree)=*(coeffs+itr); }
}

void __declspec(dllexport) roots(double const * const coeffs, double * const roots_re, double * const roots_im, const int degree)
{
}
