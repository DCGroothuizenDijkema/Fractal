
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// newton_common.cpp                                                                                                                     //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - April, 2020                                                                                               //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// Functions common to both implementations of producing Newton's fractals


#include <common.hpp>

void __declspec(dllexport) assign_roots(int * const * const index, const double * const * const re, const double * const * const im
  , const double * const roots_re, const double * const roots_im, const int degree, const int xresolution, const int yresolution)
{
  //
  // Determine the roots of a polynomial the numbers in a given subset of the complex plane converge to with Newton's method
  //
  // parameters
  // ----------
  // re,im : double **
  //  - 2D arrays to write out the root which the numbers in the subset cnoverged to
  // iterations : int **
  //  - 2D array to write out either the number of iterations needed reach a root or a marker that this root could not be reached
  // coeffs : double *
  //  - the coefficients of the polynomial given in order of the lowest degree to highest
  //    see polynomial_and_deriv() for requirements
  // max_itr : const int
  //  - the maximum number of iterations allowed
  // degree : const int
  //  - the degree of the polynomial
  // num_threads : const int
  //  - the number of threads to use in computation
  // xresolution,yresolution : const int 
  //  - the number of steps to take in the x- and y-direction (the real and imaginary components)
  // startx,endx,starty,endy : const int
  //  - the first and last values to sample at 
  // verbose : bool
  //  - flag to control logging to console
  //

  // get a list of all roots, formed from the input vectors giving the real and imaginary components of the roots
  std::vector<std::complex<double>> roots;
  zip(roots_re,roots_re+degree,roots_im,roots_im+degree,std::back_inserter(roots));
  
  for (int itr=0;itr<yresolution;++itr)
  {
    for (int jtr=0;jtr<xresolution;++jtr)
    {
      // if the current value is marked with infinity, no root was reached from it and its index is 0
      if (*(*(re+itr)+jtr)==std::numeric_limits<double>::infinity()||*(*(im+itr)+jtr)==std::numeric_limits<double>::infinity())
      {
        *(*(index+itr)+jtr)=-1;
        continue;
      }

      std::complex<double> val(*(*(re+itr)+jtr),*(*(im+itr)+jtr));
      // determine the difference between the current value and all roots
      std::vector<std::complex<double>> diffs;
      std::transform(std::begin(roots),std::end(roots),std::back_inserter(diffs)
        ,[val](std::complex<double> root) { return abs(root-val); });
      // find the argmin of the differences, and, therefore, which root was converged to
      *(*(index+itr)+jtr)=static_cast<int>(argmin(std::cbegin(diffs),std::cend(diffs)
        ,[](const std::complex<double> &x, const std::complex<double> &y){ return abs(x)<abs(y); }
      ));
    }
  }
}
