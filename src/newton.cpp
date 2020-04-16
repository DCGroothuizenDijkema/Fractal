
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// newton.cpp                                                                                                                            //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - January, 2020                                                                                             //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// Producing Newton's fractals


#include <fractal.hpp>

std::pair<std::complex<double>,std::complex<double>> polynomial_and_deriv(const std::complex<double> &x, const double * const coeffs
  , const int degree)
{
  //
  // Evaluate the value and derivative of a polynomial at a point using Horner's method
  //
  // parameters
  // ----------
  // x : const std::complex<double> &
  //  - the point to evaluate at
  // coeffs : const double * const
  //  - the coefficients of the polynomial given in order of the lowest degree to highest
  //    if the polynomial is written p(x)=a_n*x^(n-1)+...+a_(k)*^(n-k-1)+...+a_1*x+a_0, then coeffs should have the following form:
  //      *(coeffs+0)==a_0
  //      *(coeffs+1)==a_1
  //      *(coeffs+k)==a_k
  //      *(coeffs+n)==a_n
  // degree : const int
  //  - the degree of the polynomial
  //
  // returns
  // -------
  // std::pair<std::complex<double>,std::complex<double>>
  //  - the function value and derivative of the given polynomial evaluated at the given point
  //

  std::complex<double> p=*(coeffs+degree),p_prime=0;

  for (int itr=degree-1;itr>=0;--itr)
  {
    p_prime=x*p_prime+p;
    p=x*p+*(coeffs+itr);
  }

  return std::make_pair(p,p_prime);
}

std::complex<double> newton_root(const double * const coeffs, int * const itr_taken, std::complex<double> x, const int degree
  , const int max_itr, const double tol)
{
  //
  // Apply the Newton-Raphson method to a given number to find the root of a given polynomial
  //
  // parameters
  // ----------
  // coeffs : const double * const
  //  - the coefficients of the polynomial given in order of the lowest degree to highest
  //    see polynomial_and_deriv() for requirements
  // itr_taken : int * const 
  //  - a poitner to write out the number of iterations needed to reach a root
  // x : std::complex<double>
  //  - the number to a start from
  // degree : const int
  //  - the degree of the polynomial
  // max_itr : const int
  //  - the maximum number of iterations allowed
  // tol : const double
  //  - the tolerance within which a root will be deemed to have been reached
  //
  // returns
  // -------
  // std::complex<double>
  //  - the root reached from the given number
  //

  for (int itr=0;itr<max_itr;++itr)
  {
    // get the current function value and derivative
    std::complex<double> f_x,g_x;
    std::tie(f_x,g_x)=polynomial_and_deriv(x,coeffs,degree);
    // converged to a root
    if (abs(f_x)<tol)
    {
      *itr_taken=itr;
      return x;
    }
    // derivative is flat and we can't update
    if (g_x==std::complex<double>(0.,0.))
    {
      *itr_taken=std::numeric_limits<int>::max();
      return std::complex<double>(std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity());
    }
    // update
    x-=f_x/g_x;
  }
  // couldn't find a root in the given number of iterations
  *itr_taken=std::numeric_limits<int>::max();
  return std::complex<double>(std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity());
}

void compute_newton_range(double **re, double **im, int **iterations, double * const coeffs, const int max_itr, const int degree
  , const int xresolution, const int start_itr, const int end_itr, const double startx, const double starty, const double deltax
  , const double deltay, const int total, bool verbose)
{
  //
  // Find the roots of a polynomial a given set of numbers in a given subset of the complex plane converge to with Newton's method
  // The given subset belongs to a larger space to actually be computed, handled by other threads
  //
  // parameters
  // ----------
  // re,im : double **
  //  - 2D arrays to write out the root which the numbers in the subset cnoverged to
  // iterations : int **
  //  - 2D array to write out either the number of iterations needed reach a root or a marker that this root could not be reached
  // coeffs : double * const
  //  - the coefficients of the polynomial given in order of the lowest degree to highest
  //    see polynomial_and_deriv() for requirements
  // max_itr : const int
  //  - the maximum number of iterations allowed
  // degree : const int
  //  - the degree of the polynomial
  // xresolution : const int 
  //  - the number of steps to take in the x-direction (the real components)
  // start_itr,end_itr : const int 
  //  - the subset of the steps taken from the starting point in the y-direction (the complex component)
  // startx,starty : const int 
  //  - the real and imaginary components of the number defining the bottom left corner of the entire space being sampled
  // deltax,deltay : const int 
  //  - the size of the step to take in the x- and y-direction
  // total : const int 
  //  - the total number of points being processed
  //    only used if `verbose` is true.
  // verbose : bool
  //  - flag to control logging to console
  //

  // loop over the given subset of imaginary components
  for (int itr=start_itr;itr<end_itr;++itr)
  {
    double imag=starty+deltay*itr;
    // loop over all real components
    for (int jtr=0;jtr<xresolution;++jtr)
    {
      double real=startx+deltax*jtr;
      // determine the root reached and the number of iterations to get there
      std::complex<double> root=newton_root(coeffs,(*(iterations+itr)+jtr),std::complex<double>(real,imag),degree,max_itr,1e-6);
      *(*(re+itr)+jtr)=root.real();
      *(*(im+itr)+jtr)=root.imag();
    }
    if (verbose&&itr%100==0&&itr!=0) { std::cout << "Processed " << itr*xresolution << " points of " << total << "." << std::endl; }
  }
}

int __declspec(dllexport) sample_newton(double **re, double **im, int **iterations, double *coeffs, const int max_itr
  , const int num_threads, const int degree, const int xresolution, const int yresolution, const double startx
  , const double endx, const double starty, const double endy, const bool verbose)
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
  // returns
  // -------
  // int
  //  - the value which marks that a starting point could not be shown to not be in the mandelbrot set
  //

  // determine step sizing in the x- and y-direction
  const double deltax=(endx-startx)/xresolution,deltay=(endy-starty)/yresolution;
  const int total=xresolution*yresolution;

  // get the split in the y-direction to allocate to the threads
  std::vector<int> increments;
  iteration_limits(num_threads,yresolution,std::back_inserter(increments));

  // construct all threads, where each thread takes a new left-closed, right-open interval of the imaginary component of the plane
  std::vector<std::thread> threads;
  for (int itr=0;itr<increments.size()-1;++itr)
  {
    threads.push_back(std::thread(
      compute_newton_range,re,im,iterations,std::ref(coeffs),max_itr,degree,xresolution,increments[itr],increments[itr+1],startx,starty
        ,deltax,deltay,total,num_threads>1 ? false : verbose
    ));
  }

  if (verbose) { std::cout << "Processing " << total << " points." << std::endl; }

  // execute all threads
  std::chrono::time_point<std::chrono::steady_clock> start=std::chrono::high_resolution_clock::now();
  for (std::thread &th:threads) { th.join(); }
  std::chrono::time_point<std::chrono::steady_clock> finish=std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed=finish-start;
  if (verbose)
  { 
    std::cout << total << " points processed." << std::endl;
    std::cout << "Time taken: " << elapsed.count() << "s." << std::endl;
  }

  return std::numeric_limits<int>::max();
}

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
