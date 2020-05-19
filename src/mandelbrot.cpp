
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// mandelbrot.cpp                                                                                                                        //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - January, 2020                                                                                             //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// Producing fractals from the Mandelbrot Set


#include <fractal.hpp>

int iterate(std::complex<double> x, const std::complex<double> &c, const int max_itr)
{
  //
  // Iterate a given number under x^2 + c until its absolute value becomes greater than 2 or the maximum number of iterations is reached
  // If the absolute value of x does not become greater than 2, x is contained within the Mandelbrot Set
  //
  // parameters
  // ----------
  // x : std::complex<double>
  //  - the starting value of the iteration
  // c : const std::complex<double> &
  //  - the constant added at each iteration
  // max_itr : const int 
  //  - the maximum number of iterations allowed
  //
  // returns
  // -------
  // int
  //  - the number of iterations needed for the aboslute value of x became greater than 2 under iteration
  //    std::numeric_limits<int>::max() if the absolute value of x did not become greater in 2 within `max_itr` iterations
  //

  for (int itr=0;itr<max_itr;++itr)
  {
    if (std::abs(x)>2) { return itr; }
    x=x*x+c;
  }
  return std::numeric_limits<int>::max();
}

void compute_mandelbrot_range(int **iterations, const int max_itr, const int xresolution, const int start_itr, const int end_itr
  , const double startx, const double starty, const double deltax, const double deltay, const int total, bool verbose)
{
  //
  // Determine if a given set of numbers in a given subset of the complex plane are contained within the Mandelbrot Set through iteration
  // The given subset belongs to a larger space to actually be computed, handled by other threads
  //
  // parameters
  // ----------
  // iterations : int **
  //  - 2D array to write out either the number of iterations needed to show a number is not in the Mandelbrot Set, or a marker that this
  //    could not be shown
  // max_itr : const int
  //  - the maximum number of iterations allowed
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

  const std::complex<double> init(0.,0.);
  
  // loop over the given subset of imaginary components
  for (int itr=start_itr;itr<end_itr;++itr)
  {
    double imag=starty+deltay*itr;
    // loop over all real components
    for (int jtr=0;jtr<xresolution;++jtr)
    {
      double real=startx+deltax*jtr;
      // determine the number of iterations
      *(*(iterations+itr)+jtr)=iterate(init,std::complex<double>(real,imag),max_itr);
    }
    if (verbose&&itr%100==0&&itr!=0) { std::cout << "Processed " << itr*xresolution << " points of " << total << "." << std::endl; }
  }
}

int __declspec(dllexport) sample_mandelbrot(int **iterations, const int max_itr, const int num_threads, const int xresolution
  , const int yresolution, const double startx, const double endx, const double starty, const double endy
  , const bool verbose)
{
  //
  // Determine if numbers in a given subset of the complex plane are contained within the Mandelbrot Set through iteration
  //
  // parameters
  // ----------
  // iterations : int **
  //  - 2D array to write out either the number of iterations needed to show a number is not in the Mandelbrot Set, or a marker that this
  //    could not be shown
  // max_itr : const int
  //  - the maximum number of iterations allowed
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
  //  - the value which marks that a starting point could not be shown to not be in the Mandelbrot Set
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
      compute_mandelbrot_range,iterations,max_itr,xresolution,increments[itr],increments[itr+1],startx,starty,deltax,deltay,total
        ,num_threads>1 ? false : verbose
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
    std::cout << total << " points processed." << '\n'
      << "Time taken: " << elapsed.count() << "s." << std::endl;
  }

  return std::numeric_limits<int>::max();
}
