
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
  int itr;

  for (itr=0;itr<max_itr;++itr)
  {
    if (abs(x)>2) { return itr; }
    x=x*x+c;
  }
  return std::numeric_limits<int>::max();
}

void compute_mandelbrot_range(int **iterations, const int max_itr, const int xresolution, const int start_itr, const int end_itr
  , const double startx, const double starty, const double deltax, const double deltay, const int total, bool verbose)
{
  const std::complex<double> init(0.,0.);
  
  for (int itr=start_itr;itr<end_itr;++itr)
  {
    double imag=starty+deltay*itr;
    for (int jtr=0;jtr<xresolution;++jtr)
    {
      double real=startx+deltax*jtr;
      *(*(iterations+itr)+jtr)=iterate(init,std::complex<double>(real,imag),max_itr);
    }
    if (verbose&&itr%100==0&&itr!=0) { std::cout << "Processed " << itr*xresolution << " points of " << total << "." << std::endl; }
  }
}

void __declspec(dllexport) sample_mandelbrot(int **iterations, const int max_itr, const int num_threads, const int xresolution
  , const int yresolution, int * const limit, const double startx, const double endx, const double starty, const double endy
  , const bool verbose)
{
  const double deltax=(endx-startx)/xresolution,deltay=(endy-starty)/yresolution;
  const int total=xresolution*yresolution;

  const std::vector<int> increments=iteration_limits(num_threads,yresolution);

  std::vector<std::thread> threads;
  for (int itr=0;itr<increments.size()-1;++itr)
  {
    threads.push_back(std::thread(
      compute_mandelbrot_range,iterations,max_itr,xresolution,increments[itr],increments[itr+1],startx,starty,deltax,deltay,total
        ,num_threads>1 ? false : verbose
    ));
  }

  if (verbose) { std::cout << "Processing " << total << " points." << std::endl; }

  std::chrono::time_point<std::chrono::steady_clock> start=std::chrono::high_resolution_clock::now();
  for (std::thread &th:threads) { th.join(); }
  std::chrono::time_point<std::chrono::steady_clock> finish=std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed=finish-start;
  if (verbose) { std::cout << total << " points processed." << std::endl; }
  if (verbose) { std::cout << "Time taken: " << elapsed.count() << "s." << std::endl; }

  *limit=std::numeric_limits<int>::max();
}
