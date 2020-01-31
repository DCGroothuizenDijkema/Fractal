
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// mandelbrot.cpp                                                                                                                        //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - January, 2020                                                                                             //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// Producing fractals from the Mandelbrot Set


#include <fractal.hpp>

int iterate(const std::complex<double> &init, const std::complex<double> &c, const int max_itr)
{
  static int count=0;
  int itr;
  std::complex<double> x=init;

  for (itr=0;itr<max_itr;++itr)
  {
    if (abs(x)>2)
    {
      return ++itr;
    }
    x=x*x+c;
  }
  return std::numeric_limits<int>::max();
}

void __declspec(dllexport) sample_mandelbrot(int **iterations, const int max_itr, const int xresolution, const int yresolution, int * const limit,
  const double startx, const double endx, const double starty, const double endy, const bool verbose)
{
  const double deltax=(endx-startx)/xresolution,deltay=(endy-starty)/yresolution;
  const int total=xresolution*yresolution;

  if (verbose) { std::cout << "Processing " << total << " points." << std::endl; }
  const std::complex<double> init(0.,0.);
  std::chrono::time_point<std::chrono::steady_clock> start=std::chrono::high_resolution_clock::now();

  for (int itr=0;itr<yresolution;++itr)
  {
    double imag=starty+deltay*itr;
    for (int jtr=0;jtr<xresolution;++jtr)
    {
      double real=startx+deltax*jtr;
      *(*(iterations+itr)+jtr)=iterate(init,std::complex<double>(real,imag),max_itr);
    }
    if (itr%100==0&&itr!=0) if (verbose) { std::cout << "Processed " << itr*xresolution << " points of " << total << "." << std::endl; }
  }
  
  std::chrono::time_point<std::chrono::steady_clock> finish=std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed=finish-start;
  if (verbose) { std::cout << total << " points processed." << std::endl; }
  if (verbose) { std::cout << "Time taken: " << elapsed.count() << "s." << std::endl; }

  *limit=std::numeric_limits<int>::max();
}
