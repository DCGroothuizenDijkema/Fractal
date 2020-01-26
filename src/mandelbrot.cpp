
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// mandelbrot.cpp                                                                                                                        //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - January, 2020                                                                                             //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// Producing fractals from the Mandelbrot Set


#include <fractal.hpp>
#include <iostream>
#include <chrono>

int iterate(const std::complex<double> &init, const std::complex<double> &c, const int max_itr)
{
  static int count=0;
  int itr;
  std::complex<double> x=init;

  for (itr=0;itr<max_itr;++itr)
  {
    if (abs(x)>2)
    {
      if (itr>count) { count=itr; }
      break;
    }
    x=x*x+c;
  }
  return itr+1;
}

void __declspec(dllexport) sample(int **iterations, const double x, const double y, const int max_itr, const int xresolution, 
  const int yresolution, const double startx, const double endx, const double starty, const double endy)
{
  const double deltax=(endx-startx)/xresolution,deltay=(endy-starty)/yresolution;
  const int total=xresolution*yresolution;

  std::complex<double> init(x,y);

  for (int itr=0;itr<xresolution;++itr)
  {
    double real=startx+deltax*itr;
    for (int jtr=0;jtr<yresolution;++jtr)
    {
      double imag=starty+deltay*jtr;
      *(*(iterations+itr)+jtr)=iterate(init,std::complex<double>(real,imag),max_itr);
    }
  }
}
