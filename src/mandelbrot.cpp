
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
