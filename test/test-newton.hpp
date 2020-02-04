
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// test-newton.hpp                                                                                                                       //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - January, 2020                                                                                             //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// Test file for Newton's fractals


#include <random>
#include <iostream>

#include <fractal.hpp>

namespace NewtonTesting
{
BOOST_AUTO_TEST_SUITE(test_newton)

  BOOST_AUTO_TEST_CASE(horners_method)
  {
    double *coeffs=new double[1];
    *coeffs=0.;

    std::pair<std::complex<double>,std::complex<double>> result=polynomial_and_deriv(0.,coeffs,0);
    BOOST_CHECK(result.first==0.);
    BOOST_CHECK(result.second==0.);

    result=polynomial_and_deriv(std::complex<double>(4.,0.),coeffs,0);
    BOOST_CHECK(result.first==0.);
    BOOST_CHECK(result.second==0.);

    result=polynomial_and_deriv(std::complex<double>(0.,-1.),coeffs,0);
    BOOST_CHECK(result.first==0.);
    BOOST_CHECK(result.second==0.);

    result=polynomial_and_deriv(std::complex<double>(2.,3.),coeffs,0);
    BOOST_CHECK(result.first==0.);
    BOOST_CHECK(result.second==0.);
    
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> dist(-10,10);

    for (int itr=0;itr<1000;++itr)
    {
      *coeffs=dist(gen);
      result=polynomial_and_deriv(std::complex<double>(dist(gen),dist(gen)),coeffs,0);
      BOOST_CHECK(result.first==*coeffs);
      BOOST_CHECK(result.second==0.);
    }

    delete coeffs;
    coeffs=new double[2];
    *coeffs=1;
    *(coeffs+1)=1;

    result=polynomial_and_deriv(std::complex<double>(4.,0.),coeffs,1);
    BOOST_CHECK(result.first==5.);
    BOOST_CHECK(result.second==1.);

    result=polynomial_and_deriv(std::complex<double>(-2.,0.),coeffs,1);
    BOOST_CHECK(result.first==-1.);
    BOOST_CHECK(result.second==1.);

    for (int itr=0;itr<1000;++itr)
    {
      *coeffs=dist(gen);
      *(coeffs+1)=dist(gen);
      double i=dist(gen),j=dist(gen);
      result=polynomial_and_deriv(std::complex<double>(i,j),coeffs,1);
      BOOST_CHECK(result.first==*(coeffs+1)*std::complex<double>(i,j)+*coeffs);
      BOOST_CHECK(result.second==*(coeffs+1));
    }
  }

BOOST_AUTO_TEST_SUITE_END()
} // namespace NewtonTesting
