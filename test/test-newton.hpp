
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// test-newton.hpp                                                                                                                       //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - January, 2020                                                                                             //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// Test file for Newton's fractals


#include <fractal.hpp>

namespace NewtonTesting
{
BOOST_AUTO_TEST_SUITE(test_newton)

  BOOST_AUTO_TEST_CASE(horners_method_testing)
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
    *coeffs=1.;
    *(coeffs+1)=1.;

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

    delete coeffs;
    coeffs=new double[3];
    *coeffs=1.;
    *(coeffs+1)=-1.;
    *(coeffs+2)=2.;

    result=polynomial_and_deriv(std::complex<double>(1.,0.),coeffs,2);

    BOOST_CHECK(result.first==std::complex<double>(2.,0.));
    BOOST_CHECK(result.second==std::complex<double>(3.,0.));

    for (int itr=0;itr<1000;++itr)
    {
      *coeffs=dist(gen);
      *(coeffs+1)=dist(gen);
      *(coeffs+2)=dist(gen);

      double i=dist(gen),j=dist(gen);
      std::complex<double> expected_poly=*(coeffs+2)*std::complex<double>(i,j)*std::complex<double>(i,j)+
        *(coeffs+1)*std::complex<double>(i,j)+*coeffs,expected_deriv=2**(coeffs+2)*std::complex<double>(i,j)+*(coeffs+1);
      result=polynomial_and_deriv(std::complex<double>(i,j),coeffs,2);

      BOOST_CHECK_CLOSE_FRACTION(result.first.real(),expected_poly.real(),0.0001);
      BOOST_CHECK_CLOSE_FRACTION(result.first.imag(),expected_poly.imag(),0.0001);
      BOOST_CHECK_CLOSE_FRACTION(result.second.real(),expected_deriv.real(),0.0001);
      BOOST_CHECK_CLOSE_FRACTION(result.second.imag(),expected_deriv.imag(),0.0001);
    }
  }

  BOOST_AUTO_TEST_CASE(zip_testing)
  {
    const int degree=3;
    double one[degree]={1,2,3},two[degree]={10,11,12};
    std::vector<std::complex<double>> zipped;

    zip(one,one+degree,two,two+degree,std::back_inserter(zipped));

    BOOST_CHECK(zipped[0]==std::complex<double>(1,10));
    BOOST_CHECK(zipped[1]==std::complex<double>(2,11));
    BOOST_CHECK(zipped[2]==std::complex<double>(3,12));
  }

  BOOST_AUTO_TEST_CASE(argmin_testing)
  {
    std::vector<int> test_vec={4,9,1,3,2};

    BOOST_CHECK(argmin(test_vec,[](const std::complex<double> &x, const std::complex<double> &y){ return abs(x)<abs(y); })==2);
  }

BOOST_AUTO_TEST_SUITE_END()
} // namespace NewtonTesting
