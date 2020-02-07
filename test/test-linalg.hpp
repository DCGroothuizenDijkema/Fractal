
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// test-linalg.hpp                                                                                                                       //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - February, 2020                                                                                            //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// Test file for Newton's fractals


#include <random>

#include <fractal.hpp>

namespace LinalgTesting
{
BOOST_AUTO_TEST_SUITE(test_linalg)

  BOOST_AUTO_TEST_CASE(companion_matrix)
  {
    BOOST_CHECK_THROW(initialise_companion_matrix(-1),std::invalid_argument);
    BOOST_CHECK_THROW(initialise_companion_matrix(0),std::invalid_argument);
    BOOST_CHECK_THROW(initialise_companion_matrix(1),std::invalid_argument);

    for (int degree=2;degree<12;++degree)
    {
      double **mat=initialise_companion_matrix(degree);

      for (int itr=0;itr<degree;++itr)
      {
        for (int jtr=0;jtr<degree;++jtr)
        { 
          if (itr>0&&jtr==itr-1) { BOOST_CHECK(*(*(mat+itr)+jtr)==1.); }
          else { BOOST_CHECK(*(*(mat+itr)+jtr)==0.); }
        }
      }

      for (int itr=0;itr<degree;++itr) { delete[] mat[itr]; }
      delete[] mat;
    }
    
    for (int degree=2;degree<12;++degree)
    {
      double **mat=initialise_companion_matrix(degree);

      double *coeffs=new double[degree];
      std::random_device rd;
      std::mt19937_64 gen(rd());
      std::uniform_real_distribution<double> dist(-10,10);
      std::generate(coeffs,coeffs+degree,[dist,gen]() mutable { return dist(gen); });

      assign_companion_matrix(mat,coeffs,degree);
      for (int itr=0;itr<degree;++itr) { BOOST_CHECK(*(*(mat+itr)+degree-1)==*(coeffs+itr)); }

      delete[] coeffs;
      for (int itr=0;itr<degree;++itr) { delete[] mat[itr]; }
      delete[] mat;
    }
  }
BOOST_AUTO_TEST_SUITE_END()
} // namespace LinalgTesting
