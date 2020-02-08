
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
    const int degree=2;
    double **mat=new double*[degree];
    for (int itr=0;itr<degree;++itr) { *(mat+itr)=new double[degree]; }

    BOOST_CHECK_THROW(initialise_companion_matrix(mat,-1),std::invalid_argument);
    BOOST_CHECK_THROW(initialise_companion_matrix(mat,0),std::invalid_argument);
    BOOST_CHECK_THROW(initialise_companion_matrix(mat,1),std::invalid_argument);

    for (int itr=0;itr<degree;++itr) { delete[] *(mat+itr); }
    delete[] mat;

    for (int degree=2;degree<12;++degree)
    {
      double **mat=new double*[degree];
      for (int itr=0;itr<degree;++itr) { *(mat+itr)=new double[degree]; }
      initialise_companion_matrix(mat,degree);

      for (int itr=0;itr<degree;++itr)
      {
        for (int jtr=0;jtr<degree;++jtr)
        { 
          if (itr>0&&jtr==itr-1) { BOOST_CHECK(*(*(mat+itr)+jtr)==1.); }
          else { BOOST_CHECK(*(*(mat+itr)+jtr)==0.); }
        }
      }

      for (int itr=0;itr<degree;++itr) { delete[] *(mat+itr); }
      delete[] mat;
    }
    
    for (int degree=2;degree<12;++degree)
    {
      double **mat=new double*[degree];
      for (int itr=0;itr<degree;++itr) { *(mat+itr)=new double[degree]; }
      initialise_companion_matrix(mat,degree);

      double *coeffs=new double[degree];
      std::random_device rd;
      std::mt19937_64 gen(rd());
      std::uniform_real_distribution<double> dist(-10,10);
      std::generate(coeffs,coeffs+degree,[dist,gen]() mutable { return dist(gen); });

      assign_companion_matrix(mat,coeffs,degree);
      for (int itr=0;itr<degree;++itr) { BOOST_CHECK(*(*(mat+itr)+degree-1)==*(coeffs+itr)); }

      delete[] coeffs;
      for (int itr=0;itr<degree;++itr) { delete[] *(mat+itr); }
      delete[] mat;
    }
  }

  BOOST_AUTO_TEST_CASE(dot_product)
  {
    // variable setup
    const int degree=3;
    double **mat=new double*[degree];
    for (int itr=0;itr<degree;++itr) { *(mat+itr)=new double[degree]; }
    double *vec_one=new double[degree];
    double *vec_two=new double[degree];

    // make a unit matrix
    for (int itr=0;itr<degree;++itr) { std::fill(*(mat+itr),*(mat+itr)+degree,0); }
    **mat=1;
    *(*(mat+1)+1)=1;
    *(*(mat+2)+2)=1;

    std::iota(vec_one,vec_one+degree,0);
    dot(mat,vec_one,vec_two,degree);

    BOOST_CHECK(*(vec_two+0)==0);
    BOOST_CHECK(*(vec_two+1)==1);
    BOOST_CHECK(*(vec_two+2)==2);

    *(*(mat)+2)=-2;
    *(*(mat+2))=-2;
    dot(mat,vec_one,vec_two,degree);

    BOOST_CHECK(*(vec_two+0)==-4);
    BOOST_CHECK(*(vec_two+1)==1);
    BOOST_CHECK(*(vec_two+2)==2);

    for (int itr=0;itr<degree;++itr) { std::iota(*(mat+itr),*(mat+itr)+degree,itr); }

    for (int itr=0;itr<12;++itr)
    {
      std::fill(vec_one,vec_one+degree,itr);
      std::fill(vec_two,vec_two+degree,itr);
      BOOST_CHECK(dot<double>(vec_one,vec_two,degree)==3*itr*itr);

      dot(mat,vec_one,vec_two,degree);
      const int init=3*itr,step=init;
      BOOST_CHECK(*(vec_two+0)==init);
      BOOST_CHECK(*(vec_two+1)==init+step);
      BOOST_CHECK(*(vec_two+2)==init+2*step);
    }

    delete[] vec_one;
    delete[] vec_two;
    for (int itr=0;itr<degree;++itr) { delete[] *(mat+itr); }
    delete[] mat;
  }

  BOOST_AUTO_TEST_CASE(eigen)
  {
    int len=5;
    eigenpair<double> pair_one(len);

    BOOST_CHECK(pair_one()==0.);
    BOOST_CHECK(double(pair_one)==0.);

    for (int itr=0;itr<len;++itr) { BOOST_CHECK(pair_one[itr]==0.); }

    len=3;
    eigenpair<double> pair_two(3.14,len);
    BOOST_CHECK(pair_two()==3.14);
    for (int itr=0;itr<len;++itr) { BOOST_CHECK(pair_two[itr]==0.); }

    len=7;
    double *vec=new double[len];
    std::iota(vec,vec+len,0);

    eigenpair<double> pair_three(vec,1.618,len);
    BOOST_CHECK(pair_three()==1.618);
    for (int itr=0;itr<len;++itr) { BOOST_CHECK(pair_three[itr]==*(vec+itr)); }
  }
BOOST_AUTO_TEST_SUITE_END()
} // namespace LinalgTesting
