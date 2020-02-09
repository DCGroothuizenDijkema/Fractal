
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
    std::complex<double> **mat=new std::complex<double>*[degree];
    for (int itr=0;itr<degree;++itr) { *(mat+itr)=new std::complex<double>[degree]; }

    BOOST_CHECK_THROW(initialise_companion_matrix(mat,-1),std::invalid_argument);
    BOOST_CHECK_THROW(initialise_companion_matrix(mat,0),std::invalid_argument);
    BOOST_CHECK_THROW(initialise_companion_matrix(mat,1),std::invalid_argument);

    for (int itr=0;itr<degree;++itr) { delete[] *(mat+itr); }
    delete[] mat;

    for (int degree=2;degree<12;++degree)
    {
      std::complex<double> **mat=new std::complex<double>*[degree];
      for (int itr=0;itr<degree;++itr) { *(mat+itr)=new std::complex<double>[degree]; }
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
      std::complex<double> **mat=new std::complex<double>*[degree];
      for (int itr=0;itr<degree;++itr) { *(mat+itr)=new std::complex<double>[degree]; }
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
      BOOST_CHECK((dot<double,double,double>(vec_one,vec_two,degree)==3*itr*itr));

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
    double *vec_one=new double[len];
    std::iota(vec_one,vec_one+len,0);
    eigenpair<double> pair_three(vec_one,1.618,len);
    BOOST_CHECK(pair_three()==1.618);
    for (int itr=0;itr<len;++itr) { BOOST_CHECK(pair_three[itr]==*(vec_one+itr)); }

    eigenpair<double> pair_four(pair_three);
    BOOST_CHECK(pair_four()==1.618);
    for (int itr=0;itr<len;++itr) { BOOST_CHECK(pair_four[itr]==*(vec_one+itr)); }

    eigenpair<double> pair_five=pair_three;
    BOOST_CHECK(pair_five()==1.618);
    for (int itr=0;itr<len;++itr) { BOOST_CHECK(pair_five[itr]==*(vec_one+itr)); }

    eigenpair<double> pair_six=std::move(pair_three);
    BOOST_CHECK(pair_six()==1.618);
    for (int itr=0;itr<len;++itr) { BOOST_CHECK(pair_six[itr]==*(vec_one+itr)); }

    pair_six()=2.718;
    BOOST_CHECK(double(pair_six)==2.718);

    pair_six[1]=-11.924;
    BOOST_CHECK(pair_six[1]==-11.924);

    for (int itr=1;itr<12;++itr)
    {
      double *vec_two=new double[itr];
      std::fill(vec_two,vec_two+itr,1);
      eigenpair<double> pair(vec_two,0,itr);

      BOOST_CHECK(pair.norm()==std::sqrt(itr));
      BOOST_CHECK(pair.squared_norm()==itr);
    }

    len=4;
    double *vec_three=new double[len];
    *vec_three=-1;
    *(vec_three+1)=2;
    *(vec_three+2)=5;
    *(vec_three+3)=-2;
    eigenpair<double> pair_seven(vec_three,0.,len);

    BOOST_CHECK(pair_seven.norm()==std::sqrt(34));
    BOOST_CHECK(pair_seven.squared_norm()==34);

    pair_seven.normalise();
    for (int itr=0;itr<len;++itr) { BOOST_CHECK(pair_seven[itr]==*(vec_three+itr)/std::sqrt(34)); }
  }
BOOST_AUTO_TEST_SUITE_END()
} // namespace LinalgTesting
