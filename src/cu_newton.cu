
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// cu_fractal.hpp                                                                                                                        //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - April, 2020                                                                                               //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// CUDA accelerated implementation for producing Newton's fractals


#include <cu_fractal.hpp>

__device__ thrust::pair<cuDoubleComplex,cuDoubleComplex> polynomial_and_deriv(const cuDoubleComplex &x, const double * const coeffs
  , const int degree)
{
  cuDoubleComplex p,p_prime;
  p=make_cuDoubleComplex(*(coeffs+degree),0.);
  p_prime=make_cuDoubleComplex(0.,0.);

  for (int itr=degree-1;itr>=0;--itr)
  {
    p_prime=cuCadd(cuCmul(x,p_prime),p);
    p=cuCadd(cuCmul(x,p),make_cuDoubleComplex(*(coeffs+itr),0.0));
  }

  return thrust::make_pair(p,p_prime);
}

__device__ cuDoubleComplex newton_root(const double * const coeffs, int * const itr_taken, cuDoubleComplex x, const int degree
  , const int max_itr, const double tol)
{
  for (int itr=0;itr<max_itr;++itr)
  {
    // get the current function value and derivative
    cuDoubleComplex f_x,g_x;
    thrust::tie(f_x,g_x)=polynomial_and_deriv(x,coeffs,degree);
    // converged to a root
    if (cuCabs(f_x)<tol)
    {
      *itr_taken=itr;
      return x;
    }
    // derivative is flat and we can't update
    if (cuCreal(g_x)==0.&&cuCimag(g_x)==0.)
    {
      *itr_taken=NPP_MAX_32S ;
      return make_cuDoubleComplex(CUDART_NAN,CUDART_NAN);
    }
    // update
    x=cuCsub(x,cuCdiv(f_x,g_x));
  }
  // couldn't find a root in the given number of iterations
  *itr_taken=NPP_MAX_32S ;
  return make_cuDoubleComplex(CUDART_NAN,CUDART_NAN);
}
