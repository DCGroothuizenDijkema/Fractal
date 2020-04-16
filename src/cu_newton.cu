
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
      return make_cuDoubleComplex(CUDART_INF,CUDART_INF);
    }
    // update
    x=cuCsub(x,cuCdiv(f_x,g_x));
  }
  // couldn't find a root in the given number of iterations
  *itr_taken=NPP_MAX_32S ;
  return make_cuDoubleComplex(CUDART_INF,CUDART_INF);
}

__global__ void compute_newton(double *d_re, double *d_im, int *d_itr, double * const coeffs, const int max_itr, const int degree
  , const int xresolution, const int yresolution, const double startx, const double starty, const double deltax, const double deltay)
{
  const int idx=blockIdx.x*blockDim.x+threadIdx.x,idy=threadIdx.y,offset=(idx*xresolution)+idy;
  // thread check
  if (idx>=xresolution||idy>=yresolution) { return; }
  const double imag=starty+deltay*idy,real=startx+deltax*idx;

  // determine the root reached and the number of iterations to get there
  cuDoubleComplex root=newton_root(coeffs,(d_itr+offset),make_cuDoubleComplex(real,imag),degree,max_itr,1e-6);
  *(d_re+offset)=cuCreal(root);
  *(d_im+offset)=cuCimag(root);
}

void __declspec(dllexport) sample_newton(double *h_re, double *h_im, int *h_itr, double *coeffs, const int max_itr, const int degree
  , const int xresolution, const int yresolution, const double startx, const double endx, const double starty, const double endy)
{
  double *d_re=nullptr,*d_im=nullptr;
  int *d_itr=nullptr;

  const double deltax=(endx-startx)/xresolution,deltay=(endy-starty)/yresolution;
  const int total=xresolution*yresolution,d_size=total*sizeof(double),i_size=total*sizeof(int);

  cudaMalloc(reinterpret_cast<void **>(&d_re),static_cast<size_t>(d_size));
  cudaMalloc(reinterpret_cast<void **>(&d_im),static_cast<size_t>(d_size));
  cudaMalloc(reinterpret_cast<void **>(&d_itr),static_cast<size_t>(i_size));

  dim3 dim_block(32,32),dim_grid((xresolution*yresolution)/(dim_block.x*dim_block.y),1,1);

  compute_newton<<<dim_grid,dim_block>>>(d_re,d_im,d_itr,coeffs,max_itr,degree,xresolution,yresolution,startx,starty,deltax,deltay);

  cudaMemcpy(h_re,d_re,static_cast<size_t>(total),cudaMemcpyDeviceToHost);
  cudaMemcpy(h_im,d_im,static_cast<size_t>(total),cudaMemcpyDeviceToHost);
  cudaMemcpy(h_itr,d_itr,static_cast<size_t>(total),cudaMemcpyDeviceToHost);

  cudaFree(d_re);
  cudaFree(d_im);
  cudaFree(d_itr);
}

void __declspec(dllexport) assign_roots(int *index, double *re, double *im, const double * const roots_re, const double * const roots_im
  , const int degree, const int xresolution, const int yresolution)
{
  // get a list of all roots, formed from the input vectors giving the real and imaginary components of the roots
  std::vector<std::complex<double>> roots;
  zip(roots_re,roots_re+degree,roots_im,roots_im+degree,std::back_inserter(roots));
  
  for (int itr=0;itr<xresolution*yresolution;++itr)
  {
    // if the current value is marked with infinity, no root was reached from it and its index is 0
    if (*(re+itr)==std::numeric_limits<double>::infinity()||*(im+itr)==std::numeric_limits<double>::infinity())
    {
      *(index+itr)=-1;
      continue;
    }

    std::complex<double> val(*(re+itr),*(im+itr));
    // determine the difference between the current value and all roots
    std::vector<std::complex<double>> diffs;
    std::transform(std::begin(roots),std::end(roots),std::back_inserter(diffs)
      ,[val](std::complex<double> root) { return abs(root-val); });
    // find the argmin of the differences, and, therefore, which root was converged to
    *(index+itr)=static_cast<int>(argmin(std::cbegin(diffs),std::cend(diffs)
      ,[](const std::complex<double> &x, const std::complex<double> &y){ return abs(x)<abs(y); }
    ));
  }
}
