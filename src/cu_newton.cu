
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// cu_newton.cu                                                                                                                          //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - April, 2020                                                                                               //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// CUDA accelerated implementation for producing Newton's fractals


#include <cu_fractal.hpp>

__device__ thrust::pair<cuDoubleComplex,cuDoubleComplex> polynomial_and_deriv(const cuDoubleComplex &x, const double * const coeffs
  , const int degree)
{
  //
  // CUDA device to evaluate the value and derivative of a polynomial at a point using Horner's method
  //
  // parameters
  // ----------
  // x : const cuDoubleComplex &
  //  - the point to evaluate at
  // coeffs : const double * const
  //  - the coefficients of the polynomial given in order of the lowest degree to highest
  //    if the polynomial is written p(x)=a_n*x^(n-1)+...+a_k*x^(n-k-1)+...+a_1*x+a_0, then coeffs should have the following form:
  //      *(coeffs+0)==a_0
  //      *(coeffs+1)==a_1
  //      *(coeffs+k)==a_k
  //      *(coeffs+n)==a_n
  // degree : const int
  //  - the degree of the polynomial
  //
  // returns
  // -------
  // thrust::pair<cuDoubleComplex,cuDoubleComplex>
  //  - the function value and derivative of the given polynomial evaluated at the given point
  //
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
  //
  // CUDA device to apply the Newton-Raphson method to a given number to find the root of a given polynomial
  //
  // parameters
  // ----------
  // coeffs : const double * const
  //  - the coefficients of the polynomial given in order of the lowest degree to highest
  //    see polynomial_and_deriv() for requirements
  // itr_taken : int * const 
  //  - a poitner to write out the number of iterations needed to reach a root
  // x : cuDoubleComplex
  //  - the number to a start from
  // degree : const int
  //  - the degree of the polynomial
  // max_itr : const int
  //  - the maximum number of iterations allowed
  // tol : const double
  //  - the tolerance within which a root will be deemed to have been reached
  //
  // returns
  // -------
  // cuDoubleComplex
  //  - the root reached from the given number
  //

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
      *itr_taken=NPP_MAX_32S;
      return make_cuDoubleComplex(CUDART_INF,CUDART_INF);
    }
    // update
    x=cuCsub(x,cuCdiv(f_x,g_x));
  }
  // couldn't find a root in the given number of iterations
  *itr_taken=NPP_MAX_32S ;
  return make_cuDoubleComplex(CUDART_INF,CUDART_INF);
}

__global__ void compute_newton(double * const d_re, double * const d_im, int * const d_itr, const double * const d_coeffs, const int max_itr
  , const int degree, const int xresolution, const int yresolution, const double startx, const double starty, const double deltax
  , const double deltay)
{
  //
  // CUDA kernel to find the roots of a polynomial a given number in the complex plane converges to with Newton's method
  //
  // parameters
  // ----------
  // d_re,d_im : double * const
  //  - 1D flat arrays representing 2D arrays to write out the root which the numbers in the subset cnoverged to
  // d_itr : int * const
  //  - 1D flat array representing a 2D array to write out either the number of iterations needed reach a root or a marker that this root 
  //      could not be reached
  // d_coeffs : const double * const
  //  - the coefficients of the polynomial given in order of the lowest degree to highest
  //    see polynomial_and_deriv() for requirements
  // max_itr : const int
  //  - the maximum number of iterations allowed
  // degree : const int
  //  - the degree of the polynomial
  // xresolution,yresolution : const int 
  //  - the number of steps to take in the x-direction (the real components) and the y-direction (the imaginary components)
  // start_itr,end_itr : const int 
  //  - the subset of the steps taken from the starting point in the y-direction (the complex component)
  // startx,starty : const double 
  //  - the real and imaginary components of the number defining the bottom left corner of the entire space being sampled
  // deltax,deltay : const double 
  //  - the size of the step to take in the x- and y-direction
  //

  // determine where we are in memory
  const int idy=blockIdx.y*blockDim.y+threadIdx.y,idx=blockIdx.x*blockDim.x+threadIdx.x,ind=idy*xresolution+idx;
  // check we haven't gone out of bounds
  if (idx>=xresolution||idy>=yresolution) { return; }

  // determine the current point
  const double imag=starty+deltay*idy,real=startx+deltax*idx;
  // find the root we converge to from this point
  cuDoubleComplex root=newton_root(d_coeffs,(d_itr+ind),make_cuDoubleComplex(real,imag),degree,max_itr,1e-6);

  // write out
  d_re[ind]=cuCreal(root);
  d_im[ind]=cuCimag(root);
}

int __declspec(dllexport) sample_newton(double * const h_re, double * const h_im, int * const h_itr, const double * const h_coeffs
  , const int max_itr, const int degree, const int xresolution, const int yresolution, const double startx, const double endx
  , const double starty, const double endy, const bool verbose)
{
  //
  // Determine the roots of a polynomial the numbers in a given subset of the complex plane converge to with Newton's method, with CUDA
  // acceleration
  //
  // parameters
  // ----------
  // h_re,h_im : double * const
  //  - 1D flat arrays representing 2D arrays to write out the root which the numbers in the subset cnoverged to
  // h_itr : int * const
  //  - 1D flat array representing a 2D array to write out either the number of iterations needed reach a root or a marker that this root 
  //      could not be reached
  // h_coeffs : const double * const
  //  - the coefficients of the polynomial given in order of the lowest degree to highest
  //    see polynomial_and_deriv() for requirements
  // max_itr : const int
  //  - the maximum number of iterations allowed
  // degree : const int
  //  - the degree of the polynomial
  // xresolution,yresolution : const int 
  //  - the number of steps to take in the x- and y-direction (the real and imaginary components)
  // startx,endx,starty,endy : const double
  //  - the first and last values to sample at 
  // verbose : bool
  //  - flag to control logging to console
  //
  // returns
  // -------
  // int
  //  - the value which marks that a starting point could not be shown to not be in the mandelbrot set
  //

  // computation parameters
  const double deltax=(endx-startx)/xresolution,deltay=(endy-starty)/yresolution;
  const int total=xresolution*yresolution;
  // memory parameters
  const size_t c_size=(degree+1)*sizeof(double),d_size=total*sizeof(double),i_size=total*sizeof(int);

  // device memory pointers
  double *d_re=nullptr,*d_im=nullptr,*d_coeffs=nullptr;
  int *d_itr=nullptr;
  
  // allocate device memory
  CUDA_REQUIRE_SUCCESS(cudaMalloc(reinterpret_cast<void **>(&d_re),d_size));
  CUDA_REQUIRE_SUCCESS(cudaMalloc(reinterpret_cast<void **>(&d_im),d_size));
  CUDA_REQUIRE_SUCCESS(cudaMalloc(reinterpret_cast<void **>(&d_itr),i_size));

  CUDA_REQUIRE_SUCCESS(cudaMalloc(reinterpret_cast<void **>(&d_coeffs),c_size));
  CUDA_REQUIRE_SUCCESS(cudaMemcpy(d_coeffs,h_coeffs,c_size,cudaMemcpyHostToDevice));
  // copy polynomial coefficients over

  // GPU memory setup
  const dim3 dim_block(32,32),dim_grid((xresolution+dim_block.x-1)/dim_block.x,(yresolution+dim_block.y-1)/dim_block.y);
  // run and time
  float elapsed;
  cudaEvent_t start,stop;

  CUDA_REQUIRE_SUCCESS(cudaEventCreate(&start));
  CUDA_REQUIRE_SUCCESS(cudaEventCreate(&stop));
  CUDA_REQUIRE_SUCCESS(cudaEventRecord(start,0));
  
  compute_newton<<<dim_grid,dim_block>>>(d_re,d_im,d_itr,d_coeffs,max_itr,degree,xresolution,yresolution,startx,starty,deltax,deltay);
  // check for errors
  CUDA_REQUIRE_SUCCESS(cudaPeekAtLastError());
  CUDA_REQUIRE_SUCCESS(cudaDeviceSynchronize());

  CUDA_REQUIRE_SUCCESS(cudaEventRecord(stop,0));
  CUDA_REQUIRE_SUCCESS(cudaEventSynchronize(stop));
  CUDA_REQUIRE_SUCCESS(cudaEventElapsedTime(&elapsed,start,stop));

  if (verbose)
  {
    std::cout << total << " points processed." << std::endl
      << "Time taken: " << elapsed/1000 << "s." << std::endl;
  }

  CUDA_REQUIRE_SUCCESS(cudaMemcpy(h_re,d_re,d_size,cudaMemcpyDeviceToHost));
  CUDA_REQUIRE_SUCCESS(cudaMemcpy(h_im,d_im,d_size,cudaMemcpyDeviceToHost));
  CUDA_REQUIRE_SUCCESS(cudaMemcpy(h_itr,d_itr,i_size,cudaMemcpyDeviceToHost));

  // free GPU memory
  CUDA_REQUIRE_SUCCESS(cudaFree(d_re));
  CUDA_REQUIRE_SUCCESS(cudaFree(d_im));
  CUDA_REQUIRE_SUCCESS(cudaFree(d_itr));
  CUDA_REQUIRE_SUCCESS(cudaFree(d_coeffs));

  return NPP_MAX_32S;
}
