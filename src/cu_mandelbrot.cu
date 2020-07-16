
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// cu_mandelbrot.cu                                                                                                                      //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - April, 2020                                                                                               //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// CUDA accelerated implementation for producing fractals from the Mandelbrot Set


#include <cu_fractal.hpp>

__device__ int iterate(cuDoubleComplex x, const cuDoubleComplex &c, const int max_itr)
{
  //
  // CUDA device to iterate a given number under x^2 + c until its absolute value becomes greater than 2 or the maximum number of 
  // iterations is reached.
  //
  // parameters
  // ----------
  // x : cuDoubleComplex
  //  - the starting value of the iteration
  // c : const cuDoubleComplex &
  //  - the constant added at each iteration
  // max_itr : const int 
  //  - the maximum number of iterations allowed
  //
  // returns
  // -------
  // int
  //  - the number of iterations needed for the aboslute value of x became greater than 2 under iteration
  //    NPP_MAX_32S if the absolute value of x did not become greater in 2 within `max_itr` iterations
  //

  for (int itr=0;itr<max_itr;++itr)
  {
    if (cuCabs(x)>2) { return itr; }
    x=cuCadd(cuCmul(x,x),c);
  }
  return NPP_MAX_32S;
}

__global__ void compute_mandelbrot(int * const d_iterations, const int max_itr, const int xresolution, const int yresolution
  , const double startx, const double starty, const double deltax, const double deltay)
{
  //
  // CUDA kernel to determine if a given number of the complex plane is contained within the Mandelbrot Set through iteration
  //
  // parameters
  // ----------
  // d_iterations : int * const
  //  - 1D flat array representing a 2D array to write out either the number of iterations needed reach a root or a marker that this root 
  //      could not be reached
  // max_itr : const int
  //  - the maximum number of iterations allowed
  // xresolution,yresolution : const int 
  //  - the number of steps to take in the x-direction (the real components) and the y-direction (the imaginary components)
  // startx,starty : const int 
  //  - the real and imaginary components of the number defining the bottom left corner of the entire space being sampled
  // deltax,deltay : const int 
  //  - the size of the step to take in the x- and y-direction
  //

  // determine where we are in memory
  const int idy=blockIdx.y*blockDim.y+threadIdx.y,idx=blockIdx.x*blockDim.x+threadIdx.x,ind=idy*xresolution+idx;
  // check we haven't gone out of bounds
  if (idx>=xresolution||idy>=yresolution) { return; }

  // determine the current point
  const double imag=starty+deltay*idy,real=startx+deltax*idx;
  // determine the number of iterations
  d_iterations[ind]=iterate(make_cuDoubleComplex(0.,0.),make_cuDoubleComplex(real,imag),max_itr);
}

__global__ void compute_julia(int * const d_iterations, const double re, const double im, const int max_itr
  , const int xresolution, const int yresolution, const double startx, const double starty, const double deltax, const double deltay)
{
  //
  // CUDA kernel to determine if a given number of the complex plane is contained within the Julia Set of a given complex number through
  // iteration
  //
  // parameters
  // ----------
  // d_iterations : int * const
  //  - 1D flat array representing a 2D array to write out either the number of iterations needed reach a root or a marker that this root 
  //      could not be reached
  // re,im : const double
  //  - the real and imaginary parts of the complex number to find the Julia Set of
  // max_itr : const int
  //  - the maximum number of iterations allowed
  // xresolution,yresolution : const int 
  //  - the number of steps to take in the x-direction (the real components) and the y-direction (the imaginary components)
  // startx,starty : const int 
  //  - the real and imaginary components of the number defining the bottom left corner of the entire space being sampled
  // deltax,deltay : const int 
  //  - the size of the step to take in the x- and y-direction
  //

  // determine where we are in memory
  const int idy=blockIdx.y*blockDim.y+threadIdx.y,idx=blockIdx.x*blockDim.x+threadIdx.x,ind=idy*xresolution+idx;
  // check we haven't gone out of bounds
  if (idx>=xresolution||idy>=yresolution) { return; }

  // determine the current point
  const double imag=starty+deltay*idy,real=startx+deltax*idx;
  // determine the number of iterations
  d_iterations[ind]=iterate(make_cuDoubleComplex(real,imag),make_cuDoubleComplex(re,im),max_itr);
}

int sample_mandelbrot(int * const h_itr, const int max_itr, const int xresolution, const int yresolution
  , const double startx, const double endx, const double starty, const double endy, const bool verbose)
{
  //
  // Determine if numbers in a given subset of the complex plane are contained within the Mandelbrot Set through iteration, with CUDA
  // acceleration
  //
  // parameters
  // ----------
  // h_itr : int * const
  //  - 1D flat array representing a 2D array to write out either the number of iterations needed to diverge or a marker that the number was
  //      stable
  // max_itr : const int
  //  - the maximum number of iterations allowed
  // xresolution,yresolution : const int 
  //  - the number of steps to take in the x- and y-direction (the real and imaginary components)
  // startx,endx,starty,endy : const int
  //  - the first and last values to sample at 
  // verbose : bool
  //  - flag to control logging to console
  //
  // returns
  // -------
  // int
  //  - the value which marks that a starting point could not be shown to not be in the Mandelbrot Set
  //

  // computation parameters
  const double deltax=(endx-startx)/xresolution,deltay=(endy-starty)/yresolution;
  const int total=xresolution*yresolution;
  // memory parameters
  const int i_size=total*sizeof(int);

  // device memory pointers
  int *d_itr=nullptr;
  // allocate device memory
  CUDA_REQUIRE_SUCCESS(cudaMalloc(reinterpret_cast<void **>(&d_itr),static_cast<size_t>(i_size)));

  // GPU memory setup
  const dim3 dim_block(32,32),dim_grid((xresolution+dim_block.x-1)/dim_block.x,(yresolution+dim_block.y-1)/dim_block.y);
  // run and time
  float elapsed;
  cudaEvent_t start,stop;

  CUDA_REQUIRE_SUCCESS(cudaEventCreate(&start));
  CUDA_REQUIRE_SUCCESS(cudaEventCreate(&stop));
  CUDA_REQUIRE_SUCCESS(cudaEventRecord(start,0));

  compute_mandelbrot<<<dim_grid,dim_block>>>(d_itr,max_itr,xresolution,yresolution,startx,starty,deltax,deltay);
  // check for errors
  CUDA_REQUIRE_SUCCESS(cudaPeekAtLastError());
  CUDA_REQUIRE_SUCCESS(cudaDeviceSynchronize());

  CUDA_REQUIRE_SUCCESS(cudaEventRecord(stop,0));
  CUDA_REQUIRE_SUCCESS(cudaEventSynchronize(stop));
  CUDA_REQUIRE_SUCCESS(cudaEventElapsedTime(&elapsed,start,stop));

  if (verbose)
  {
    std::cout << total << " points processed." << '\n'
      << "Time taken: " << elapsed/1000 << "s." << std::endl;
  }

  // copy back to host
  CUDA_REQUIRE_SUCCESS(cudaMemcpy(h_itr,d_itr,static_cast<size_t>(i_size),cudaMemcpyDeviceToHost));

  // free GPU memory
  CUDA_REQUIRE_SUCCESS(cudaFree(d_itr));

  return NPP_MAX_32S;
}

int sample_julia(int * const h_itr, const double re, const double im, const int max_itr
  , const int xresolution, const int yresolution, const double startx, const double endx, const double starty, const double endy
  , const bool verbose)
{
  //
  // Determine if numbers in a given subset of the complex plane are contained within the Julia Set of a given complex number through 
  // iteration, with CUDA acceleration
  //
  // parameters
  // ----------
  // h_itr : int * const
  //  - 1D flat array representing a 2D array to write out either the number of iterations needed to diverge or a marker that the number was
  //      stable
  // re,im : const double
  //  - the real and imaginary parts of the complex number to find the Julia Set of
  // max_itr : const int
  //  - the maximum number of iterations allowed
  // xresolution,yresolution : const int 
  //  - the number of steps to take in the x- and y-direction (the real and imaginary components)
  // startx,endx,starty,endy : const int
  //  - the first and last values to sample at 
  // verbose : bool
  //  - flag to control logging to console
  //
  // returns
  // -------
  // int
  //  - the value which marks that a starting point could not be shown to not be in the Mandelbrot Set
  //

  // computation parameters
  const double deltax=(endx-startx)/xresolution,deltay=(endy-starty)/yresolution;
  const int total=xresolution*yresolution;
  // memory parameters
  const int i_size=total*sizeof(int);

  // device memory pointers
  int *d_itr=nullptr;
  // allocate device memory
  CUDA_REQUIRE_SUCCESS(cudaMalloc(reinterpret_cast<void **>(&d_itr),static_cast<size_t>(i_size)));

  // GPU memory setup
  const dim3 dim_block(32,32),dim_grid((xresolution+dim_block.x-1)/dim_block.x,(yresolution+dim_block.y-1)/dim_block.y);
  // run and time
  float elapsed;
  cudaEvent_t start,stop;

  CUDA_REQUIRE_SUCCESS(cudaEventCreate(&start));
  CUDA_REQUIRE_SUCCESS(cudaEventCreate(&stop));
  CUDA_REQUIRE_SUCCESS(cudaEventRecord(start,0));

  compute_julia<<<dim_grid,dim_block>>>(d_itr,re,im,max_itr,xresolution,yresolution,startx,starty,deltax,deltay);
  // check for errors
  CUDA_REQUIRE_SUCCESS(cudaPeekAtLastError());
  CUDA_REQUIRE_SUCCESS(cudaDeviceSynchronize());

  CUDA_REQUIRE_SUCCESS(cudaEventRecord(stop,0));
  CUDA_REQUIRE_SUCCESS(cudaEventSynchronize(stop));
  CUDA_REQUIRE_SUCCESS(cudaEventElapsedTime(&elapsed,start,stop));

  if (verbose)
  {
    std::cout << total << " points processed." << '\n'
      << "Time taken: " << elapsed/1000 << "s." << std::endl;
  }

  // copy back to host
  CUDA_REQUIRE_SUCCESS(cudaMemcpy(h_itr,d_itr,static_cast<size_t>(i_size),cudaMemcpyDeviceToHost));

  // free GPU memory
  CUDA_REQUIRE_SUCCESS(cudaFree(d_itr));

  return NPP_MAX_32S;
}
