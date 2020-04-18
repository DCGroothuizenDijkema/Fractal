
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
  // determine where we are in memory
  const int idy=blockIdx.y*blockDim.y+threadIdx.y,idx=blockIdx.x*blockDim.x+threadIdx.x,ind=idy*xresolution+idx;
  // check we haven't gone out of bounds
  if (idx>=xresolution||idy>=yresolution) { return; }

  // determine the current point
  const double imag=starty+deltay*idy,real=startx+deltax*idx;
  // determine the number of iterations
  d_iterations[ind]=iterate(make_cuDoubleComplex(0.,0.),make_cuDoubleComplex(real,imag),max_itr);
}

__global__ void compute_julia(int * const d_iterations, const cuDoubleComplex &c, const int max_itr
  , const int xresolution, const int yresolution, const double startx, const double starty, const double deltax, const double deltay)
{
  // determine where we are in memory
  const int idy=blockIdx.y*blockDim.y+threadIdx.y,idx=blockIdx.x*blockDim.x+threadIdx.x,ind=idy*xresolution+idx;
  // check we haven't gone out of bounds
  if (idx>=xresolution||idy>=yresolution) { return; }

  // determine the current point
  const double imag=starty+deltay*idy,real=startx+deltax*idx;
  // determine the number of iterations
  d_iterations[ind]=iterate(make_cuDoubleComplex(real,imag),c,max_itr);
}

int __declspec(dllexport) sample_mandelbrot(int * const h_itr, const int max_itr, const int xresolution, const int yresolution
  , const double startx, const double endx, const double starty, const double endy, const bool verbose)
{
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
    std::cout << total << " points processed." << std::endl
      << "Time taken: " << elapsed/1000 << "s." << std::endl;
  }

  // copy back to host
  CUDA_REQUIRE_SUCCESS(cudaMemcpy(h_itr,d_itr,static_cast<size_t>(i_size),cudaMemcpyDeviceToHost));

  // free GPU memory
  CUDA_REQUIRE_SUCCESS(cudaFree(d_itr));

  return NPP_MAX_32S;
}

int __declspec(dllexport) sample_julia(int * const h_itr, const int max_itr, const int xresolution, const int yresolution
  , const double startx, const double endx, const double starty, const double endy, const bool verbose)
{
  return NPP_MAX_32S;
}
