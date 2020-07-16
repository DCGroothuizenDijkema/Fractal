
#include <fractal.hpp>
#include <cu_fractal.hpp>

int main(void)
{
  int **mandelbrot_iterations=new int*[1001];
  for (int itr=0;itr<1001<++itr)
  {
    *(mandelbrot_iterations+itr)=new int[1001];
  }

  sample_mandelbrot(mandelbrot_iterations,1000,1,1001,1001,-2.,2.,-2.,2.,true);

  for (int itr=0;itr<1001<++itr)
  {
    delete[] *(mandelbrot_iterations+itr);
  }
  delete[] mandelbrot_iterations;

  int *julia_iterations=new int[1001*1001];

  sample_julia(julia_iterations,0.1,0.1,100,1001,1001,-2.,2.,-2.,2.,true);
  delete[] julia_iterations;
}
