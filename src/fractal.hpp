
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// fractal.hpp                                                                                                                           //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - January, 2020                                                                                             //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// Header file for Fractal project


#pragma once

#ifndef FRACTAL_H__
#define FRACTAL_H__

#include <cmath>

#include <chrono>
#include <complex>
#include <iostream>
#include <limits>
#include <thread>
#include <tuple>
#include <vector>

class convergence_error : public std::exception {};

template <typename T>
struct eigenpair
{
  using scalar_type=T;

  eigenpair(const size_t size) : size(size), vector(new T[size]), value(T())
  {
    std::fill(this->vector,this->vector+this->size,T());
  }
  eigenpair(const T value, const size_t size) : eigenpair(size)
  {
    this->value=value;
  }
  eigenpair(const T * const vector, const T value, const size_t size) : eigenpair(value,size)
  {
    std::copy(vector,vector+size,this->vector);
  }

  eigenpair(const eigenpair &other) : eigenpair(other.value,other.size)
  {
    std::copy(other.vector,other.vector+other.size,this->vector);
  }
  eigenpair(eigenpair &&other) : eigenpair(other.value,other.size)
  {
    swap(*this,other);
  }

  ~eigenpair()
  {
    delete[] vector;
  }

  T norm(void)
  {
    return std::sqrt(this->squared_norm());
  }
  T squared_norm(void)
  {
    return dot<T,T,T>(this->vector,this->vector,this->size);
  }
  void normalise(void)
  {
    T norm=this->norm();
    for (int itr=0;itr<this->size;++itr) { *(this->vector+itr)/=norm; }
  }

  // assignment
  eigenpair &operator=(eigenpair other)
  {
    swap(*this,other);
    return *this;
  }
  // getters and setters
  [[nodiscard]] T &operator[](size_t idx)
  {
    if (idx>=size)
    {
      std::ostringstream err_message;
      err_message << idx << " is not in range for an eigenvector of size " << this->size;
      throw std::out_of_range(err_message.str());
    }
    return *(vector+idx);
  }
  operator scalar_type() const noexcept
  {
    return this->value;
  }
  [[nodiscard]] scalar_type &operator()()
  {
    return this->value;
  }
  [[nodiscard]] scalar_type operator()() const noexcept
  {
    return this->value;
  }
  // indirection
  const scalar_type * const operator*()
  {
    return vector;
  }

  // swap
  template<typename U>
  friend inline void swap(eigenpair<U> &first, eigenpair<U> &second) noexcept;

private:
  T value;
  T *vector;
  size_t size;
};

template <typename T>
inline void swap(eigenpair<T> &first, eigenpair<T> &second) noexcept
{
  std::swap(first.value,second.value);
  std::swap(first.vector,second.vector);
  std::swap(first.size,second.size);
}


template <typename T>
inline void initialise_companion_matrix(T * const * const mat, const int degree)
{
  if (degree<2) { throw std::invalid_argument("`degree` must be greater than or equal to 2"); }
  
  for (int itr=0;itr<degree;++itr) { std::fill(*(mat+itr),*(mat+itr)+degree,0.); }
  for (int itr=1,jtr=0;itr<degree,jtr<degree-1;++itr,++jtr) { *(*(mat+itr)+jtr)=1.; }
}

template <typename T, typename U>
inline void assign_companion_matrix(T * const * const mat, const U * const coeffs, const int degree)
{
  for (int itr=0;itr<degree;++itr) { *(*(mat+itr)+degree-1)=*(coeffs+itr); }
}


template <typename VecT, typename VecU, typename Scalar>
inline Scalar dot(const VecT * const vector_one, const VecU * const vector_two, const size_t size)
{
  Scalar dot_product=Scalar();
  for (int itr=0;itr<size;++itr) { dot_product+=*(vector_one+itr)**(vector_two+itr); }
  return dot_product;
}

template <typename Mat, typename VecIn, typename VecOut>
inline void dot(const Mat * const * const mat, const VecIn * const vector, VecOut * const out, const size_t size)
{
  for (int itr=0;itr<size;++itr) { *(out+itr)=dot<Mat,VecIn,VecOut>(*(mat+itr),vector,size); }
}


inline std::vector<int> iteration_limits(const int num_threads, const int yresolution)
{
  const int span=yresolution/num_threads;
  std::vector<int> increments;

  for (int itr=0;itr<num_threads;++itr) { increments.push_back(span*itr); }
  increments.push_back(yresolution);

  return increments;
}

int iterate(std::complex<double> x, const std::complex<double> &c, const int max_itr);
void compute_mandelbrot_range(int **iterations, const int max_itr, const int xresolution, const int start_itr, const int end_itr
  , const double startx, const double starty, const double deltax, const double deltay, const int total, bool verbose);
void __declspec(dllexport) sample_mandelbrot(int **iterations, const int max_itr, const int num_threads, const int xresolution
  , const int yresolution, int * const limit, const double startx, const double endx, const double starty, const double endy
  , const bool verbose);

std::pair<std::complex<double>,std::complex<double>> polynomial_and_deriv(const std::complex<double> &x, double * const coeffs
  , const int degree);
std::complex<double> newton_root(double * const coeffs, int * const itr_taken, const std::complex<double> x, const int degree
  , const int max_itr, const double tol);
void compute_newton_range(double **re, double **im, int **iterations, double * coeffs, const int max_itr, const int degree
  , const int xresolution, const int start_itr, const int end_itr, const double startx, const double starty, const double deltax
  , const double deltay, const int total, bool verbose);
void __declspec(dllexport) sample_newton(double **re, double **im, int **iterations, double * coeffs, const int max_itr
  , const int num_threads, const int degree, const int xresolution, const int yresolution, int * const limit, const double startx
  , const double endx, const double starty, const double endy, const bool verbose);

eigenpair<std::complex<double>> power_iteration(std::complex<double> const * const * const mat, std::complex<double> const * const vec
  , const size_t size, const double tol, const int max_itr);
void deflate(std::complex<double> * const * const mat, const eigenpair<std::complex<double>> &pair);
void __declspec(dllexport) roots(double const * const coeffs, double * const roots_re, double * const roots_im, const int degree);

#endif // FRACTAL_H__
