
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// linalg.cpp                                                                                                                            //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - February, 2020                                                                                            //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// Linear algebra routines for fractals


#include <fractal.hpp>

template <typename T>
eigenpair<T>::eigenpair(const size_t size) : size(size), vector(new T[size]), value(T())
{
}

template <typename T>
eigenpair<T>::eigenpair(const T value, const size_t size) : eigenpair(size), value(value)
{
}

template <typename T>
eigenpair<T>::eigenpair(const eigenpair<T> &other) : eigenpair(other.value,other.size)
{
  std::copy(other.vector,other.vector+other.size,this->vector);
}

template <typename T>
eigenpair<T>::eigenpair(eigenpair<T> &&other) : eigenpair(other.value,other.size)
{
  swap(*this,other)
}

template <typename T>
eigenpair<T>::~eigenpair()
{
  delete[] vector;
}

template <typename T>
T eigenpair<T>::norm(void)
{
  return std::sqrt(this->squared_norm());
}

template <typename T>
T eigenpair<T>::squared_norm(void)
{
  return dot(this->vector,this->vector,this->size);
}

template <typename T>
void eigenpair<T>::normalise(void)
{
  T norm=this->norm();
  for (int itr=0;itr<this->size;++itr) { *(this->vector+itr)/=norm; }
}

template <typename T>
eigenpair<T> &eigenpair<T>::operator=(eigenpair<T> other)
{
  swap(*this,other);
  return *this;
}

template <typename T>
[[nodiscard]] T &eigenpair<T>::operator[](size_t idx)
{
  if (idx>=size)
  {
    std::ostringstream err_message;
    err_message << idx << " is not in range for an eigenvector of size " << this->size;
    throw std::out_of_range(err_message.str());
  }
  return *(vector+idx)
}

template <typename T>
[[nodiscard]] const T &eigenpair<T>::operator[](const size_t idx)
{
  return (*this)[const_cast<size_t>(idx)]
}

template <typename T>
eigenpair<T>::operator scalar_type() const noexcept
{
  return this->value;
}

template <typename T>
[[nodiscard]] typename eigenpair<T>::scalar_type &eigenpair<T>::operator()()
{
  return this->value;
}

template <typename T>
[[nodiscard]] typename eigenpair<T>::scalar_type eigenpair<T>::operator()() const noexcept
{
  return this->value;
}

template <typename T>
void swap(eigenpair<T> &first, eigenpair<T> &second) noexcept
{
  std::swap(first.value,second.value);
  std::swap(first.vector,second.vector);
  std::swap(first.size,second.size);
}


void initialise_companion_matrix(double * const * const mat, const int degree)
{
  if (degree<2) { throw std::invalid_argument("`degree` must be greater than or equal to 2"); }
  
  for (int itr=0;itr<degree;++itr) { std::fill(*(mat+itr),*(mat+itr)+degree,0.); }
  for (int itr=1,jtr=0;itr<degree,jtr<degree-1;++itr,++jtr) { *(*(mat+itr)+jtr)=1.; }
}

void assign_companion_matrix(double * const * const mat, double const * const coeffs, const int degree)
{
  for (int itr=0;itr<degree;++itr) { *(*(mat+itr)+degree-1)=*(coeffs+itr); }
}

void __declspec(dllexport) roots(double const * const coeffs, double * const roots_re, double * const roots_im, const int degree)
{
}
