
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
//                                                                                                                                       //
// cu_fractal.hpp                                                                                                                        //
//                                                                                                                                       //
// D. C. Groothuizen Dijkema - April, 2020                                                                                               //
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

// Header file for common template functions of the Fractal project


#pragma once

#ifndef COMMON_H__
#define COMMON_H__

#include <algorithm>
#include <complex>
#include <limits>
#include <vector>

template <typename InputIt, typename OutputIt>
inline OutputIt zip(InputIt first1, InputIt last1, InputIt first2, InputIt last2, OutputIt out)
{
  //
  // Produces a container of one type constructed from two ranges
  // The constructor of the value type of the container type of the output should take two parameters, the first of which comes from the 
  // first range and the second of which comes from the second range
  //
  // parameters
  // ----------
  // first1,last1 : InputIt
  //  - the range of the first set of items to zip
  // first2,last2 : InputIt
  //  - the range of the second set of items to zip
  // out : OutputIt
  //  - output iterator where the zipped items are written
  //
  // returns
  // -------
  // OutputIt
  //  - an iterator to the element past the last element added to out
  //
  // InputIt must meet the requirements of LegacyInputIterator
  // OutputIt must meet the requirements of LegacyOutputIterator
  //

  using scalar=std::decay_t<decltype(*first1)>;
  using zipped=typename OutputIt::container_type::value_type;

  // iteratre across both ranges, zip, and write out
  while (first1!=last1&&first2!=last2)
  {
    *out++=zipped(static_cast<scalar>(*first1++),static_cast<scalar>(*first2++));
  }

  return out;
}

template <typename It, typename Compare>
inline std::size_t argmin(It first, It last, Compare comp)
{
  //
  // Find the index of the minimum element of a range
  //
  // parameters
  // ----------
  // first,last : It
  //  - the range over which to find the argmin
  // comep : OutputIt
  //  - output iterator where the complex numbers are written
  //
  // returns
  // -------
  // size_t
  //  - the index of the minimum element of a range
  //
  // It must meet the requirements of LegacyForwardterator
  //
  
  return std::min_element(first,last,comp)-first;
}

#endif // !COMMON_H__
