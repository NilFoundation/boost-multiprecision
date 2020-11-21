///////////////////////////////////////////////////////////////
//  Copyright 2020 Mikhail Komarov. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt

#ifdef _MSC_VER
#define _SCL_SECURE_NO_WARNINGS
#endif

#include "test.hpp"

#if !defined(TEST_MPZ) && !defined(TEST_TOMMATH) && !defined(TEST_CPP_INT)
#define TEST_TOMMATH
#define TEST_MPZ
#define TEST_CPP_INT

#ifdef _MSC_VER
#pragma message("CAUTION!!: No backend type specified so testing everything.... this will take some time!!")
#endif
#ifdef __GNUC__
#pragma warning "CAUTION!!: No backend type specified so testing everything.... this will take some time!!"
#endif

#endif

#if defined(TEST_MPZ)
#include <boost/multiprecision/gmp.hpp>
#endif
#if defined(TEST_TOMMATH)
#include <boost/multiprecision/tommath.hpp>
#endif
#if defined(TEST_CPP_INT)
#include <boost/multiprecision/cpp_int.hpp>
#endif

#include <boost/multiprecision/wnaf.hpp>

int main()
{
   using namespace boost::multiprecision;

#if defined(TEST_CPP_INT)
#endif
#if defined(TEST_MPZ)
#endif
#if defined(TEST_TOMMATH)
#endif

   return boost::report_errors();
}