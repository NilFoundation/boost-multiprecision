///////////////////////////////////////////////////////////////
//  Copyright 2020 Pavel Kharitonov Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt


#ifdef _MSC_VER
#define _SCL_SECURE_NO_WARNINGS
#endif

#include "test.hpp"

#if !defined(TEST_CPP_INT)
#define TEST_MPZ
#define TEST_CPP_INT

#ifdef _MSC_VER
#pragma message("CAUTION!!: No backend type specified so testing everything.... this will take some time!!")
#endif
#ifdef __GNUC__
#pragma warning "CAUTION!!: No backend type specified so testing everything.... this will take some time!!"
#endif

#endif

#ifdef TEST_CPP_INT
#include <boost/multiprecision/cpp_int.hpp>
#endif

#include <boost/multiprecision/multiexp.hpp>

template <typename T>
void test()
{
   using namespace boost::multiprecision;

   BOOST_CHECK_EQUAL(multiexp({T(2), T(3)}, {T(5), T(3)}, T(1), T(2), T(2), T(2)), T(19));
}

int main()
{
   using namespace boost::multiprecision;

   #ifdef TEST_CPP_INT
      test<cpp_int>();
   #endif

   return boost::report_errors();
}