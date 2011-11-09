///////////////////////////////////////////////////////////////
//  Copyright Christopher Kormanyos 2002 - 2011.
//  Copyright 2011 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_
//
// This work is based on an earlier work:
// "Algorithm 910: A Portable C++ Multiple-Precision System for Special-Function Calculations",
// in ACM TOMS, {VOL 37, ISSUE 4, (February 2011)} (C) ACM, 2011. http://doi.acm.org/10.1145/1916461.1916469

#include <boost/detail/lightweight_test.hpp>
#include <boost/array.hpp>
#include "test.hpp"

#if !defined(TEST_MPF_50) && !defined(TEST_MPF) && !defined(TEST_BACKEND) && !defined(TEST_MPZ) && !defined(TEST_MP_FLOAT) && !defined(TEST_MPFR) && !defined(TEST_MPFR_50) && !defined(TEST_MPQ)
#  define TEST_MPF_50
//#  define TEST_MPF
#  define TEST_BACKEND
#  define TEST_MP_FLOAT

#ifdef _MSC_VER
#pragma message("CAUTION!!: No backend type specified so testing everything.... this will take some time!!")
#endif
#ifdef __GNUC__
#pragma warning "CAUTION!!: No backend type specified so testing everything.... this will take some time!!"
#endif

#endif

#if defined(TEST_MPF_50)
#include <boost/multiprecision/gmp.hpp>
#endif
#if defined(TEST_MPFR_50)
#include <boost/multiprecision/mpfr.hpp>
#endif
#ifdef TEST_BACKEND
#include <boost/multiprecision/concepts/mp_number_architypes.hpp>
#endif
#ifdef TEST_MP_FLOAT
#include <boost/multiprecision/mp_float.hpp>
#endif

template <class T>
void test()
{
   //
   // Test with some exact binary values as input - this tests our code
   // rather than the test data:
   //
   static const boost::array<boost::array<T, 2>, 6> exact_data =
   {{
      {{ 0.5, "0.523598775598298873077107230546583814032861566562517636829157432051302734381034833104672470890352844663691347752213717775" }},
      {{ 0.25, "0.252680255142078653485657436993710972252193733096838193633923778740575060481021222411748742228014601605092602909414066566" }},
      {{0.75, "0.848062078981481008052944338998418080073366213263112642860718163570200821228474234349189801731957230300995227265307531834" }},
      {{std::ldexp(1.0, -20), "9.53674316406394560289664793089102218648031077292419572854816420395098616062014311172490017625353237219958438022056661501e-7" }},
      {{ 1 - std::ldexp(1.0, -20), "1.56941525875313420204921285316218397515809899320201864334535204504240776023375739189119474528488143494473216475057072728" }},
      {{ 1, "1.57079632679489661923132169163975144209858469968755291048747229615390820314310449931401741267105853399107404325664115332354692230477529111586267970406424055872514205135096926055277982231147447746519098221440548783296672306423782411689339158263560095457282428346173017430522716332410669680363012457064" }},
   }};
   unsigned max_err = 0;
   for(unsigned k = 0; k < exact_data.size(); k++)
   {
      T val = asin(exact_data[k][0]);
      T e = relative_error(val, exact_data[k][1]);
      unsigned err = e.template convert_to<unsigned>();
      if(err > max_err)
         max_err = err;
      val = asin(-exact_data[k][0]);
      e = relative_error(val, T(-exact_data[k][1]));
      err = e.template convert_to<unsigned>();
      if(err > max_err)
      {
         std::cout << val << std::endl;
         std::cout << -exact_data[k][0] << std::endl;
         max_err = err;
      }
   }
   std::cout << "Max error was: " << max_err << std::endl;
   BOOST_TEST(max_err < 20);
   BOOST_TEST(asin(T(0)) == 0);
}


int main()
{
#ifdef TEST_BACKEND
   test<boost::multiprecision::mp_number<boost::multiprecision::concepts::mp_number_backend_float_architype> >();
#endif
#ifdef TEST_MPF_50
   test<boost::multiprecision::mpf_float_50>();
   test<boost::multiprecision::mpf_float_100>();
#endif
#ifdef TEST_MPFR_50
   test<boost::multiprecision::mpfr_float_50>();
   test<boost::multiprecision::mpfr_float_100>();
#endif
#ifdef TEST_MP_FLOAT
   test<boost::multiprecision::mp_float_50>();
   test<boost::multiprecision::mp_float_100>();
#endif
   return boost::report_errors();
}


