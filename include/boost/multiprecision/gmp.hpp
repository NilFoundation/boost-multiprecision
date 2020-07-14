///////////////////////////////////////////////////////////////////////////////
//  Copyright 2011 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_ER_GMP_BACKEND_HPP
#define BOOST_MATH_ER_GMP_BACKEND_HPP

#include <boost/multiprecision/number.hpp>
#include <boost/multiprecision/debug_adaptor.hpp>
#include <boost/multiprecision/detail/integer_ops.hpp>
#include <boost/multiprecision/detail/big_lanczos.hpp>
#include <boost/multiprecision/detail/digits.hpp>
#include <boost/multiprecision/detail/atomic.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/cstdint.hpp>
#include <boost/functional/hash_fwd.hpp>
//
// Some includes we need from Boost.Math, since we rely on that library to provide these functions:
//
#include <boost/math/special_functions/asinh.hpp>
#include <boost/math/special_functions/acosh.hpp>
#include <boost/math/special_functions/atanh.hpp>
#include <boost/math/special_functions/cbrt.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/gamma.hpp>

#ifdef BOOST_MSVC
#pragma warning(push)
#pragma warning(disable : 4127)
#endif
#include <gmp.h>
#ifdef BOOST_MSVC
#pragma warning(pop)
#endif

#if defined(__MPIR_VERSION) && defined(__MPIR_VERSION_MINOR) && defined(__MPIR_VERSION_PATCHLEVEL)
#define BOOST_MP_MPIR_VERSION (__MPIR_VERSION * 10000 + __MPIR_VERSION_MINOR * 100 + __MPIR_VERSION_PATCHLEVEL)
#else
#define BOOST_MP_MPIR_VERSION 0
#endif

#include <cctype>
#include <cmath>
#include <limits>
#include <climits>

namespace boost {
namespace multiprecision {
namespace backends {

#ifdef BOOST_MSVC
// warning C4127: conditional expression is constant
#pragma warning(push)
#pragma warning(disable : 4127)
#endif

template <unsigned digits10>
struct gmp_float;
struct gmp_int;
struct gmp_rational;

} // namespace backends

template <>
struct number_category<backends::gmp_int> : public mpl::int_<number_kind_integer>
{};
template <>
struct number_category<backends::gmp_rational> : public mpl::int_<number_kind_rational>
{};
template <unsigned digits10>
struct number_category<backends::gmp_float<digits10> > : public mpl::int_<number_kind_floating_point>
{};

namespace backends {
//
// Within this file, the only functions we mark as noexcept are those that manipulate
// (but don't create) an mpf_t.  All other types may allocate at pretty much any time
// via a user-supplied allocator, and therefore throw.
//
namespace detail {

template <unsigned digits10>
struct gmp_float_imp
{
#ifdef BOOST_HAS_LONG_LONG
   typedef mpl::list<long, boost::long_long_type>           signed_types;
   typedef mpl::list<unsigned long, boost::ulong_long_type> unsigned_types;
#else
   typedef mpl::list<long>          signed_types;
   typedef mpl::list<unsigned long> unsigned_types;
#endif
   typedef mpl::list<double, long double> float_types;
   typedef long                           exponent_type;

   gmp_float_imp() BOOST_NOEXCEPT
   {
      m_data[0]._mp_d = 0; // uninitialized m_data
   }

   gmp_float_imp(const gmp_float_imp& o)
   {
      //
      // We have to do an init followed by a set here, otherwise *this may be at
      // a lower precision than o: seems like mpf_init_set copies just enough bits
      // to get the right value, but if it's then used in further calculations
      // things go badly wrong!!
      //
      mpf_init2(m_data, mpf_get_prec(o.data()));
      if (o.m_data[0]._mp_d)
         mpf_set(m_data, o.m_data);
   }
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
   gmp_float_imp(gmp_float_imp&& o) BOOST_NOEXCEPT
   {
      m_data[0]         = o.m_data[0];
      o.m_data[0]._mp_d = 0;
   }
#endif
   gmp_float_imp& operator=(const gmp_float_imp& o)
   {
      if (m_data[0]._mp_d == 0)
         mpf_init2(m_data, mpf_get_prec(o.data()));
      if (mpf_get_prec(data()) != mpf_get_prec(o.data()))
      {
         mpf_t t;
         mpf_init2(t, mpf_get_prec(o.data()));
         mpf_set(t, o.data());
         mpf_swap(data(), t);
         mpf_clear(t);
      }
      else
      {
         if (o.m_data[0]._mp_d)
            mpf_set(m_data, o.m_data);
      }
      return *this;
   }
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
   gmp_float_imp& operator=(gmp_float_imp&& o) BOOST_NOEXCEPT
   {
      mpf_swap(m_data, o.m_data);
      return *this;
   }
#endif

#ifdef BOOST_HAS_LONG_LONG
#if defined(ULLONG_MAX) && (ULLONG_MAX == ULONG_MAX)
   gmp_float_imp& operator=(boost::ulong_long_type i)
   {
      *this = static_cast<unsigned long>(i);
      return *this;
   }
#else
   gmp_float_imp& operator=(boost::ulong_long_type i)
   {
      if (m_data[0]._mp_d == 0)
         mpf_init2(m_data, multiprecision::detail::digits10_2_2(digits10 ? digits10 : (unsigned)get_default_precision()));
      boost::ulong_long_type mask  = ((((1uLL << (std::numeric_limits<unsigned long>::digits - 1)) - 1) << 1) | 1uLL);
      unsigned               shift = 0;
      mpf_t                  t;
      mpf_init2(t, multiprecision::detail::digits10_2_2(digits10 ? digits10 : (unsigned)get_default_precision()));
      mpf_set_ui(m_data, 0);
      while (i)
      {
         mpf_set_ui(t, static_cast<unsigned long>(i & mask));
         if (shift)
            mpf_mul_2exp(t, t, shift);
         mpf_add(m_data, m_data, t);
         shift += std::numeric_limits<unsigned long>::digits;
         i >>= std::numeric_limits<unsigned long>::digits;
      }
      mpf_clear(t);
      return *this;
   }
#endif
   gmp_float_imp& operator=(boost::long_long_type i)
   {
      if (m_data[0]._mp_d == 0)
         mpf_init2(m_data, multiprecision::detail::digits10_2_2(digits10 ? digits10 : (unsigned)get_default_precision()));
      bool neg = i < 0;
      *this    = static_cast<boost::ulong_long_type>(boost::multiprecision::detail::unsigned_abs(i));
      if (neg)
         mpf_neg(m_data, m_data);
      return *this;
   }
#endif
   gmp_float_imp& operator=(unsigned long i)
   {
      if (m_data[0]._mp_d == 0)
         mpf_init2(m_data, multiprecision::detail::digits10_2_2(digits10 ? digits10 : (unsigned)get_default_precision()));
      mpf_set_ui(m_data, i);
      return *this;
   }
   gmp_float_imp& operator=(long i)
   {
      if (m_data[0]._mp_d == 0)
         mpf_init2(m_data, multiprecision::detail::digits10_2_2(digits10 ? digits10 : (unsigned)get_default_precision()));
      mpf_set_si(m_data, i);
      return *this;
   }
   gmp_float_imp& operator=(double d)
   {
      if (m_data[0]._mp_d == 0)
         mpf_init2(m_data, multiprecision::detail::digits10_2_2(digits10 ? digits10 : (unsigned)get_default_precision()));
      mpf_set_d(m_data, d);
      return *this;
   }
   gmp_float_imp& operator=(long double a)
   {
      using std::floor;
      using std::frexp;
      using std::ldexp;

      if (m_data[0]._mp_d == 0)
         mpf_init2(m_data, multiprecision::detail::digits10_2_2(digits10 ? digits10 : (unsigned)get_default_precision()));

      if (a == 0)
      {
         mpf_set_si(m_data, 0);
         return *this;
      }

      if (a == 1)
      {
         mpf_set_si(m_data, 1);
         return *this;
      }

      BOOST_ASSERT(!(boost::math::isinf)(a));
      BOOST_ASSERT(!(boost::math::isnan)(a));

      int         e;
      long double f, term;
      mpf_set_ui(m_data, 0u);

      f = frexp(a, &e);

      static const int shift = std::numeric_limits<int>::digits - 1;

      while (f)
      {
         // extract int sized bits from f:
         f    = ldexp(f, shift);
         term = floor(f);
         e -= shift;
         mpf_mul_2exp(m_data, m_data, shift);
         if (term > 0)
            mpf_add_ui(m_data, m_data, static_cast<unsigned>(term));
         else
            mpf_sub_ui(m_data, m_data, static_cast<unsigned>(-term));
         f -= term;
      }
      if (e > 0)
         mpf_mul_2exp(m_data, m_data, e);
      else if (e < 0)
         mpf_div_2exp(m_data, m_data, -e);
      return *this;
   }
   gmp_float_imp& operator=(const char* s)
   {
      if (m_data[0]._mp_d == 0)
         mpf_init2(m_data, multiprecision::detail::digits10_2_2(digits10 ? digits10 : (unsigned)get_default_precision()));
      if (s && (*s == '+'))
         ++s;  // Leading "+" sign not supported by mpf_set_str:
      if (0 != mpf_set_str(m_data, s, 10))
         BOOST_THROW_EXCEPTION(std::runtime_error(std::string("The string \"") + s + std::string("\"could not be interpreted as a valid floating point number.")));
      return *this;
   }
   void swap(gmp_float_imp& o) BOOST_NOEXCEPT
   {
      mpf_swap(m_data, o.m_data);
   }
   std::string str(std::streamsize digits, std::ios_base::fmtflags f) const
   {
      BOOST_ASSERT(m_data[0]._mp_d);

      bool            scientific = (f & std::ios_base::scientific) == std::ios_base::scientific;
      bool            fixed      = (f & std::ios_base::fixed) == std::ios_base::fixed;
      std::streamsize org_digits(digits);

      if (scientific && digits)
         ++digits;

      std::string result;
      mp_exp_t    e;
      void* (*alloc_func_ptr)(size_t);
      void* (*realloc_func_ptr)(void*, size_t, size_t);
      void (*free_func_ptr)(void*, size_t);
      mp_get_memory_functions(&alloc_func_ptr, &realloc_func_ptr, &free_func_ptr);

      if (mpf_sgn(m_data) == 0)
      {
         e      = 0;
         result = "0";
         if (fixed && digits)
            ++digits;
      }
      else
      {
         char* ps = mpf_get_str(0, &e, 10, static_cast<std::size_t>(digits), m_data);
         --e; // To match with what our formatter expects.
         if (fixed && e != -1)
         {
            // Oops we actually need a different number of digits to what we asked for:
            (*free_func_ptr)((void*)ps, std::strlen(ps) + 1);
            digits += e + 1;
            if (digits == 0)
            {
               // We need to get *all* the digits and then possibly round up,
               // we end up with either "0" or "1" as the result.
               ps = mpf_get_str(0, &e, 10, 0, m_data);
               --e;
               unsigned offset = *ps == '-' ? 1 : 0;
               if (ps[offset] > '5')
               {
                  ++e;
                  ps[offset]     = '1';
                  ps[offset + 1] = 0;
               }
               else if (ps[offset] == '5')
               {
                  unsigned i        = offset + 1;
                  bool     round_up = false;
                  while (ps[i] != 0)
                  {
                     if (ps[i] != '0')
                     {
                        round_up = true;
                        break;
                     }
                     ++i;
                  }
                  if (round_up)
                  {
                     ++e;
                     ps[offset]     = '1';
                     ps[offset + 1] = 0;
                  }
                  else
                  {
                     ps[offset]     = '0';
                     ps[offset + 1] = 0;
                  }
               }
               else
               {
                  ps[offset]     = '0';
                  ps[offset + 1] = 0;
               }
            }
            else if (digits > 0)
            {
               mp_exp_t old_e = e;
               ps             = mpf_get_str(0, &e, 10, static_cast<std::size_t>(digits), m_data);
               --e; // To match with what our formatter expects.
               if (old_e > e)
               {
                  // in some cases, when we ask for more digits of precision, it will
                  // change the number of digits to the left of the decimal, if that
                  // happens, account for it here.
                  // example: cout << fixed << setprecision(3) << mpf_float_50("99.9809")
                  digits -= old_e - e;
                  ps = mpf_get_str(0, &e, 10, static_cast<std::size_t>(digits), m_data);
                  --e; // To match with what our formatter expects.
               }
            }
            else
            {
               ps = mpf_get_str(0, &e, 10, 1, m_data);
               --e;
               unsigned offset = *ps == '-' ? 1 : 0;
               ps[offset]      = '0';
               ps[offset + 1]  = 0;
            }
         }
         result = ps;
         (*free_func_ptr)((void*)ps, std::strlen(ps) + 1);
      }
      boost::multiprecision::detail::format_float_string(result, e, org_digits, f, mpf_sgn(m_data) == 0);
      return result;
   }
   ~gmp_float_imp() BOOST_NOEXCEPT
   {
      if (m_data[0]._mp_d)
         mpf_clear(m_data);
   }
   void negate() BOOST_NOEXCEPT
   {
      BOOST_ASSERT(m_data[0]._mp_d);
      mpf_neg(m_data, m_data);
   }
   int compare(const gmp_float<digits10>& o) const BOOST_NOEXCEPT
   {
      BOOST_ASSERT(m_data[0]._mp_d && o.m_data[0]._mp_d);
      return mpf_cmp(m_data, o.m_data);
   }
   int compare(long i) const BOOST_NOEXCEPT
   {
      BOOST_ASSERT(m_data[0]._mp_d);
      return mpf_cmp_si(m_data, i);
   }
   int compare(unsigned long i) const BOOST_NOEXCEPT
   {
      BOOST_ASSERT(m_data[0]._mp_d);
      return mpf_cmp_ui(m_data, i);
   }
   template <class V>
   typename enable_if<is_arithmetic<V>, int>::type compare(V v) const
   {
      gmp_float<digits10> d;
      d = v;
      return compare(d);
   }
   mpf_t& data() BOOST_NOEXCEPT
   {
      BOOST_ASSERT(m_data[0]._mp_d);
      return m_data;
   }
   const mpf_t& data() const BOOST_NOEXCEPT
   {
      BOOST_ASSERT(m_data[0]._mp_d);
      return m_data;
   }

 protected:
   mpf_t            m_data;
   static boost::multiprecision::detail::precision_type& get_default_precision() BOOST_NOEXCEPT
   {
      static boost::multiprecision::detail::precision_type val(50);
      return val;
   }
};

} // namespace detail

struct gmp_int;
struct gmp_rational;

template <unsigned digits10>
struct gmp_float : public detail::gmp_float_imp<digits10>
{
   gmp_float()
   {
      mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(digits10));
   }
   gmp_float(const gmp_float& o) : detail::gmp_float_imp<digits10>(o) {}
   template <unsigned D>
   gmp_float(const gmp_float<D>& o, typename enable_if_c<D <= digits10>::type* = 0);
   template <unsigned D>
   explicit gmp_float(const gmp_float<D>& o, typename disable_if_c<D <= digits10>::type* = 0);
   gmp_float(const gmp_int& o);
   gmp_float(const gmp_rational& o);
   gmp_float(const mpf_t val)
   {
      mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(digits10));
      mpf_set(this->m_data, val);
   }
   gmp_float(const mpz_t val)
   {
      mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(digits10));
      mpf_set_z(this->m_data, val);
   }
   gmp_float(const mpq_t val)
   {
      mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(digits10));
      mpf_set_q(this->m_data, val);
   }
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
   gmp_float(gmp_float&& o) BOOST_NOEXCEPT : detail::gmp_float_imp<digits10>(static_cast<detail::gmp_float_imp<digits10>&&>(o))
   {}
#endif
   gmp_float& operator=(const gmp_float& o)
   {
      *static_cast<detail::gmp_float_imp<digits10>*>(this) = static_cast<detail::gmp_float_imp<digits10> const&>(o);
      return *this;
   }
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
   gmp_float& operator=(gmp_float&& o) BOOST_NOEXCEPT
   {
      *static_cast<detail::gmp_float_imp<digits10>*>(this) = static_cast<detail::gmp_float_imp<digits10>&&>(o);
      return *this;
   }
#endif
   template <unsigned D>
   gmp_float& operator=(const gmp_float<D>& o);
   gmp_float& operator=(const gmp_int& o);
   gmp_float& operator=(const gmp_rational& o);
   gmp_float& operator=(const mpf_t val)
   {
      if (this->m_data[0]._mp_d == 0)
         mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(digits10));
      mpf_set(this->m_data, val);
      return *this;
   }
   gmp_float& operator=(const mpz_t val)
   {
      if (this->m_data[0]._mp_d == 0)
         mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(digits10));
      mpf_set_z(this->m_data, val);
      return *this;
   }
   gmp_float& operator=(const mpq_t val)
   {
      if (this->m_data[0]._mp_d == 0)
         mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(digits10));
      mpf_set_q(this->m_data, val);
      return *this;
   }
   template <class V>
   gmp_float& operator=(const V& v)
   {
      *static_cast<detail::gmp_float_imp<digits10>*>(this) = v;
      return *this;
   }
};

template <>
struct gmp_float<0> : public detail::gmp_float_imp<0>
{
   //
   // We have a problem with mpf_t in that the precision we request isn't what we get.
   // As a result the front end can end up chasing it's tail trying to create a variable
   // with the the correct precision to hold the result of an expression.
   // See: https://github.com/boostorg/multiprecision/issues/164
   // The problem is made worse by the fact that our conversions from base10 to 2 and
   // vice-versa do not exactly round trip (and probably never will).
   // The workaround is to keep track of the precision requested, and always return
   // that as the current actual precision.
   //
 private:
   unsigned requested_precision;

 public:
   gmp_float() : requested_precision(get_default_precision())
   {
      mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(requested_precision));
   }
   gmp_float(const mpf_t val) : requested_precision(get_default_precision())
   {
      mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(requested_precision));
      mpf_set(this->m_data, val);
   }
   gmp_float(const mpz_t val) : requested_precision(get_default_precision())
   {
      mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(requested_precision));
      mpf_set_z(this->m_data, val);
   }
   gmp_float(const mpq_t val) : requested_precision(get_default_precision())
   {
      mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(requested_precision));
      mpf_set_q(this->m_data, val);
   }
   gmp_float(const gmp_float& o) : detail::gmp_float_imp<0>(o), requested_precision(o.requested_precision) {}
   template <unsigned D>
   gmp_float(const gmp_float<D>& o)
   {
      mpf_init2(this->m_data, mpf_get_prec(o.data()));
      mpf_set(this->m_data, o.data());
      requested_precision = D;
   }
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
   gmp_float(gmp_float&& o) BOOST_NOEXCEPT : detail::gmp_float_imp<0>(static_cast<detail::gmp_float_imp<0>&&>(o)), requested_precision(o.requested_precision)
   {}
#endif
   gmp_float(const gmp_int& o);
   gmp_float(const gmp_rational& o);
   gmp_float(const gmp_float& o, unsigned digits10) : requested_precision(digits10)
   {
      mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(digits10));
      mpf_set(this->m_data, o.data());
   }
   template <class V>
   gmp_float(const V& o, unsigned digits10) : requested_precision(digits10)
   {
      mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(digits10));
      *this = o;
   }

#ifndef BOOST_NO_CXX17_HDR_STRING_VIEW
   //
   // Support for new types in C++17
   //
   template <class Traits>
   gmp_float(const std::basic_string_view<char, Traits>& o, unsigned digits10) : requested_precision(digits10)
   {
      using default_ops::assign_from_string_view;
      mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(digits10));
      assign_from_string_view(*this, o);
   }
#endif
   gmp_float& operator=(const gmp_float& o)
   {
      *static_cast<detail::gmp_float_imp<0>*>(this) = static_cast<detail::gmp_float_imp<0> const&>(o);
      requested_precision                           = o.requested_precision;
      return *this;
   }
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
   gmp_float& operator=(gmp_float&& o) BOOST_NOEXCEPT
   {
      *static_cast<detail::gmp_float_imp<0>*>(this) = static_cast<detail::gmp_float_imp<0>&&>(o);
      requested_precision                           = o.requested_precision;
      return *this;
   }
#endif
   template <unsigned D>
   gmp_float& operator=(const gmp_float<D>& o)
   {
      if (this->m_data[0]._mp_d == 0)
      {
         mpf_init2(this->m_data, mpf_get_prec(o.data()));
      }
      else
      {
         mpf_set_prec(this->m_data, mpf_get_prec(o.data()));
      }
      mpf_set(this->m_data, o.data());
      requested_precision = D;
      return *this;
   }
   gmp_float& operator=(const gmp_int& o);
   gmp_float& operator=(const gmp_rational& o);
   gmp_float& operator=(const mpf_t val)
   {
      if (this->m_data[0]._mp_d == 0)
      {
         requested_precision = get_default_precision();
         mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(requested_precision));
      }
      mpf_set(this->m_data, val);
      return *this;
   }
   gmp_float& operator=(const mpz_t val)
   {
      if (this->m_data[0]._mp_d == 0)
      {
         requested_precision = get_default_precision();
         mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(requested_precision));
      }
      mpf_set_z(this->m_data, val);
      return *this;
   }
   gmp_float& operator=(const mpq_t val)
   {
      if (this->m_data[0]._mp_d == 0)
      {
         requested_precision = get_default_precision();
         mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(requested_precision));
      }
      mpf_set_q(this->m_data, val);
      return *this;
   }
   template <class V>
   gmp_float& operator=(const V& v)
   {
      *static_cast<detail::gmp_float_imp<0>*>(this) = v;
      return *this;
   }
   static unsigned default_precision() BOOST_NOEXCEPT
   {
      return get_default_precision();
   }
   static void default_precision(unsigned v) BOOST_NOEXCEPT
   {
      get_default_precision() = v;
   }
   unsigned precision() const BOOST_NOEXCEPT
   {
      return requested_precision;
   }
   void precision(unsigned digits10) BOOST_NOEXCEPT
   {
      requested_precision = digits10;
      mpf_set_prec(this->m_data, multiprecision::detail::digits10_2_2(requested_precision));
   }
   void swap(gmp_float& o)
   {
      std::swap(requested_precision, o.requested_precision);
      gmp_float_imp<0>::swap(o);
   }
};

template <unsigned digits10, class T>
inline typename enable_if_c<is_arithmetic<T>::value, bool>::type eval_eq(const gmp_float<digits10>& a, const T& b) BOOST_NOEXCEPT
{
   return a.compare(b) == 0;
}
template <unsigned digits10, class T>
inline typename enable_if_c<is_arithmetic<T>::value, bool>::type eval_lt(const gmp_float<digits10>& a, const T& b) BOOST_NOEXCEPT
{
   return a.compare(b) < 0;
}
template <unsigned digits10, class T>
inline typename enable_if_c<is_arithmetic<T>::value, bool>::type eval_gt(const gmp_float<digits10>& a, const T& b) BOOST_NOEXCEPT
{
   return a.compare(b) > 0;
}

template <unsigned D1, unsigned D2>
inline void eval_add(gmp_float<D1>& result, const gmp_float<D2>& o)
{
   mpf_add(result.data(), result.data(), o.data());
}
template <unsigned D1, unsigned D2>
inline void eval_subtract(gmp_float<D1>& result, const gmp_float<D2>& o)
{
   mpf_sub(result.data(), result.data(), o.data());
}
template <unsigned D1, unsigned D2>
inline void eval_multiply(gmp_float<D1>& result, const gmp_float<D2>& o)
{
   mpf_mul(result.data(), result.data(), o.data());
}
template <unsigned digits10>
inline bool eval_is_zero(const gmp_float<digits10>& val) BOOST_NOEXCEPT
{
   return mpf_sgn(val.data()) == 0;
}
template <unsigned D1, unsigned D2>
inline void eval_divide(gmp_float<D1>& result, const gmp_float<D2>& o)
{
   if (eval_is_zero(o))
      BOOST_THROW_EXCEPTION(std::overflow_error("Division by zero."));
   mpf_div(result.data(), result.data(), o.data());
}
template <unsigned digits10>
inline void eval_add(gmp_float<digits10>& result, unsigned long i)
{
   mpf_add_ui(result.data(), result.data(), i);
}
template <unsigned digits10>
inline void eval_subtract(gmp_float<digits10>& result, unsigned long i)
{
   mpf_sub_ui(result.data(), result.data(), i);
}
template <unsigned digits10>
inline void eval_multiply(gmp_float<digits10>& result, unsigned long i)
{
   mpf_mul_ui(result.data(), result.data(), i);
}
template <unsigned digits10>
inline void eval_divide(gmp_float<digits10>& result, unsigned long i)
{
   if (i == 0)
      BOOST_THROW_EXCEPTION(std::overflow_error("Division by zero."));
   mpf_div_ui(result.data(), result.data(), i);
}
template <unsigned digits10>
inline void eval_add(gmp_float<digits10>& result, long i)
{
   if (i > 0)
      mpf_add_ui(result.data(), result.data(), i);
   else
      mpf_sub_ui(result.data(), result.data(), boost::multiprecision::detail::unsigned_abs(i));
}
template <unsigned digits10>
inline void eval_subtract(gmp_float<digits10>& result, long i)
{
   if (i > 0)
      mpf_sub_ui(result.data(), result.data(), i);
   else
      mpf_add_ui(result.data(), result.data(), boost::multiprecision::detail::unsigned_abs(i));
}
template <unsigned digits10>
inline void eval_multiply(gmp_float<digits10>& result, long i)
{
   mpf_mul_ui(result.data(), result.data(), boost::multiprecision::detail::unsigned_abs(i));
   if (i < 0)
      mpf_neg(result.data(), result.data());
}
template <unsigned digits10>
inline void eval_divide(gmp_float<digits10>& result, long i)
{
   if (i == 0)
      BOOST_THROW_EXCEPTION(std::overflow_error("Division by zero."));
   mpf_div_ui(result.data(), result.data(), boost::multiprecision::detail::unsigned_abs(i));
   if (i < 0)
      mpf_neg(result.data(), result.data());
}
//
// Specialised 3 arg versions of the basic operators:
//
template <unsigned D1, unsigned D2, unsigned D3>
inline void eval_add(gmp_float<D1>& a, const gmp_float<D2>& x, const gmp_float<D3>& y)
{
   mpf_add(a.data(), x.data(), y.data());
}
template <unsigned D1, unsigned D2>
inline void eval_add(gmp_float<D1>& a, const gmp_float<D2>& x, unsigned long y)
{
   mpf_add_ui(a.data(), x.data(), y);
}
template <unsigned D1, unsigned D2>
inline void eval_add(gmp_float<D1>& a, const gmp_float<D2>& x, long y)
{
   if (y < 0)
      mpf_sub_ui(a.data(), x.data(), boost::multiprecision::detail::unsigned_abs(y));
   else
      mpf_add_ui(a.data(), x.data(), y);
}
template <unsigned D1, unsigned D2>
inline void eval_add(gmp_float<D1>& a, unsigned long x, const gmp_float<D2>& y)
{
   mpf_add_ui(a.data(), y.data(), x);
}
template <unsigned D1, unsigned D2>
inline void eval_add(gmp_float<D1>& a, long x, const gmp_float<D2>& y)
{
   if (x < 0)
   {
      mpf_ui_sub(a.data(), boost::multiprecision::detail::unsigned_abs(x), y.data());
      mpf_neg(a.data(), a.data());
   }
   else
      mpf_add_ui(a.data(), y.data(), x);
}
template <unsigned D1, unsigned D2, unsigned D3>
inline void eval_subtract(gmp_float<D1>& a, const gmp_float<D2>& x, const gmp_float<D3>& y)
{
   mpf_sub(a.data(), x.data(), y.data());
}
template <unsigned D1, unsigned D2>
inline void eval_subtract(gmp_float<D1>& a, const gmp_float<D2>& x, unsigned long y)
{
   mpf_sub_ui(a.data(), x.data(), y);
}
template <unsigned D1, unsigned D2>
inline void eval_subtract(gmp_float<D1>& a, const gmp_float<D2>& x, long y)
{
   if (y < 0)
      mpf_add_ui(a.data(), x.data(), boost::multiprecision::detail::unsigned_abs(y));
   else
      mpf_sub_ui(a.data(), x.data(), y);
}
template <unsigned D1, unsigned D2>
inline void eval_subtract(gmp_float<D1>& a, unsigned long x, const gmp_float<D2>& y)
{
   mpf_ui_sub(a.data(), x, y.data());
}
template <unsigned D1, unsigned D2>
inline void eval_subtract(gmp_float<D1>& a, long x, const gmp_float<D2>& y)
{
   if (x < 0)
   {
      mpf_add_ui(a.data(), y.data(), boost::multiprecision::detail::unsigned_abs(x));
      mpf_neg(a.data(), a.data());
   }
   else
      mpf_ui_sub(a.data(), x, y.data());
}

template <unsigned D1, unsigned D2, unsigned D3>
inline void eval_multiply(gmp_float<D1>& a, const gmp_float<D2>& x, const gmp_float<D3>& y)
{
   mpf_mul(a.data(), x.data(), y.data());
}
template <unsigned D1, unsigned D2>
inline void eval_multiply(gmp_float<D1>& a, const gmp_float<D2>& x, unsigned long y)
{
   mpf_mul_ui(a.data(), x.data(), y);
}
template <unsigned D1, unsigned D2>
inline void eval_multiply(gmp_float<D1>& a, const gmp_float<D2>& x, long y)
{
   if (y < 0)
   {
      mpf_mul_ui(a.data(), x.data(), boost::multiprecision::detail::unsigned_abs(y));
      a.negate();
   }
   else
      mpf_mul_ui(a.data(), x.data(), y);
}
template <unsigned D1, unsigned D2>
inline void eval_multiply(gmp_float<D1>& a, unsigned long x, const gmp_float<D2>& y)
{
   mpf_mul_ui(a.data(), y.data(), x);
}
template <unsigned D1, unsigned D2>
inline void eval_multiply(gmp_float<D1>& a, long x, const gmp_float<D2>& y)
{
   if (x < 0)
   {
      mpf_mul_ui(a.data(), y.data(), boost::multiprecision::detail::unsigned_abs(x));
      mpf_neg(a.data(), a.data());
   }
   else
      mpf_mul_ui(a.data(), y.data(), x);
}

template <unsigned D1, unsigned D2, unsigned D3>
inline void eval_divide(gmp_float<D1>& a, const gmp_float<D2>& x, const gmp_float<D3>& y)
{
   if (eval_is_zero(y))
      BOOST_THROW_EXCEPTION(std::overflow_error("Division by zero."));
   mpf_div(a.data(), x.data(), y.data());
}
template <unsigned D1, unsigned D2>
inline void eval_divide(gmp_float<D1>& a, const gmp_float<D2>& x, unsigned long y)
{
   if (y == 0)
      BOOST_THROW_EXCEPTION(std::overflow_error("Division by zero."));
   mpf_div_ui(a.data(), x.data(), y);
}
template <unsigned D1, unsigned D2>
inline void eval_divide(gmp_float<D1>& a, const gmp_float<D2>& x, long y)
{
   if (y == 0)
      BOOST_THROW_EXCEPTION(std::overflow_error("Division by zero."));
   if (y < 0)
   {
      mpf_div_ui(a.data(), x.data(), boost::multiprecision::detail::unsigned_abs(y));
      a.negate();
   }
   else
      mpf_div_ui(a.data(), x.data(), y);
}
template <unsigned D1, unsigned D2>
inline void eval_divide(gmp_float<D1>& a, unsigned long x, const gmp_float<D2>& y)
{
   if (eval_is_zero(y))
      BOOST_THROW_EXCEPTION(std::overflow_error("Division by zero."));
   mpf_ui_div(a.data(), x, y.data());
}
template <unsigned D1, unsigned D2>
inline void eval_divide(gmp_float<D1>& a, long x, const gmp_float<D2>& y)
{
   if (eval_is_zero(y))
      BOOST_THROW_EXCEPTION(std::overflow_error("Division by zero."));
   if (x < 0)
   {
      mpf_ui_div(a.data(), boost::multiprecision::detail::unsigned_abs(x), y.data());
      mpf_neg(a.data(), a.data());
   }
   else
      mpf_ui_div(a.data(), x, y.data());
}

template <unsigned digits10>
inline int eval_get_sign(const gmp_float<digits10>& val) BOOST_NOEXCEPT
{
   return mpf_sgn(val.data());
}

template <unsigned digits10>
inline void eval_convert_to(unsigned long* result, const gmp_float<digits10>& val) BOOST_NOEXCEPT
{
   if (0 == mpf_fits_ulong_p(val.data()))
      *result = (std::numeric_limits<unsigned long>::max)();
   else
      *result = (unsigned long)mpf_get_ui(val.data());
}
template <unsigned digits10>
inline void eval_convert_to(long* result, const gmp_float<digits10>& val) BOOST_NOEXCEPT
{
   if (0 == mpf_fits_slong_p(val.data()))
   {
      *result = (std::numeric_limits<long>::max)();
      *result *= mpf_sgn(val.data());
   }
   else
      *result = (long)mpf_get_si(val.data());
}
template <unsigned digits10>
inline void eval_convert_to(double* result, const gmp_float<digits10>& val) BOOST_NOEXCEPT
{
   *result = mpf_get_d(val.data());
}
#ifdef BOOST_HAS_LONG_LONG
template <unsigned digits10>
inline void eval_convert_to(boost::long_long_type* result, const gmp_float<digits10>& val)
{
   gmp_float<digits10> t(val);
   if (eval_get_sign(t) < 0)
      t.negate();

   long digits = std::numeric_limits<boost::long_long_type>::digits - std::numeric_limits<long>::digits;

   if (digits > 0)
      mpf_div_2exp(t.data(), t.data(), digits);

   if (!mpf_fits_slong_p(t.data()))
   {
      if (eval_get_sign(val) < 0)
         *result = (std::numeric_limits<boost::long_long_type>::min)();
      else
         *result = (std::numeric_limits<boost::long_long_type>::max)();
      return;
   };

   *result = mpf_get_si(t.data());
   while (digits > 0)
   {
      *result <<= digits;
      digits -= std::numeric_limits<unsigned long>::digits;
      mpf_mul_2exp(t.data(), t.data(), digits >= 0 ? std::numeric_limits<unsigned long>::digits : std::numeric_limits<unsigned long>::digits + digits);
      unsigned long l = (unsigned long)mpf_get_ui(t.data());
      if (digits < 0)
         l >>= -digits;
      *result |= l;
   }
   if (eval_get_sign(val) < 0)
      *result = -*result;
}
template <unsigned digits10>
inline void eval_convert_to(boost::ulong_long_type* result, const gmp_float<digits10>& val)
{
   gmp_float<digits10> t(val);

   long digits = std::numeric_limits<boost::long_long_type>::digits - std::numeric_limits<long>::digits;

   if (digits > 0)
      mpf_div_2exp(t.data(), t.data(), digits);

   if (!mpf_fits_ulong_p(t.data()))
   {
      *result = (std::numeric_limits<boost::long_long_type>::max)();
      return;
   }

   *result = mpf_get_ui(t.data());
   while (digits > 0)
   {
      *result <<= digits;
      digits -= std::numeric_limits<unsigned long>::digits;
      mpf_mul_2exp(t.data(), t.data(), digits >= 0 ? std::numeric_limits<unsigned long>::digits : std::numeric_limits<unsigned long>::digits + digits);
      unsigned long l = (unsigned long)mpf_get_ui(t.data());
      if (digits < 0)
         l >>= -digits;
      *result |= l;
   }
}
#endif

std::vector<long> mpf_wnaf(const size_t window_size, const number<gmp_int> &scalar) {
   const size_t length = scalar.max_bits();    // upper bound
   std::vector<long> res(length + 1);
   number_type<n> c = scalar;
   long j = 0;
   while (!c.is_zero()) {
      long u;
      if ((c.data[0] & 1) == 1) {
         u = c.data[0] % (1u << (window_size + 1));
         if (u > (1 << window_size)) {
            u = u - (1 << (window_size + 1));
         }

         if (u > 0) {
            mpn_sub_1(c.data, c.data, n, u);
         } else {
            mpn_add_1(c.data, c.data, n, -u);
         }
      } else {
         u = 0;
      }
         res[j] = u;
         ++j;

         mpn_rshift(c.data, c.data, n, 1);    // c = c/2
      }

   return res;
}

        /**
         * In additive notation, use wNAF exponentiation (with the given window size) to compute scalar * base.
         */
template<typename T>
   T fixed_window_wnaf_exp(const size_t window_size, const T &base, const number<gmp_int> &scalar) {
      std::vector<long> naf = mpf_wnaf(window_size, scalar);
      std::vector<T> table(1ul << (window_size - 1));
      T tmp = base;
      T dbl = base.dbl();
      for (size_t i = 0; i < 1ul << (window_size - 1); ++i) {
         table[i] = tmp;
         tmp = tmp + dbl;
      }

      T res = T::zero();
      bool found_nonzero = false;
      for (long i = naf.size() - 1; i >= 0; --i) {
         if (found_nonzero) {
            res = res.dbl();
         }

         if (naf[i] != 0) {
            found_nonzero = true;
            if (naf[i] > 0) {
               res = res + table[naf[i] / 2];
            } else {
               res = res - table[(-naf[i]) / 2];
            }
         }
      }

   return res;
}

        /**
         * In additive notation, use wNAF exponentiation (with the window size determined by T) to compute scalar *
         * base.
         */
template<typename T>
   T opt_window_wnaf_exp(const T &base, const number<gmp_int> &scalar, const size_t scalar_bits) {
      size_t best = 0;
      for (long i = T::wnaf_window_table.size() - 1; i >= 0; --i) {
         if (scalar_bits >= T::wnaf_window_table[i]) {
            best = i + 1;
            break;
         }
      }

      if (best > 0) {
         return fixed_window_wnaf_exp(best, base, scalar);
      } else {
         return scalar * base;
      }
   }
//
// Native non-member operations:
//
enum multi_exp_method {
   /**
   * Naive multi-exponentiation individually multiplies each base by the
   * corresponding scalar and adds up the results.
   * multi_exp_method_naive uses opt_window_wnaf_exp for exponentiation,
   * while multi_exp_method_plain uses operator *.
   */
   multi_exp_method_naive,
   multi_exp_method_naive_plain,
   /**
   * A variant of the Bos-Coster algorithm [1],
   * with implementation suggestions from [2].
   *
   * [1] = Bos and Coster, "Addition chain heuristics", CRYPTO '89
   * [2] = Bernstein, Duif, Lange, Schwabe, and Yang, "High-speed high-security signatures", CHES '11
   */
   multi_exp_method_bos_coster,
   /**
   * A special case of Pippenger's algorithm from Page 15 of
   * Bernstein, Doumen, Lange, Oosterwijk,
   * "Faster batch forgery identification", INDOCRYPT 2012
   * (https://eprint.iacr.org/2012/549.pdf)
   * When compiled with USE_MIXED_ADDITION, assumes input is in special form.
   * Requires that T implements .dbl() (and, if USE_MIXED_ADDITION is defined,
   * .to_special(), .mixed_add(), and batch_to_special()).
   */
   multi_exp_method_BDLO12
};

/**
* A window table stores window sizes for different instance sizes for fixed-base multi-scalar
* multiplications.
*/
template<typename T>
using window_table = std::vector<std::vector<T>>;

template<typename number<Backend>>
class ordered_exponent {
// to use std::push_heap and friends later
typedef number NumberType;

   public:
      size_t idx;
      NumberType r;

      ordered_exponent(const size_t idx, const NumberType &r) : idx(idx), r(r) {};

      bool operator<(const ordered_exponent<n> &other) const;
   };

   /**
   * multi_exp_inner<T, OT, Method>() implementes the specified
   * multiexponentiation method.
   * this implementation relies on some rather arcane template magic:
   * function templates cannot be partially specialized, so we cannot just write
   *     template<typename T, typename OT>
   *     T multi_exp_inner<T, OT, multi_exp_method_naive>
   * thus we resort to using std::enable_if. the basic idea is that *overloading*
   * is what's actually happening here, it's just that, for any given value of
   * Method, only one of the templates will be valid, and thus the correct
   * implementation will be used.
   */

   template<typename T, typename OT, multi_exp_method Method,
            typename std::enable_if<(Method == multi_exp_method_naive), int>::type = 0>
   T multi_exp_inner(typename std::vector<T>::const_iterator vec_start,
                     typename std::vector<T>::const_iterator vec_end,
                     typename std::vector<OT>::const_iterator scalar_start,
                     typename std::vector<OT>::const_iterator scalar_end) {
            T result(T::zero());

      typename std::vector<T>::const_iterator vec_it;
      typename std::vector<OT>::const_iterator scalar_it;

      for (vec_it = vec_start, scalar_it = scalar_start; vec_it != vec_end; ++vec_it, ++scalar_it) {
         NumberType scalar_bigint = scalar_it->as_bigint();
         result = result + opt_window_wnaf_exp(*vec_it, scalar_bigint, scalar_bigint.num_bits());
      }
      assert(scalar_it == scalar_end);

      return result;
   }

   template<typename T, typename OT, multi_exp_method Method,
            typename std::enable_if<(Method == multi_exp_method_naive_plain), int>::type = 0>
   T multi_exp_inner(typename std::vector<T>::const_iterator vec_start,
                     typename std::vector<T>::const_iterator vec_end,
                     typename std::vector<OT>::const_iterator scalar_start,
                     typename std::vector<OT>::const_iterator scalar_end) {
      T result(T::zero());

      typename std::vector<T>::const_iterator vec_it;
      typename std::vector<OT>::const_iterator scalar_it;

      for (vec_it = vec_start, scalar_it = scalar_start; vec_it != vec_end; ++vec_it, ++scalar_it) {
         result = result + (*scalar_it) * (*vec_it);
      }
      assert(scalar_it == scalar_end);

      return result;
   }

   template<typename T, typename OT, multi_exp_method Method,
            typename std::enable_if<(Method == multi_exp_method_BDLO12), int>::type = 0>
   T multi_exp_inner(typename std::vector<T>::const_iterator bases,
                     typename std::vector<T>::const_iterator bases_end,
                     typename std::vector<OT>::const_iterator exponents,
                     typename std::vector<OT>::const_iterator exponents_end) {
      UNUSED(exponents_end);
      size_t length = bases_end - bases;

      // empirically, this seems to be a decent estimate of the optimal value of c
      size_t log2_length = get_constant_ln2(length);
      size_t c = log2_length - (log2_length / 3 - 2);

      const mp_size_t exp_num_limbs = std::remove_reference<decltype(*exponents)>::type::num_limbs;
      std::vector<number_type<exp_num_limbs>> bn_exponents(length);
      size_t num_bits = 0;

      for (size_t i = 0; i < length; i++) {
         bn_exponents[i] = exponents[i].as_bigint();
         num_bits = std::max(num_bits, bn_exponents[i].num_bits());
      }

      size_t num_groups = (num_bits + c - 1) / c;

      T result;
      bool result_nonzero = false;

      for (size_t k = num_groups - 1; k <= num_groups; k--) {
         if (result_nonzero) {
            for (size_t i = 0; i < c; i++) {
               result = result.dbl();
            }
         }

         std::vector<T> buckets(1 << c);
         std::vector<bool> bucket_nonzero(1 << c);

         for (size_t i = 0; i < length; i++) {
            size_t id = 0;
            for (size_t j = 0; j < c; j++) {
               if (bn_exponents[i].test_bit(k * c + j)) {
                  id |= 1 << j;
               }
            }

            if (id == 0) {
               continue;
            }

            if (bucket_nonzero[id]) {
#ifdef USE_MIXED_ADDITION
               buckets[id] = buckets[id].mixed_add(bases[i]);
#else
               buckets[id] = buckets[id] + bases[i];
#endif
               } else {
                  buckets[id] = bases[i];
                  bucket_nonzero[id] = true;
               }
            }

#ifdef USE_MIXED_ADDITION
            batch_to_special(buckets);
#endif

            T running_sum;
            bool running_sum_nonzero = false;

            for (size_t i = (1u << c) - 1; i > 0; i--) {
               if (bucket_nonzero[i]) {
                  if (running_sum_nonzero) {
#ifdef USE_MIXED_ADDITION
                     running_sum = running_sum.mixed_add(buckets[i]);
#else
                     running_sum = running_sum + buckets[i];
#endif
                  } else {
                     running_sum = buckets[i];
                     running_sum_nonzero = true;
                  }
               }

               if (running_sum_nonzero) {
                  if (result_nonzero) {
                     result = result + running_sum;
                  } else {
                     result = running_sum;
                     result_nonzero = true;
                  }
               }
            }
         }

      return result;
   }

   template<typename T, typename OT, multi_exp_method Method,
            typename std::enable_if<(Method == multi_exp_method_bos_coster), int>::type = 0>
   T multi_exp_inner(typename std::vector<T>::const_iterator vec_start,
                     typename std::vector<T>::const_iterator vec_end,
                     typename std::vector<OT>::const_iterator scalar_start,
                     typename std::vector<OT>::const_iterator scalar_end) {
         const mp_size_t n = std::remove_reference<decltype(*scalar_start)>::type::num_limbs;

         if (vec_start == vec_end) {
            return T::zero();
         }

         if (vec_start + 1 == vec_end) {
            return (*scalar_start) * (*vec_start);
         }

         std::vector<ordered_exponent<n>> opt_q;
         const size_t vec_len = scalar_end - scalar_start;
         const size_t odd_vec_len = (vec_len % 2 == 1 ? vec_len : vec_len + 1);
         opt_q.reserve(odd_vec_len);
         std::vector<T> g;
         g.reserve(odd_vec_len);

         typename std::vector<T>::const_iterator vec_it;
         typename std::vector<OT>::const_iterator scalar_it;
         size_t i;
         for (i = 0, vec_it = vec_start, scalar_it = scalar_start; vec_it != vec_end; ++vec_it, ++scalar_it, ++i) {
            g.emplace_back(*vec_it);

            opt_q.emplace_back(ordered_exponent<n>(i, scalar_it->as_bigint()));
         }
         std::make_heap(opt_q.begin(), opt_q.end());
         assert(scalar_it == scalar_end);

         if (vec_len != odd_vec_len) {
            g.emplace_back(T::zero());
            opt_q.emplace_back(ordered_exponent<n>(odd_vec_len - 1, NumberType 0ul);
         }
         assert(g.size() % 2 == 1);
         assert(opt_q.size() == g.size());

         T opt_result = T::zero();

         while (true) {
            ordered_exponent<n> &a = opt_q[0];
            ordered_exponent<n> &b = (opt_q[1] < opt_q[2] ? opt_q[2] : opt_q[1]);

            const size_t abits = a.r.num_bits();

            if (b.r.is_zero()) {
               // opt_result = opt_result + (a.r * g[a.idx]);
               opt_result = opt_result + opt_window_wnaf_exp(g[a.idx], a.r, abits);
               break;
            }

            const size_t bbits = b.r.num_bits();
            const size_t limit = (abits - bbits >= 20 ? 20 : abits - bbits);

            if (bbits < 1ul << limit) {
               /*
               In this case, exponentiating to the power of a is cheaper than
               subtracting b from a multiple times, so let's do it directly
               */
               // opt_result = opt_result + (a.r * g[a.idx]);
               opt_result = opt_result + opt_window_wnaf_exp(g[a.idx], a.r, abits);
#ifdef DEBUG
               printf("Skipping the following pair (%zu bit number vs %zu bit):\n", abits, bbits);
               a.r.print();
               b.r.print();
#endif
               a.r.clear();
            } else {
               // x A + y B => (x-y) A + y (B+A)
               mpn_sub_n(a.r.data, a.r.data, b.r.data, n);
               g[b.idx] = g[b.idx] + g[a.idx];
            }

            // regardless of whether a was cleared or subtracted from we push it down, then take back up

            /* heapify A down */
            size_t a_pos = 0;
            while (2 * a_pos + 2 < odd_vec_len) {
               // this is a max-heap so to maintain a heap property we swap with the largest of the two
               if (opt_q[2 * a_pos + 1] < opt_q[2 * a_pos + 2]) {
                  std::swap(opt_q[a_pos], opt_q[2 * a_pos + 2]);
                  a_pos = 2 * a_pos + 2;
               } else {
                  std::swap(opt_q[a_pos], opt_q[2 * a_pos + 1]);
                  a_pos = 2 * a_pos + 1;
                  }
            }

            /* now heapify A up appropriate amount of times */
            while (a_pos > 0 && opt_q[(a_pos - 1) / 2] < opt_q[a_pos]) {
               std::swap(opt_q[a_pos], opt_q[(a_pos - 1) / 2]);
               a_pos = (a_pos - 1) / 2;
            }
         }

      return opt_result;
   }

   /**
   * Computes the sum
   * \sum_i scalar_start[i] * vec_start[i]
   * using the selected method.
   * Input is split into the given number of chunks, and, when compiled with
   * MULTICORE, the chunks are processed in parallel.
   */
   template<typename T, typename OT, multi_exp_method Method>
   T multi_exp(typename std::vector<T>::const_iterator vec_start,
               typename std::vector<T>::const_iterator vec_end,
               typename std::vector<OT>::const_iterator scalar_start,
               typename std::vector<OT>::const_iterator scalar_end,
               const size_t chunks) {
      const size_t total = vec_end - vec_start;
      if ((total < chunks) || (chunks == 1)) {
         // no need to split into "chunks", can call implementation directly
         return multi_exp_inner<T, OT, Method>(vec_start, vec_end, scalar_start, scalar_end);
      }

      const size_t one = total / chunks;

      std::vector<T> partial(chunks, T::zero());

#ifdef MULTICORE
#pragma omp parallel for
#endif
      for (size_t i = 0; i < chunks; ++i) {
         partial[i] =
         multi_exp_inner<T, OT, Method>(vec_start + i * one,
                                             (i == chunks - 1 ? vec_end : vec_start + (i + 1) * one),
                                             scalar_start + i * one,
                                             (i == chunks - 1 ? scalar_end : scalar_start + (i + 1) * one));
         }

         T final = T::zero();

         for (size_t i = 0; i < chunks; ++i) {
            final = final + partial[i];
         }

      return final;
   }

   /**
   * A variant of multi_exp that takes advantage of the method mixed_add (instead
   * of the operator '+').
   * Assumes input is in special form, and includes special pre-processing for
   * scalars equal to 0 or 1.
   */
   template<typename T, typename OT, multi_exp_method Method>
   T multi_exp_with_mixed_addition(typename std::vector<T>::const_iterator vec_start,
                                    typename std::vector<T>::const_iterator vec_end,
                                    typename std::vector<OT>::const_iterator scalar_start,
                                    typename std::vector<OT>::const_iterator scalar_end,
                                    const size_t chunks) {
      assert(std::distance(vec_start, vec_end) == std::distance(scalar_start, scalar_end));
      enter_block("Process scalar vector");
      auto value_it = vec_start;
      auto scalar_it = scalar_start;

      const OT zero = OT::zero();
      const OT one = OT::one();
      std::vector<OT> p;
      std::vector<T> g;

      T acc = T::zero();

      size_t num_skip = 0;
      size_t num_add = 0;
      size_t num_other = 0;

      for (; scalar_it != scalar_end; ++scalar_it, ++value_it) {
            if (*scalar_it == zero) {
               // do nothing
               ++num_skip;
            } else if (*scalar_it == one) {
#ifdef USE_MIXED_ADDITION
               acc = acc.mixed_add(*value_it);
#else
               acc = acc + (*value_it);
#endif
               ++num_add;
            } else {
               p.emplace_back(*scalar_it);
               g.emplace_back(*value_it);
               ++num_other;
            }
         }
      print_indent();
      printf("* Elements of w skipped: %zu (%0.2f%%)\n", num_skip,
               100. * num_skip / (num_skip + num_add + num_other));
      print_indent();
      printf("* Elements of w processed with special addition: %zu (%0.2f%%)\n", num_add,
               100. * num_add / (num_skip + num_add + num_other));
      print_indent();
      printf("* Elements of w remaining: %zu (%0.2f%%)\n", num_other,
               100. * num_other / (num_skip + num_add + num_other));

      leave_block("Process scalar vector");

      return acc + multi_exp<T, OT, Method>(g.begin(), g.end(), p.begin(), p.end(), chunks);
   }

   /**
   * A convenience function for calculating a pure inner product, where the
   * more complicated methods are not required.
   */
   template<typename InputIterator>
   typename std::iterator_traits<InputIterator>::value_type
      inner_product(InputIterator a_start, InputIterator a_end, InputIterator b_start, InputIterator b_end) {
      return multi_exp<typename std::iterator_traits<InputIterator>::value_type,
                        typename std::iterator_traits<InputIterator>::value_type, multi_exp_method_naive_plain>(
      a_start, a_end, b_start, b_end, 1);
   }

        /**
         * Compute window size for the given number of scalars.
         */
   template<typename T>
   size_t get_exp_window_size(const size_t num_scalars) {
      if (T::fixed_base_exp_window_table.empty()) {
#ifdef LOWMEM
      return 14;
#else
      return 17;
#endif
      }
      size_t window = 1;
      for (long i = T::fixed_base_exp_window_table.size() - 1; i >= 0; --i) {
#ifdef DEBUG
      if (!inhibit_profiling_info) {
         printf("%ld %zu %zu\n", i, num_scalars, T::fixed_base_exp_window_table[i]);
      }
#endif
      if (T::fixed_base_exp_window_table[i] != 0 && num_scalars >= T::fixed_base_exp_window_table[i]) {
         window = i + 1;
         break;
         }
      }

      if (!inhibit_profiling_info) {
         print_indent();
         printf("Choosing window size %zu for %zu elements\n", window, num_scalars);
      }

#ifdef LOWMEM
      window = std::min((size_t)14, window);
#endif
      return window;
   }

   /**
   * Compute table of window sizes.
   */
   template<typename T>
   window_table<T> get_window_table(const size_t scalar_size, const size_t window, const T &g) {
      const size_t in_window = 1ul << window;
      const size_t outerc = (scalar_size + window - 1) / window;
      const size_t last_in_window = 1ul << (scalar_size - (outerc - 1) * window);
#ifdef DEBUG
      if (!inhibit_profiling_info) {
         print_indent();
         printf("* scalar_size=%zu; window=%zu; in_window=%zu; outerc=%zu\n", scalar_size, window, in_window,
               outerc);
      }
#endif

      window_table<T> powers_of_g(outerc, std::vector<T>(in_window, T::zero()));

      T gouter = g;

      for (size_t outer = 0; outer < outerc; ++outer) {
         T ginner = T::zero();
         size_t cur_in_window = outer == outerc - 1 ? last_in_window : in_window;
         for (size_t inner = 0; inner < cur_in_window; ++inner) {
            powers_of_g[outer][inner] = ginner;
            ginner = ginner + gouter;
         }

         for (size_t i = 0; i < window; ++i) {
            gouter = gouter + gouter;
            }
         }

      return powers_of_g;
   }

   template<typename T, typename OT>
   T windowed_exp(const size_t scalar_size,
                  const size_t window,
                  const window_table<T> &powers_of_g,
                  const OT &pow) {
      const size_t outerc = (scalar_size + window - 1) / window;
      const NumberType pow_val = pow.as_bigint();

      /* exp */
      T res = powers_of_g[0][0];

      for (size_t outer = 0; outer < outerc; ++outer) {
            size_t inner = 0;
            for (size_t i = 0; i < window; ++i) {
               if (pow_val.test_bit(outer * window + i)) {
                  inner |= 1u << i;
               }
            }

         res = res + powers_of_g[outer][inner];
         }

      return res;
   }

   template<typename T, typename OT>
   std::vector<T> batch_exp(const size_t scalar_size,
                            const size_t window,
                            const window_table<T> &table,
                            const std::vector<OT> &v) {
      if (!inhibit_profiling_info) {
         print_indent();
      }
      std::vector<T> res(v.size(), table[0][0]);

#ifdef MULTICORE
#pragma omp parallel for
#endif
      for (size_t i = 0; i < v.size(); ++i) {
         res[i] = windowed_exp(scalar_size, window, table, v[i]);

         if (!inhibit_profiling_info && (i % 10000 == 0)) {
            printf(".");
            fflush(stdout);
            }
         }

         if (!inhibit_profiling_info) {
            printf(" DONE!\n");
         }

      return res;
   }

   template<typename T, typename OT>
   std::vector<T> batch_exp_with_coeff(const size_t scalar_size,
                                       const size_t window,
                                       const window_table<T> &table,
                                       const OT &coeff,
                                       const std::vector<OT> &v) {
      if (!inhibit_profiling_info) {
         print_indent();
      }
      std::vector<T> res(v.size(), table[0][0]);

#ifdef MULTICORE
#pragma omp parallel for
#endif
      for (size_t i = 0; i < v.size(); ++i) {
            res[i] = windowed_exp(scalar_size, window, table, coeff * v[i]);

            if (!inhibit_profiling_info && (i % 10000 == 0)) {
               printf(".");
               fflush(stdout);
            }
         }

      if (!inhibit_profiling_info) {
         printf(" DONE!\n");
      }

      return res;
   }

   template<typename T>
   void batch_to_special(std::vector<T> &vec) {
      enter_block("Batch-convert elements to special form");

      std::vector<T> non_zero_vec;
      for (size_t i = 0; i < vec.size(); ++i) {
         if (!vec[i].is_zero()) {
               non_zero_vec.emplace_back(vec[i]);
            }
         }

      T::batch_to_special_all_non_zeros(non_zero_vec);
      auto it = non_zero_vec.begin();
      T zero_special = T::zero();
      zero_special.to_special();

      for (size_t i = 0; i < vec.size(); ++i) {
         if (!vec[i].is_zero()) {
            vec[i] = *it;
             ++it;
         } else {
            vec[i] = zero_special;
         }
      }
      leave_block("Batch-convert elements to special form");
   }

template <unsigned Digits10>
inline void eval_sqrt(gmp_float<Digits10>& result, const gmp_float<Digits10>& val)
{
   mpf_sqrt(result.data(), val.data());
}

template <unsigned Digits10>
inline void eval_abs(gmp_float<Digits10>& result, const gmp_float<Digits10>& val)
{
   mpf_abs(result.data(), val.data());
}

template <unsigned Digits10>
inline void eval_fabs(gmp_float<Digits10>& result, const gmp_float<Digits10>& val)
{
   mpf_abs(result.data(), val.data());
}
template <unsigned Digits10>
inline void eval_ceil(gmp_float<Digits10>& result, const gmp_float<Digits10>& val)
{
   mpf_ceil(result.data(), val.data());
}
template <unsigned Digits10>
inline void eval_floor(gmp_float<Digits10>& result, const gmp_float<Digits10>& val)
{
   mpf_floor(result.data(), val.data());
}
template <unsigned Digits10>
inline void eval_trunc(gmp_float<Digits10>& result, const gmp_float<Digits10>& val)
{
   mpf_trunc(result.data(), val.data());
}
template <unsigned Digits10>
inline void eval_ldexp(gmp_float<Digits10>& result, const gmp_float<Digits10>& val, long e)
{
   if (e > 0)
      mpf_mul_2exp(result.data(), val.data(), e);
   else if (e < 0)
      mpf_div_2exp(result.data(), val.data(), -e);
   else
      result = val;
}
template <unsigned Digits10>
inline void eval_frexp(gmp_float<Digits10>& result, const gmp_float<Digits10>& val, int* e)
{
#if (BOOST_MP_MPIR_VERSION >= 20600) && (BOOST_MP_MPIR_VERSION < 30000)
   mpir_si v;
   mpf_get_d_2exp(&v, val.data());
#else
   long                             v;
   mpf_get_d_2exp(&v, val.data());
#endif
   *e = v;
   eval_ldexp(result, val, -v);
}
template <unsigned Digits10>
inline void eval_frexp(gmp_float<Digits10>& result, const gmp_float<Digits10>& val, long* e)
{
#if (BOOST_MP_MPIR_VERSION >= 20600) && (BOOST_MP_MPIR_VERSION < 30000)
   mpir_si v;
   mpf_get_d_2exp(&v, val.data());
   *e = v;
   eval_ldexp(result, val, -v);
#else
   mpf_get_d_2exp(e, val.data());
   eval_ldexp(result, val, -*e);
#endif
}

template <unsigned Digits10>
inline std::size_t hash_value(const gmp_float<Digits10>& val)
{
   std::size_t result = 0;
   for (int i = 0; i < std::abs(val.data()[0]._mp_size); ++i)
      boost::hash_combine(result, val.data()[0]._mp_d[i]);
   boost::hash_combine(result, val.data()[0]._mp_exp);
   boost::hash_combine(result, val.data()[0]._mp_size);
   return result;
}

struct gmp_int
{
#ifdef BOOST_HAS_LONG_LONG
   typedef mpl::list<long, boost::long_long_type>           signed_types;
   typedef mpl::list<unsigned long, boost::ulong_long_type> unsigned_types;
#else
   typedef mpl::list<long>          signed_types;
   typedef mpl::list<unsigned long> unsigned_types;
#endif
   typedef mpl::list<double, long double> float_types;

   gmp_int()
   {
      mpz_init(this->m_data);
   }
   gmp_int(const gmp_int& o)
   {
      if (o.m_data[0]._mp_d)
         mpz_init_set(m_data, o.m_data);
      else
         mpz_init(this->m_data);
   }
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
   gmp_int(gmp_int&& o) BOOST_NOEXCEPT
   {
      m_data[0]         = o.m_data[0];
      o.m_data[0]._mp_d = 0;
   }
#endif
   explicit gmp_int(const mpf_t val)
   {
      mpz_init(this->m_data);
      mpz_set_f(this->m_data, val);
   }
   gmp_int(const mpz_t val)
   {
      mpz_init_set(this->m_data, val);
   }
   explicit gmp_int(const mpq_t val)
   {
      mpz_init(this->m_data);
      mpz_set_q(this->m_data, val);
   }
   template <unsigned Digits10>
   explicit gmp_int(const gmp_float<Digits10>& o)
   {
      mpz_init(this->m_data);
      mpz_set_f(this->m_data, o.data());
   }
   explicit gmp_int(const gmp_rational& o);
   gmp_int& operator=(const gmp_int& o)
   {
      if (m_data[0]._mp_d == 0)
         mpz_init(this->m_data);
      mpz_set(m_data, o.m_data);
      return *this;
   }
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
   gmp_int& operator=(gmp_int&& o) BOOST_NOEXCEPT
   {
      mpz_swap(m_data, o.m_data);
      return *this;
   }
#endif
#ifdef BOOST_HAS_LONG_LONG
#if defined(ULLONG_MAX) && (ULLONG_MAX == ULONG_MAX)
   gmp_int& operator=(boost::ulong_long_type i)
   {
      *this = static_cast<unsigned long>(i);
      return *this;
   }
#else
   gmp_int& operator=(boost::ulong_long_type i)
   {
      if (m_data[0]._mp_d == 0)
         mpz_init(this->m_data);
      boost::ulong_long_type mask  = ((((1uLL << (std::numeric_limits<unsigned long>::digits - 1)) - 1) << 1) | 1uLL);
      unsigned               shift = 0;
      mpz_t                  t;
      mpz_set_ui(m_data, 0);
      mpz_init_set_ui(t, 0);
      while (i)
      {
         mpz_set_ui(t, static_cast<unsigned long>(i & mask));
         if (shift)
            mpz_mul_2exp(t, t, shift);
         mpz_add(m_data, m_data, t);
         shift += std::numeric_limits<unsigned long>::digits;
         i >>= std::numeric_limits<unsigned long>::digits;
      }
      mpz_clear(t);
      return *this;
   }
#endif
   gmp_int& operator=(boost::long_long_type i)
   {
      if (m_data[0]._mp_d == 0)
         mpz_init(this->m_data);
      bool neg = i < 0;
      *this    = boost::multiprecision::detail::unsigned_abs(i);
      if (neg)
         mpz_neg(m_data, m_data);
      return *this;
   }
#endif
#ifdef BOOST_HAS_INT128
   gmp_int& operator=(unsigned __int128 i)
   {
      if (m_data[0]._mp_d == 0)
         mpz_init(this->m_data);
      unsigned __int128 mask  = ((((1uLL << (std::numeric_limits<unsigned long>::digits - 1)) - 1) << 1) | 1uLL);
      unsigned               shift = 0;
      mpz_t                  t;
      mpz_set_ui(m_data, 0);
      mpz_init_set_ui(t, 0);
      while (i)
      {
         mpz_set_ui(t, static_cast<unsigned long>(i & mask));
         if (shift)
            mpz_mul_2exp(t, t, shift);
         mpz_add(m_data, m_data, t);
         shift += std::numeric_limits<unsigned long>::digits;
         i >>= std::numeric_limits<unsigned long>::digits;
      }
      mpz_clear(t);
      return *this;
   }
   gmp_int& operator=(__int128 i)
   {
      if (m_data[0]._mp_d == 0)
         mpz_init(this->m_data);
      bool neg = i < 0;
      *this    = boost::multiprecision::detail::unsigned_abs(i);
      if (neg)
         mpz_neg(m_data, m_data);
      return *this;
   }
#endif
   gmp_int& operator=(unsigned long i)
   {
      if (m_data[0]._mp_d == 0)
         mpz_init(this->m_data);
      mpz_set_ui(m_data, i);
      return *this;
   }
   gmp_int& operator=(long i)
   {
      if (m_data[0]._mp_d == 0)
         mpz_init(this->m_data);
      mpz_set_si(m_data, i);
      return *this;
   }
   gmp_int& operator=(double d)
   {
      if (m_data[0]._mp_d == 0)
         mpz_init(this->m_data);
      mpz_set_d(m_data, d);
      return *this;
   }
   gmp_int& operator=(long double a)
   {
      using std::floor;
      using std::frexp;
      using std::ldexp;

      if (m_data[0]._mp_d == 0)
         mpz_init(this->m_data);

      if (a == 0)
      {
         mpz_set_si(m_data, 0);
         return *this;
      }

      if (a == 1)
      {
         mpz_set_si(m_data, 1);
         return *this;
      }

      BOOST_ASSERT(!(boost::math::isinf)(a));
      BOOST_ASSERT(!(boost::math::isnan)(a));

      int         e;
      long double f, term;
      mpz_set_ui(m_data, 0u);

      f = frexp(a, &e);

      static const int shift = std::numeric_limits<int>::digits - 1;

      while (f)
      {
         // extract int sized bits from f:
         f    = ldexp(f, shift);
         term = floor(f);
         e -= shift;
         mpz_mul_2exp(m_data, m_data, shift);
         if (term > 0)
            mpz_add_ui(m_data, m_data, static_cast<unsigned>(term));
         else
            mpz_sub_ui(m_data, m_data, static_cast<unsigned>(-term));
         f -= term;
      }
      if (e > 0)
         mpz_mul_2exp(m_data, m_data, e);
      else if (e < 0)
         mpz_div_2exp(m_data, m_data, -e);
      return *this;
   }
   gmp_int& operator=(const char* s)
   {
      if (m_data[0]._mp_d == 0)
         mpz_init(this->m_data);
      std::size_t n     = s ? std::strlen(s) : 0;
      int         radix = 10;
      if (n && (*s == '0'))
      {
         if ((n > 1) && ((s[1] == 'x') || (s[1] == 'X')))
         {
            radix = 16;
            s += 2;
            n -= 2;
         }
         else
         {
            radix = 8;
            n -= 1;
         }
      }
      if (n)
      {
         if (0 != mpz_set_str(m_data, s, radix))
            BOOST_THROW_EXCEPTION(std::runtime_error(std::string("The string \"") + s + std::string("\"could not be interpreted as a valid integer.")));
      }
      else
         mpz_set_ui(m_data, 0);
      return *this;
   }
   gmp_int& operator=(const mpf_t val)
   {
      if (m_data[0]._mp_d == 0)
         mpz_init(this->m_data);
      mpz_set_f(this->m_data, val);
      return *this;
   }
   gmp_int& operator=(const mpz_t val)
   {
      if (m_data[0]._mp_d == 0)
         mpz_init(this->m_data);
      mpz_set(this->m_data, val);
      return *this;
   }
   gmp_int& operator=(const mpq_t val)
   {
      if (m_data[0]._mp_d == 0)
         mpz_init(this->m_data);
      mpz_set_q(this->m_data, val);
      return *this;
   }
   template <unsigned Digits10>
   gmp_int& operator=(const gmp_float<Digits10>& o)
   {
      if (m_data[0]._mp_d == 0)
         mpz_init(this->m_data);
      mpz_set_f(this->m_data, o.data());
      return *this;
   }
   gmp_int& operator=(const gmp_rational& o);
   void     swap(gmp_int& o)
   {
      mpz_swap(m_data, o.m_data);
   }
   std::string str(std::streamsize /*digits*/, std::ios_base::fmtflags f) const
   {
      BOOST_ASSERT(m_data[0]._mp_d);

      int base = 10;
      if ((f & std::ios_base::oct) == std::ios_base::oct)
         base = 8;
      else if ((f & std::ios_base::hex) == std::ios_base::hex)
         base = 16;
      //
      // sanity check, bases 8 and 16 are only available for positive numbers:
      //
      if ((base != 10) && (mpz_sgn(m_data) < 0))
         BOOST_THROW_EXCEPTION(std::runtime_error("Formatted output in bases 8 or 16 is only available for positive numbers"));
      void* (*alloc_func_ptr)(size_t);
      void* (*realloc_func_ptr)(void*, size_t, size_t);
      void (*free_func_ptr)(void*, size_t);
      const char* ps = mpz_get_str(0, base, m_data);
      std::string s  = ps;
      mp_get_memory_functions(&alloc_func_ptr, &realloc_func_ptr, &free_func_ptr);
      (*free_func_ptr)((void*)ps, std::strlen(ps) + 1);
      if (f & std::ios_base::uppercase)
         for (size_t i = 0; i < s.length(); ++i)
            s[i] = std::toupper(s[i]);
      if ((base != 10) && (f & std::ios_base::showbase))
      {
         int         pos = s[0] == '-' ? 1 : 0;
         const char* pp  = base == 8 ? "0" : (f & std::ios_base::uppercase) ? "0X" : "0x";
         s.insert(static_cast<std::string::size_type>(pos), pp);
      }
      if ((f & std::ios_base::showpos) && (s[0] != '-'))
         s.insert(static_cast<std::string::size_type>(0), 1, '+');

      return s;
   }
   ~gmp_int() BOOST_NOEXCEPT
   {
      if (m_data[0]._mp_d)
         mpz_clear(m_data);
   }
   void negate() BOOST_NOEXCEPT
   {
      BOOST_ASSERT(m_data[0]._mp_d);
      mpz_neg(m_data, m_data);
   }
   int compare(const gmp_int& o) const BOOST_NOEXCEPT
   {
      BOOST_ASSERT(m_data[0]._mp_d && o.m_data[0]._mp_d);
      return mpz_cmp(m_data, o.m_data);
   }
   int compare(long i) const BOOST_NOEXCEPT
   {
      BOOST_ASSERT(m_data[0]._mp_d);
      return mpz_cmp_si(m_data, i);
   }
   int compare(unsigned long i) const BOOST_NOEXCEPT
   {
      BOOST_ASSERT(m_data[0]._mp_d);
      return mpz_cmp_ui(m_data, i);
   }
   template <class V>
   int compare(V v) const
   {
      gmp_int d;
      d = v;
      return compare(d);
   }
   mpz_t& data() BOOST_NOEXCEPT
   {
      BOOST_ASSERT(m_data[0]._mp_d);
      return m_data;
   }
   const mpz_t& data() const BOOST_NOEXCEPT
   {
      BOOST_ASSERT(m_data[0]._mp_d);
      return m_data;
   }

 protected:
   mpz_t m_data;
};

template <class T>
inline typename enable_if<is_arithmetic<T>, bool>::type eval_eq(const gmp_int& a, const T& b)
{
   return a.compare(b) == 0;
}
template <class T>
inline typename enable_if<is_arithmetic<T>, bool>::type eval_lt(const gmp_int& a, const T& b)
{
   return a.compare(b) < 0;
}
template <class T>
inline typename enable_if<is_arithmetic<T>, bool>::type eval_gt(const gmp_int& a, const T& b)
{
   return a.compare(b) > 0;
}

inline bool eval_is_zero(const gmp_int& val)
{
   return mpz_sgn(val.data()) == 0;
}
inline void eval_add(gmp_int& t, const gmp_int& o)
{
   mpz_add(t.data(), t.data(), o.data());
}
inline void eval_multiply_add(gmp_int& t, const gmp_int& a, const gmp_int& b)
{
   mpz_addmul(t.data(), a.data(), b.data());
}
inline void eval_multiply_subtract(gmp_int& t, const gmp_int& a, const gmp_int& b)
{
   mpz_submul(t.data(), a.data(), b.data());
}
inline void eval_subtract(gmp_int& t, const gmp_int& o)
{
   mpz_sub(t.data(), t.data(), o.data());
}
inline void eval_multiply(gmp_int& t, const gmp_int& o)
{
   mpz_mul(t.data(), t.data(), o.data());
}
inline void eval_divide(gmp_int& t, const gmp_int& o)
{
   if (eval_is_zero(o))
      BOOST_THROW_EXCEPTION(std::overflow_error("Division by zero."));
   mpz_tdiv_q(t.data(), t.data(), o.data());
}
inline void eval_modulus(gmp_int& t, const gmp_int& o)
{
   mpz_tdiv_r(t.data(), t.data(), o.data());
}
inline void eval_add(gmp_int& t, unsigned long i)
{
   mpz_add_ui(t.data(), t.data(), i);
}
inline void eval_multiply_add(gmp_int& t, const gmp_int& a, unsigned long i)
{
   mpz_addmul_ui(t.data(), a.data(), i);
}
inline void eval_multiply_subtract(gmp_int& t, const gmp_int& a, unsigned long i)
{
   mpz_submul_ui(t.data(), a.data(), i);
}
inline void eval_subtract(gmp_int& t, unsigned long i)
{
   mpz_sub_ui(t.data(), t.data(), i);
}
inline void eval_multiply(gmp_int& t, unsigned long i)
{
   mpz_mul_ui(t.data(), t.data(), i);
}
inline void eval_modulus(gmp_int& t, unsigned long i)
{
   mpz_tdiv_r_ui(t.data(), t.data(), i);
}
inline void eval_divide(gmp_int& t, unsigned long i)
{
   if (i == 0)
      BOOST_THROW_EXCEPTION(std::overflow_error("Division by zero."));
   mpz_tdiv_q_ui(t.data(), t.data(), i);
}
inline void eval_add(gmp_int& t, long i)
{
   if (i > 0)
      mpz_add_ui(t.data(), t.data(), i);
   else
      mpz_sub_ui(t.data(), t.data(), boost::multiprecision::detail::unsigned_abs(i));
}
inline void eval_multiply_add(gmp_int& t, const gmp_int& a, long i)
{
   if (i > 0)
      mpz_addmul_ui(t.data(), a.data(), i);
   else
      mpz_submul_ui(t.data(), a.data(), boost::multiprecision::detail::unsigned_abs(i));
}
inline void eval_multiply_subtract(gmp_int& t, const gmp_int& a, long i)
{
   if (i > 0)
      mpz_submul_ui(t.data(), a.data(), i);
   else
      mpz_addmul_ui(t.data(), a.data(), boost::multiprecision::detail::unsigned_abs(i));
}
inline void eval_subtract(gmp_int& t, long i)
{
   if (i > 0)
      mpz_sub_ui(t.data(), t.data(), i);
   else
      mpz_add_ui(t.data(), t.data(), boost::multiprecision::detail::unsigned_abs(i));
}
inline void eval_multiply(gmp_int& t, long i)
{
   mpz_mul_ui(t.data(), t.data(), boost::multiprecision::detail::unsigned_abs(i));
   if (i < 0)
      mpz_neg(t.data(), t.data());
}
inline void eval_modulus(gmp_int& t, long i)
{
   mpz_tdiv_r_ui(t.data(), t.data(), boost::multiprecision::detail::unsigned_abs(i));
}
inline void eval_divide(gmp_int& t, long i)
{
   if (i == 0)
      BOOST_THROW_EXCEPTION(std::overflow_error("Division by zero."));
   mpz_tdiv_q_ui(t.data(), t.data(), boost::multiprecision::detail::unsigned_abs(i));
   if (i < 0)
      mpz_neg(t.data(), t.data());
}
template <class UI>
inline void eval_left_shift(gmp_int& t, UI i)
{
   mpz_mul_2exp(t.data(), t.data(), static_cast<unsigned long>(i));
}
template <class UI>
inline void eval_right_shift(gmp_int& t, UI i)
{
   mpz_fdiv_q_2exp(t.data(), t.data(), static_cast<unsigned long>(i));
}
template <class UI>
inline void eval_left_shift(gmp_int& t, const gmp_int& v, UI i)
{
   mpz_mul_2exp(t.data(), v.data(), static_cast<unsigned long>(i));
}
template <class UI>
inline void eval_right_shift(gmp_int& t, const gmp_int& v, UI i)
{
   mpz_fdiv_q_2exp(t.data(), v.data(), static_cast<unsigned long>(i));
}

inline void eval_bitwise_and(gmp_int& result, const gmp_int& v)
{
   mpz_and(result.data(), result.data(), v.data());
}

inline void eval_bitwise_or(gmp_int& result, const gmp_int& v)
{
   mpz_ior(result.data(), result.data(), v.data());
}

inline void eval_bitwise_xor(gmp_int& result, const gmp_int& v)
{
   mpz_xor(result.data(), result.data(), v.data());
}

inline void eval_add(gmp_int& t, const gmp_int& p, const gmp_int& o)
{
   mpz_add(t.data(), p.data(), o.data());
}
inline void eval_subtract(gmp_int& t, const gmp_int& p, const gmp_int& o)
{
   mpz_sub(t.data(), p.data(), o.data());
}
inline void eval_multiply(gmp_int& t, const gmp_int& p, const gmp_int& o)
{
   mpz_mul(t.data(), p.data(), o.data());
}
inline void eval_divide(gmp_int& t, const gmp_int& p, const gmp_int& o)
{
   if (eval_is_zero(o))
      BOOST_THROW_EXCEPTION(std::overflow_error("Division by zero."));
   mpz_tdiv_q(t.data(), p.data(), o.data());
}
inline void eval_modulus(gmp_int& t, const gmp_int& p, const gmp_int& o)
{
   mpz_tdiv_r(t.data(), p.data(), o.data());
}
inline void eval_add(gmp_int& t, const gmp_int& p, unsigned long i)
{
   mpz_add_ui(t.data(), p.data(), i);
}
inline void eval_subtract(gmp_int& t, const gmp_int& p, unsigned long i)
{
   mpz_sub_ui(t.data(), p.data(), i);
}
inline void eval_multiply(gmp_int& t, const gmp_int& p, unsigned long i)
{
   mpz_mul_ui(t.data(), p.data(), i);
}
inline void eval_modulus(gmp_int& t, const gmp_int& p, unsigned long i)
{
   mpz_tdiv_r_ui(t.data(), p.data(), i);
}
inline void eval_divide(gmp_int& t, const gmp_int& p, unsigned long i)
{
   if (i == 0)
      BOOST_THROW_EXCEPTION(std::overflow_error("Division by zero."));
   mpz_tdiv_q_ui(t.data(), p.data(), i);
}
inline void eval_add(gmp_int& t, const gmp_int& p, long i)
{
   if (i > 0)
      mpz_add_ui(t.data(), p.data(), i);
   else
      mpz_sub_ui(t.data(), p.data(), boost::multiprecision::detail::unsigned_abs(i));
}
inline void eval_subtract(gmp_int& t, const gmp_int& p, long i)
{
   if (i > 0)
      mpz_sub_ui(t.data(), p.data(), i);
   else
      mpz_add_ui(t.data(), p.data(), boost::multiprecision::detail::unsigned_abs(i));
}
inline void eval_multiply(gmp_int& t, const gmp_int& p, long i)
{
   mpz_mul_ui(t.data(), p.data(), boost::multiprecision::detail::unsigned_abs(i));
   if (i < 0)
      mpz_neg(t.data(), t.data());
}
inline void eval_modulus(gmp_int& t, const gmp_int& p, long i)
{
   mpz_tdiv_r_ui(t.data(), p.data(), boost::multiprecision::detail::unsigned_abs(i));
}
inline void eval_divide(gmp_int& t, const gmp_int& p, long i)
{
   if (i == 0)
      BOOST_THROW_EXCEPTION(std::overflow_error("Division by zero."));
   mpz_tdiv_q_ui(t.data(), p.data(), boost::multiprecision::detail::unsigned_abs(i));
   if (i < 0)
      mpz_neg(t.data(), t.data());
}

inline void eval_bitwise_and(gmp_int& result, const gmp_int& u, const gmp_int& v)
{
   mpz_and(result.data(), u.data(), v.data());
}

inline void eval_bitwise_or(gmp_int& result, const gmp_int& u, const gmp_int& v)
{
   mpz_ior(result.data(), u.data(), v.data());
}

inline void eval_bitwise_xor(gmp_int& result, const gmp_int& u, const gmp_int& v)
{
   mpz_xor(result.data(), u.data(), v.data());
}

inline void eval_complement(gmp_int& result, const gmp_int& u)
{
   mpz_com(result.data(), u.data());
}

inline int eval_get_sign(const gmp_int& val)
{
   return mpz_sgn(val.data());
}
inline void eval_convert_to(unsigned long* result, const gmp_int& val)
{
   if (mpz_sgn(val.data()) < 0)
   {
      BOOST_THROW_EXCEPTION(std::range_error("Conversion from negative integer to an unsigned type results in undefined behaviour"));
   }
   else
      *result = (unsigned long)mpz_get_ui(val.data());
}
inline void eval_convert_to(long* result, const gmp_int& val)
{
   if (0 == mpz_fits_slong_p(val.data()))
   {
      *result = mpz_sgn(val.data()) < 0 ? (std::numeric_limits<long>::min)() : (std::numeric_limits<long>::max)();
   }
   else
      *result = (signed long)mpz_get_si(val.data());
}
inline void eval_convert_to(double* result, const gmp_int& val)
{
   *result = mpz_get_d(val.data());
}
#ifdef BOOST_HAS_LONG_LONG
inline void eval_convert_to(boost::ulong_long_type* result, const gmp_int& val)
{
   if (mpz_sgn(val.data()) < 0)
   {
      BOOST_THROW_EXCEPTION(std::range_error("Conversion from negative integer to an unsigned type results in undefined behaviour"));
   }
   *result = 0;
   gmp_int t(val);
   unsigned parts = sizeof(boost::ulong_long_type) / sizeof(unsigned long);

   for (unsigned i = 0; i < parts; ++i)
   {
      boost::ulong_long_type part = mpz_get_ui(t.data());
      if (i)
         *result |= part << (i * sizeof(unsigned long) * CHAR_BIT);
      else
         *result = part;
      mpz_tdiv_q_2exp(t.data(), t.data(), sizeof(unsigned long) * CHAR_BIT);
   }
}
inline void eval_convert_to(boost::long_long_type* result, const gmp_int& val)
{
   int s = mpz_sgn(val.data());
   *result = 0;
   gmp_int t(val);
   unsigned parts = sizeof(boost::ulong_long_type) / sizeof(unsigned long);
   boost::ulong_long_type unsigned_result = 0;

   for (unsigned i = 0; i < parts; ++i)
   {
      boost::ulong_long_type part = mpz_get_ui(t.data());
      if (i)
         unsigned_result |= part << (i * sizeof(unsigned long) * CHAR_BIT);
      else
         unsigned_result = part;
      mpz_tdiv_q_2exp(t.data(), t.data(), sizeof(unsigned long) * CHAR_BIT);
   }
   //
   // Overflow check:
   //
   bool overflow = false;
   if (mpz_sgn(t.data()))
   {
      overflow = true;
   }
   if ((s > 0) && (unsigned_result > static_cast<boost::ulong_long_type>((std::numeric_limits<boost::long_long_type>::max)())))
      overflow = true;
   if((s < 0) && (unsigned_result > 1u - static_cast<boost::ulong_long_type>((std::numeric_limits<boost::long_long_type>::min)() + 1)))
      overflow = true;
   if(overflow)
      *result = s < 0 ? (std::numeric_limits<boost::long_long_type>::min)() : (std::numeric_limits<boost::long_long_type>::max)();
   else
      *result = s < 0 ? -boost::long_long_type(unsigned_result - 1) - 1 : unsigned_result;
}
#endif
#ifdef BOOST_HAS_INT128
inline void eval_convert_to(unsigned __int128* result, const gmp_int& val)
{
   if (mpz_sgn(val.data()) < 0)
   {
      BOOST_THROW_EXCEPTION(std::range_error("Conversion from negative integer to an unsigned type results in undefined behaviour"));
   }
   *result = 0;
   gmp_int t(val);
   unsigned parts = sizeof(unsigned __int128) / sizeof(unsigned long);

   for (unsigned i = 0; i < parts; ++i)
   {
      unsigned __int128 part = mpz_get_ui(t.data());
      if (i)
         *result |= part << (i * sizeof(unsigned long) * CHAR_BIT);
      else
         *result = part;
      mpz_tdiv_q_2exp(t.data(), t.data(), sizeof(unsigned long) * CHAR_BIT);
   }
}
inline void eval_convert_to(__int128* result, const gmp_int& val)
{
   int s = mpz_sgn(val.data());
   *result = 0;
   gmp_int t(val);
   unsigned parts = sizeof(unsigned __int128) / sizeof(unsigned long);
   unsigned __int128 unsigned_result = 0;

   for (unsigned i = 0; i < parts; ++i)
   {
      unsigned __int128 part = mpz_get_ui(t.data());
      if (i)
         unsigned_result |= part << (i * sizeof(unsigned long) * CHAR_BIT);
      else
         unsigned_result = part;
      mpz_tdiv_q_2exp(t.data(), t.data(), sizeof(unsigned long) * CHAR_BIT);
   }
   //
   // Overflow check:
   //
   static const __int128 int128_max = static_cast<__int128>((static_cast<unsigned __int128>(1u) << 127) - 1);
   static const __int128 int128_min = (static_cast<unsigned __int128>(1u) << 127);
   bool overflow = false;
   if (mpz_sgn(t.data()))
   {
      overflow = true;
   }
   if ((s > 0) && (unsigned_result > static_cast<unsigned __int128>(int128_max)))
      overflow = true;
   if ((s < 0) && (unsigned_result > 1u - static_cast<unsigned __int128>(int128_min + 1)))
      overflow = true;
   if (overflow)
      *result = s < 0 ? int128_min : int128_max;
   else
      *result = s < 0 ? -__int128(unsigned_result - 1) - 1 : unsigned_result;
}
#endif
inline void eval_abs(gmp_int& result, const gmp_int& val)
{
   mpz_abs(result.data(), val.data());
}

inline void eval_gcd(gmp_int& result, const gmp_int& a, const gmp_int& b)
{
   mpz_gcd(result.data(), a.data(), b.data());
}
inline void eval_lcm(gmp_int& result, const gmp_int& a, const gmp_int& b)
{
   mpz_lcm(result.data(), a.data(), b.data());
}
template <class I>
inline typename enable_if_c<(is_unsigned<I>::value && (sizeof(I) <= sizeof(unsigned long)))>::type eval_gcd(gmp_int& result, const gmp_int& a, const I b)
{
   mpz_gcd_ui(result.data(), a.data(), b);
}
template <class I>
inline typename enable_if_c<(is_unsigned<I>::value && (sizeof(I) <= sizeof(unsigned long)))>::type eval_lcm(gmp_int& result, const gmp_int& a, const I b)
{
   mpz_lcm_ui(result.data(), a.data(), b);
}
template <class I>
inline typename enable_if_c<(is_signed<I>::value && (sizeof(I) <= sizeof(long)))>::type eval_gcd(gmp_int& result, const gmp_int& a, const I b)
{
   mpz_gcd_ui(result.data(), a.data(), boost::multiprecision::detail::unsigned_abs(b));
}
template <class I>
inline typename enable_if_c<is_signed<I>::value && ((sizeof(I) <= sizeof(long)))>::type eval_lcm(gmp_int& result, const gmp_int& a, const I b)
{
   mpz_lcm_ui(result.data(), a.data(), boost::multiprecision::detail::unsigned_abs(b));
}

inline void eval_integer_sqrt(gmp_int& s, gmp_int& r, const gmp_int& x)
{
   mpz_sqrtrem(s.data(), r.data(), x.data());
}

inline unsigned eval_lsb(const gmp_int& val)
{
   int c = eval_get_sign(val);
   if (c == 0)
   {
      BOOST_THROW_EXCEPTION(std::range_error("No bits were set in the operand."));
   }
   if (c < 0)
   {
      BOOST_THROW_EXCEPTION(std::range_error("Testing individual bits in negative values is not supported - results are undefined."));
   }
   return static_cast<unsigned>(mpz_scan1(val.data(), 0));
}

inline unsigned eval_msb(const gmp_int& val)
{
   int c = eval_get_sign(val);
   if (c == 0)
   {
      BOOST_THROW_EXCEPTION(std::range_error("No bits were set in the operand."));
   }
   if (c < 0)
   {
      BOOST_THROW_EXCEPTION(std::range_error("Testing individual bits in negative values is not supported - results are undefined."));
   }
   return static_cast<unsigned>(mpz_sizeinbase(val.data(), 2) - 1);
}

inline bool eval_bit_test(const gmp_int& val, unsigned index)
{
   return mpz_tstbit(val.data(), index) ? true : false;
}

inline void eval_bit_set(gmp_int& val, unsigned index)
{
   mpz_setbit(val.data(), index);
}

inline void eval_bit_unset(gmp_int& val, unsigned index)
{
   mpz_clrbit(val.data(), index);
}

inline void eval_bit_flip(gmp_int& val, unsigned index)
{
   mpz_combit(val.data(), index);
}

inline void eval_qr(const gmp_int& x, const gmp_int& y,
                    gmp_int& q, gmp_int& r)
{
   mpz_tdiv_qr(q.data(), r.data(), x.data(), y.data());
}

template <class Integer>
inline typename enable_if<is_unsigned<Integer>, Integer>::type eval_integer_modulus(const gmp_int& x, Integer val)
{
#if defined(__MPIR_VERSION) && (__MPIR_VERSION >= 3)
   if ((sizeof(Integer) <= sizeof(mpir_ui)) || (val <= (std::numeric_limits<mpir_ui>::max)()))
#else
   if ((sizeof(Integer) <= sizeof(long)) || (val <= (std::numeric_limits<unsigned long>::max)()))
#endif
   {
      return static_cast<Integer>(mpz_tdiv_ui(x.data(), val));
   }
   else
   {
      return default_ops::eval_integer_modulus(x, val);
   }
}
template <class Integer>
inline typename enable_if<is_signed<Integer>, Integer>::type eval_integer_modulus(const gmp_int& x, Integer val)
{
   return eval_integer_modulus(x, boost::multiprecision::detail::unsigned_abs(val));
}
inline void eval_powm(gmp_int& result, const gmp_int& base, const gmp_int& p, const gmp_int& m)
{
   if (eval_get_sign(p) < 0)
   {
      BOOST_THROW_EXCEPTION(std::runtime_error("powm requires a positive exponent."));
   }
   mpz_powm(result.data(), base.data(), p.data(), m.data());
}

template <class Integer>
inline typename enable_if<
    mpl::and_<
        is_unsigned<Integer>,
        mpl::bool_<sizeof(Integer) <= sizeof(unsigned long)> > >::type
eval_powm(gmp_int& result, const gmp_int& base, Integer p, const gmp_int& m)
{
   mpz_powm_ui(result.data(), base.data(), p, m.data());
}
template <class Integer>
inline typename enable_if<
    mpl::and_<
        is_signed<Integer>,
        mpl::bool_<sizeof(Integer) <= sizeof(unsigned long)> > >::type
eval_powm(gmp_int& result, const gmp_int& base, Integer p, const gmp_int& m)
{
   if (p < 0)
   {
      BOOST_THROW_EXCEPTION(std::runtime_error("powm requires a positive exponent."));
   }
   mpz_powm_ui(result.data(), base.data(), p, m.data());
}

inline std::size_t hash_value(const gmp_int& val)
{
   // We should really use mpz_limbs_read here, but that's unsupported on older versions:
   std::size_t result = 0;
   for (int i = 0; i < std::abs(val.data()[0]._mp_size); ++i)
      boost::hash_combine(result, val.data()[0]._mp_d[i]);
   boost::hash_combine(result, val.data()[0]._mp_size);
   return result;
}

struct gmp_rational;
void eval_add(gmp_rational& t, const gmp_rational& o);

struct gmp_rational
{
#ifdef BOOST_HAS_LONG_LONG
   typedef mpl::list<long, boost::long_long_type>           signed_types;
   typedef mpl::list<unsigned long, boost::ulong_long_type> unsigned_types;
#else
   typedef mpl::list<long>          signed_types;
   typedef mpl::list<unsigned long> unsigned_types;
#endif
   typedef mpl::list<double, long double> float_types;

   gmp_rational()
   {
      mpq_init(this->m_data);
   }
   gmp_rational(const gmp_rational& o)
   {
      mpq_init(m_data);
      if (o.m_data[0]._mp_num._mp_d)
         mpq_set(m_data, o.m_data);
   }
   gmp_rational(const gmp_int& o)
   {
      mpq_init(m_data);
      mpq_set_z(m_data, o.data());
   }
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
   gmp_rational(gmp_rational&& o) BOOST_NOEXCEPT
   {
      m_data[0]                 = o.m_data[0];
      o.m_data[0]._mp_num._mp_d = 0;
      o.m_data[0]._mp_den._mp_d = 0;
   }
#endif
   gmp_rational(const mpq_t o)
   {
      mpq_init(m_data);
      mpq_set(m_data, o);
   }
   gmp_rational(const mpz_t o)
   {
      mpq_init(m_data);
      mpq_set_z(m_data, o);
   }
   gmp_rational& operator=(const gmp_rational& o)
   {
      if (m_data[0]._mp_den._mp_d == 0)
         mpq_init(m_data);
      mpq_set(m_data, o.m_data);
      return *this;
   }
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
   gmp_rational& operator=(gmp_rational&& o) BOOST_NOEXCEPT
   {
      mpq_swap(m_data, o.m_data);
      return *this;
   }
#endif
#ifdef BOOST_HAS_LONG_LONG
#if defined(ULLONG_MAX) && (ULLONG_MAX == ULONG_MAX)
   gmp_rational& operator=(boost::ulong_long_type i)
   {
      *this = static_cast<unsigned long>(i);
      return *this;
   }
#else
   gmp_rational& operator=(boost::ulong_long_type i)
   {
      if (m_data[0]._mp_den._mp_d == 0)
         mpq_init(m_data);
      gmp_int zi;
      zi = i;
      mpq_set_z(m_data, zi.data());
      return *this;
   }
   gmp_rational& operator=(boost::long_long_type i)
   {
      if (m_data[0]._mp_den._mp_d == 0)
         mpq_init(m_data);
      bool neg = i < 0;
      *this    = boost::multiprecision::detail::unsigned_abs(i);
      if (neg)
         mpq_neg(m_data, m_data);
      return *this;
   }
#endif
#endif
   gmp_rational& operator=(unsigned long i)
   {
      if (m_data[0]._mp_den._mp_d == 0)
         mpq_init(m_data);
      mpq_set_ui(m_data, i, 1);
      return *this;
   }
   gmp_rational& operator=(long i)
   {
      if (m_data[0]._mp_den._mp_d == 0)
         mpq_init(m_data);
      mpq_set_si(m_data, i, 1);
      return *this;
   }
   gmp_rational& operator=(double d)
   {
      if (m_data[0]._mp_den._mp_d == 0)
         mpq_init(m_data);
      mpq_set_d(m_data, d);
      return *this;
   }
   gmp_rational& operator=(long double a)
   {
      using default_ops::eval_add;
      using default_ops::eval_subtract;
      using std::floor;
      using std::frexp;
      using std::ldexp;

      if (m_data[0]._mp_den._mp_d == 0)
         mpq_init(m_data);

      if (a == 0)
      {
         mpq_set_si(m_data, 0, 1);
         return *this;
      }

      if (a == 1)
      {
         mpq_set_si(m_data, 1, 1);
         return *this;
      }

      BOOST_ASSERT(!(boost::math::isinf)(a));
      BOOST_ASSERT(!(boost::math::isnan)(a));

      int         e;
      long double f, term;
      mpq_set_ui(m_data, 0, 1);
      mpq_set_ui(m_data, 0u, 1);
      gmp_rational t;

      f = frexp(a, &e);

      static const int shift = std::numeric_limits<int>::digits - 1;

      while (f)
      {
         // extract int sized bits from f:
         f    = ldexp(f, shift);
         term = floor(f);
         e -= shift;
         mpq_mul_2exp(m_data, m_data, shift);
         t = static_cast<long>(term);
         eval_add(*this, t);
         f -= term;
      }
      if (e > 0)
         mpq_mul_2exp(m_data, m_data, e);
      else if (e < 0)
         mpq_div_2exp(m_data, m_data, -e);
      return *this;
   }
   gmp_rational& operator=(const char* s)
   {
      if (m_data[0]._mp_den._mp_d == 0)
         mpq_init(m_data);
      if (0 != mpq_set_str(m_data, s, 10))
         BOOST_THROW_EXCEPTION(std::runtime_error(std::string("The string \"") + s + std::string("\"could not be interpreted as a valid rational number.")));
      return *this;
   }
   gmp_rational& operator=(const gmp_int& o)
   {
      if (m_data[0]._mp_den._mp_d == 0)
         mpq_init(m_data);
      mpq_set_z(m_data, o.data());
      return *this;
   }
   gmp_rational& operator=(const mpq_t o)
   {
      if (m_data[0]._mp_den._mp_d == 0)
         mpq_init(m_data);
      mpq_set(m_data, o);
      return *this;
   }
   gmp_rational& operator=(const mpz_t o)
   {
      if (m_data[0]._mp_den._mp_d == 0)
         mpq_init(m_data);
      mpq_set_z(m_data, o);
      return *this;
   }
   void swap(gmp_rational& o)
   {
      mpq_swap(m_data, o.m_data);
   }
   std::string str(std::streamsize /*digits*/, std::ios_base::fmtflags /*f*/) const
   {
      BOOST_ASSERT(m_data[0]._mp_num._mp_d);
      // TODO make a better job of this including handling of f!!
      void* (*alloc_func_ptr)(size_t);
      void* (*realloc_func_ptr)(void*, size_t, size_t);
      void (*free_func_ptr)(void*, size_t);
      const char* ps = mpq_get_str(0, 10, m_data);
      std::string s  = ps;
      mp_get_memory_functions(&alloc_func_ptr, &realloc_func_ptr, &free_func_ptr);
      (*free_func_ptr)((void*)ps, std::strlen(ps) + 1);
      return s;
   }
   ~gmp_rational()
   {
      if (m_data[0]._mp_num._mp_d || m_data[0]._mp_den._mp_d)
         mpq_clear(m_data);
   }
   void negate()
   {
      BOOST_ASSERT(m_data[0]._mp_num._mp_d);
      mpq_neg(m_data, m_data);
   }
   int compare(const gmp_rational& o) const
   {
      BOOST_ASSERT(m_data[0]._mp_num._mp_d && o.m_data[0]._mp_num._mp_d);
      return mpq_cmp(m_data, o.m_data);
   }
   template <class V>
   int compare(V v) const
   {
      gmp_rational d;
      d = v;
      return compare(d);
   }
   int compare(unsigned long v) const
   {
      BOOST_ASSERT(m_data[0]._mp_num._mp_d);
      return mpq_cmp_ui(m_data, v, 1);
   }
   int compare(long v) const
   {
      BOOST_ASSERT(m_data[0]._mp_num._mp_d);
      return mpq_cmp_si(m_data, v, 1);
   }
   mpq_t& data()
   {
      BOOST_ASSERT(m_data[0]._mp_num._mp_d);
      return m_data;
   }
   const mpq_t& data() const
   {
      BOOST_ASSERT(m_data[0]._mp_num._mp_d);
      return m_data;
   }

 protected:
   mpq_t m_data;
};

inline bool eval_is_zero(const gmp_rational& val)
{
   return mpq_sgn(val.data()) == 0;
}
template <class T>
inline bool eval_eq(gmp_rational& a, const T& b)
{
   return a.compare(b) == 0;
}
template <class T>
inline bool eval_lt(gmp_rational& a, const T& b)
{
   return a.compare(b) < 0;
}
template <class T>
inline bool eval_gt(gmp_rational& a, const T& b)
{
   return a.compare(b) > 0;
}

inline void eval_add(gmp_rational& t, const gmp_rational& o)
{
   mpq_add(t.data(), t.data(), o.data());
}
inline void eval_subtract(gmp_rational& t, const gmp_rational& o)
{
   mpq_sub(t.data(), t.data(), o.data());
}
inline void eval_multiply(gmp_rational& t, const gmp_rational& o)
{
   mpq_mul(t.data(), t.data(), o.data());
}
inline void eval_divide(gmp_rational& t, const gmp_rational& o)
{
   if (eval_is_zero(o))
      BOOST_THROW_EXCEPTION(std::overflow_error("Division by zero."));
   mpq_div(t.data(), t.data(), o.data());
}
inline void eval_add(gmp_rational& t, const gmp_rational& p, const gmp_rational& o)
{
   mpq_add(t.data(), p.data(), o.data());
}
inline void eval_subtract(gmp_rational& t, const gmp_rational& p, const gmp_rational& o)
{
   mpq_sub(t.data(), p.data(), o.data());
}
inline void eval_multiply(gmp_rational& t, const gmp_rational& p, const gmp_rational& o)
{
   mpq_mul(t.data(), p.data(), o.data());
}
inline void eval_divide(gmp_rational& t, const gmp_rational& p, const gmp_rational& o)
{
   if (eval_is_zero(o))
      BOOST_THROW_EXCEPTION(std::overflow_error("Division by zero."));
   mpq_div(t.data(), p.data(), o.data());
}

inline int eval_get_sign(const gmp_rational& val)
{
   return mpq_sgn(val.data());
}
template <class R>
inline typename enable_if_c<number_category<R>::value == number_kind_floating_point>::type eval_convert_to(R* result, const gmp_rational& backend)
{
   //
   // The generic conversion is as good as anything we can write here:
   //
   // This does not round correctly:
   //
   //*result = mpq_get_d(val.data());
   //
   // This does:
   //
   ::boost::multiprecision::detail::generic_convert_rational_to_float(*result, backend);
}
#ifdef BOOST_HAS_FLOAT128
inline void eval_convert_to(__float128* result, const gmp_rational& val)
{
   using default_ops::eval_convert_to;

   gmp_int n, d;
   __float128 fn, fd;
   mpz_set(n.data(), mpq_numref(val.data()));
   mpz_set(d.data(), mpq_denref(val.data()));

   eval_convert_to(&fn, n);
   eval_convert_to(&fd, d);
   
   *result = fn / fd;
}
#endif

inline void eval_convert_to(long* result, const gmp_rational& val)
{
   double r;
   eval_convert_to(&r, val);
   *result = static_cast<long>(r);
}

inline void eval_convert_to(unsigned long* result, const gmp_rational& val)
{
   double r;
   eval_convert_to(&r, val);
   *result = static_cast<long>(r);
}

inline void eval_abs(gmp_rational& result, const gmp_rational& val)
{
   mpq_abs(result.data(), val.data());
}

inline void assign_components(gmp_rational& result, unsigned long v1, unsigned long v2)
{
   mpq_set_ui(result.data(), v1, v2);
   mpq_canonicalize(result.data());
}
inline void assign_components(gmp_rational& result, long v1, long v2)
{
   mpq_set_si(result.data(), v1, v2);
   mpq_canonicalize(result.data());
}
inline void assign_components(gmp_rational& result, gmp_int const& v1, gmp_int const& v2)
{
   mpz_set(mpq_numref(result.data()), v1.data());
   mpz_set(mpq_denref(result.data()), v2.data());
   mpq_canonicalize(result.data());
}

inline std::size_t hash_value(const gmp_rational& val)
{
   std::size_t result = 0;
   for (int i = 0; i < std::abs(val.data()[0]._mp_num._mp_size); ++i)
      boost::hash_combine(result, val.data()[0]._mp_num._mp_d[i]);
   for (int i = 0; i < std::abs(val.data()[0]._mp_den._mp_size); ++i)
      boost::hash_combine(result, val.data()[0]._mp_den._mp_d[i]);
   boost::hash_combine(result, val.data()[0]._mp_num._mp_size);
   return result;
}

//
// Some member functions that are dependent upon previous code go here:
//
template <unsigned Digits10>
template <unsigned D>
inline gmp_float<Digits10>::gmp_float(const gmp_float<D>& o, typename enable_if_c<D <= Digits10>::type*)
{
   mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(Digits10 ? Digits10 : (unsigned)this->get_default_precision()));
   mpf_set(this->m_data, o.data());
}
template <unsigned Digits10>
template <unsigned D>
inline gmp_float<Digits10>::gmp_float(const gmp_float<D>& o, typename disable_if_c<D <= Digits10>::type*)
{
   mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(Digits10 ? Digits10 : (unsigned)this->get_default_precision()));
   mpf_set(this->m_data, o.data());
}
template <unsigned Digits10>
inline gmp_float<Digits10>::gmp_float(const gmp_int& o)
{
   mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(Digits10 ? Digits10 : (unsigned)this->get_default_precision()));
   mpf_set_z(this->data(), o.data());
}
template <unsigned Digits10>
inline gmp_float<Digits10>::gmp_float(const gmp_rational& o)
{
   mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(Digits10 ? Digits10 : (unsigned)this->get_default_precision()));
   mpf_set_q(this->data(), o.data());
}
template <unsigned Digits10>
template <unsigned D>
inline gmp_float<Digits10>& gmp_float<Digits10>::operator=(const gmp_float<D>& o)
{
   if (this->m_data[0]._mp_d == 0)
      mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(Digits10 ? Digits10 : (unsigned)this->get_default_precision()));
   mpf_set(this->m_data, o.data());
   return *this;
}
template <unsigned Digits10>
inline gmp_float<Digits10>& gmp_float<Digits10>::operator=(const gmp_int& o)
{
   if (this->m_data[0]._mp_d == 0)
      mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(Digits10 ? Digits10 : (unsigned)this->get_default_precision()));
   mpf_set_z(this->data(), o.data());
   return *this;
}
template <unsigned Digits10>
inline gmp_float<Digits10>& gmp_float<Digits10>::operator=(const gmp_rational& o)
{
   if (this->m_data[0]._mp_d == 0)
      mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(Digits10 ? Digits10 : (unsigned)this->get_default_precision()));
   mpf_set_q(this->data(), o.data());
   return *this;
}
inline gmp_float<0>::gmp_float(const gmp_int& o) : requested_precision(get_default_precision())
{
   mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(requested_precision));
   mpf_set_z(this->data(), o.data());
}
inline gmp_float<0>::gmp_float(const gmp_rational& o) : requested_precision(get_default_precision())
{
   mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(requested_precision));
   mpf_set_q(this->data(), o.data());
}
inline gmp_float<0>& gmp_float<0>::operator=(const gmp_int& o)
{
   if (this->m_data[0]._mp_d == 0)
   {
      requested_precision = this->get_default_precision();
      mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(requested_precision));
   }
   mpf_set_z(this->data(), o.data());
   return *this;
}
inline gmp_float<0>& gmp_float<0>::operator=(const gmp_rational& o)
{
   if (this->m_data[0]._mp_d == 0)
   {
      requested_precision = this->get_default_precision();
      mpf_init2(this->m_data, multiprecision::detail::digits10_2_2(requested_precision));
   }
   mpf_set_q(this->data(), o.data());
   return *this;
}
inline gmp_int::gmp_int(const gmp_rational& o)
{
   mpz_init(this->m_data);
   mpz_set_q(this->m_data, o.data());
}
inline gmp_int& gmp_int::operator=(const gmp_rational& o)
{
   if (this->m_data[0]._mp_d == 0)
      mpz_init(this->m_data);
   mpz_set_q(this->m_data, o.data());
   return *this;
}

} //namespace backends

using boost::multiprecision::backends::gmp_float;
using boost::multiprecision::backends::gmp_int;
using boost::multiprecision::backends::gmp_rational;

template <expression_template_option ExpressionTemplates>
struct component_type<number<gmp_rational, ExpressionTemplates> >
{
   typedef number<gmp_int, ExpressionTemplates> type;
};

template <expression_template_option ET>
inline number<gmp_int, ET> numerator(const number<gmp_rational, ET>& val)
{
   number<gmp_int, ET> result;
   mpz_set(result.backend().data(), (mpq_numref(val.backend().data())));
   return result;
}
template <expression_template_option ET>
inline number<gmp_int, ET> denominator(const number<gmp_rational, ET>& val)
{
   number<gmp_int, ET> result;
   mpz_set(result.backend().data(), (mpq_denref(val.backend().data())));
   return result;
}

namespace detail {

#ifdef BOOST_NO_SFINAE_EXPR

template <>
struct is_explicitly_convertible<canonical<mpf_t, gmp_int>::type, gmp_int> : public mpl::true_
{};
template <>
struct is_explicitly_convertible<canonical<mpq_t, gmp_int>::type, gmp_int> : public mpl::true_
{};
template <unsigned Digits10>
struct is_explicitly_convertible<gmp_float<Digits10>, gmp_int> : public mpl::true_
{};
template <>
struct is_explicitly_convertible<gmp_rational, gmp_int> : public mpl::true_
{};
template <unsigned D1, unsigned D2>
struct is_explicitly_convertible<gmp_float<D1>, gmp_float<D2> > : public mpl::true_
{};

#endif

template <>
struct digits2<number<gmp_float<0>, et_on> >
{
   static long value()
   {
      return multiprecision::detail::digits10_2_2(gmp_float<0>::default_precision());
   }
};

template <>
struct digits2<number<gmp_float<0>, et_off> >
{
   static long value()
   {
      return multiprecision::detail::digits10_2_2(gmp_float<0>::default_precision());
   }
};

template <>
struct digits2<number<debug_adaptor<gmp_float<0> >, et_on> >
{
   static long value()
   {
      return multiprecision::detail::digits10_2_2(gmp_float<0>::default_precision());
   }
};

template <>
struct digits2<number<debug_adaptor<gmp_float<0> >, et_off> >
{
   static long value()
   {
      return multiprecision::detail::digits10_2_2(gmp_float<0>::default_precision());
   }
};

} // namespace detail

template <>
struct number_category<detail::canonical<mpz_t, gmp_int>::type> : public mpl::int_<number_kind_integer>
{};
template <>
struct number_category<detail::canonical<mpq_t, gmp_rational>::type> : public mpl::int_<number_kind_rational>
{};
template <>
struct number_category<detail::canonical<mpf_t, gmp_float<0> >::type> : public mpl::int_<number_kind_floating_point>
{};

namespace detail {
template <>
struct is_variable_precision<backends::gmp_float<0> > : public true_type
{};
} // namespace detail

typedef number<gmp_float<50> >   mpf_float_50;
typedef number<gmp_float<100> >  mpf_float_100;
typedef number<gmp_float<500> >  mpf_float_500;
typedef number<gmp_float<1000> > mpf_float_1000;
typedef number<gmp_float<0> >    mpf_float;
typedef number<gmp_int>          mpz_int;
typedef number<gmp_rational>     mpq_rational;

} // namespace multiprecision

namespace math { namespace tools {

inline void set_output_precision(const boost::multiprecision::mpf_float& val, std::ostream& os)
{
   os << std::setprecision(val.precision());
}

template <>
inline int digits<boost::multiprecision::mpf_float>()
#ifdef BOOST_MATH_NOEXCEPT
    BOOST_NOEXCEPT
#endif
{
   return multiprecision::detail::digits10_2_2(boost::multiprecision::mpf_float::default_precision());
}
template <>
inline int digits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, boost::multiprecision::et_off> >()
#ifdef BOOST_MATH_NOEXCEPT
    BOOST_NOEXCEPT
#endif
{
   return multiprecision::detail::digits10_2_2(boost::multiprecision::mpf_float::default_precision());
}

template <>
inline boost::multiprecision::mpf_float
max_value<boost::multiprecision::mpf_float>()
{
   boost::multiprecision::mpf_float result(0.5);
   mpf_mul_2exp(result.backend().data(), result.backend().data(), (std::numeric_limits<mp_exp_t>::max)() / 64 + 1);
   return result;
}

template <>
inline boost::multiprecision::mpf_float
min_value<boost::multiprecision::mpf_float>()
{
   boost::multiprecision::mpf_float result(0.5);
   mpf_div_2exp(result.backend().data(), result.backend().data(), (std::numeric_limits<mp_exp_t>::max)() / 64 + 1);
   return result;
}

template <>
inline boost::multiprecision::number<boost::multiprecision::gmp_float<0>, boost::multiprecision::et_off>
max_value<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, boost::multiprecision::et_off> >()
{
   boost::multiprecision::number<boost::multiprecision::gmp_float<0>, boost::multiprecision::et_off> result(0.5);
   mpf_mul_2exp(result.backend().data(), result.backend().data(), (std::numeric_limits<mp_exp_t>::max)() / 64 + 1);
   return result;
}

template <>
inline boost::multiprecision::number<boost::multiprecision::gmp_float<0>, boost::multiprecision::et_off>
min_value<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, boost::multiprecision::et_off> >()
{
   boost::multiprecision::number<boost::multiprecision::gmp_float<0>, boost::multiprecision::et_off> result(0.5);
   mpf_div_2exp(result.backend().data(), result.backend().data(), (std::numeric_limits<mp_exp_t>::max)() / 64 + 1);
   return result;
}

template <>
inline int digits<boost::multiprecision::number<boost::multiprecision::debug_adaptor<boost::multiprecision::mpf_float::backend_type> > >()
#ifdef BOOST_MATH_NOEXCEPT
    BOOST_NOEXCEPT
#endif
{
   return multiprecision::detail::digits10_2_2(boost::multiprecision::number<boost::multiprecision::debug_adaptor<boost::multiprecision::mpf_float::backend_type> >::default_precision());
}
template <>
inline int digits<boost::multiprecision::number<boost::multiprecision::debug_adaptor<boost::multiprecision::gmp_float<0> >, boost::multiprecision::et_off> >()
#ifdef BOOST_MATH_NOEXCEPT
    BOOST_NOEXCEPT
#endif
{
   return multiprecision::detail::digits10_2_2(boost::multiprecision::number<boost::multiprecision::debug_adaptor<boost::multiprecision::mpf_float::backend_type> >::default_precision());
}

template <>
inline boost::multiprecision::number<boost::multiprecision::debug_adaptor<boost::multiprecision::mpf_float::backend_type> >
max_value<boost::multiprecision::number<boost::multiprecision::debug_adaptor<boost::multiprecision::mpf_float::backend_type> > >()
{
   return max_value<boost::multiprecision::mpf_float>().backend();
}

template <>
inline boost::multiprecision::number<boost::multiprecision::debug_adaptor<boost::multiprecision::mpf_float::backend_type> >
min_value<boost::multiprecision::number<boost::multiprecision::debug_adaptor<boost::multiprecision::mpf_float::backend_type> > >()
{
   return min_value<boost::multiprecision::mpf_float>().backend();
}

template <>
inline boost::multiprecision::number<boost::multiprecision::debug_adaptor<boost::multiprecision::gmp_float<0> >, boost::multiprecision::et_off>
max_value<boost::multiprecision::number<boost::multiprecision::debug_adaptor<boost::multiprecision::gmp_float<0> >, boost::multiprecision::et_off> >()
{
   return max_value<boost::multiprecision::mpf_float>().backend();
}

template <>
inline boost::multiprecision::number<boost::multiprecision::debug_adaptor<boost::multiprecision::gmp_float<0> >, boost::multiprecision::et_off>
min_value<boost::multiprecision::number<boost::multiprecision::debug_adaptor<boost::multiprecision::gmp_float<0> >, boost::multiprecision::et_off> >()
{
   return min_value<boost::multiprecision::mpf_float>().backend();
}

}} // namespace math::tools

} // namespace boost

namespace std {

//
// numeric_limits [partial] specializations for the types declared in this header:
//
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
class numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >
{
   typedef boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> number_type;

 public:
   BOOST_STATIC_CONSTEXPR bool is_specialized = true;
   //
   // min and max values chosen so as to not cause segfaults when calling
   // mpf_get_str on 64-bit Linux builds.  Possibly we could use larger
   // exponent values elsewhere.
   //
   static number_type(min)()
   {
      initializer.do_nothing();
      static std::pair<bool, number_type> value;
      if (!value.first)
      {
         value.first  = true;
         value.second = 1;
         mpf_div_2exp(value.second.backend().data(), value.second.backend().data(), (std::numeric_limits<mp_exp_t>::max)() / 64 + 1);
      }
      return value.second;
   }
   static number_type(max)()
   {
      initializer.do_nothing();
      static std::pair<bool, number_type> value;
      if (!value.first)
      {
         value.first  = true;
         value.second = 1;
         mpf_mul_2exp(value.second.backend().data(), value.second.backend().data(), (std::numeric_limits<mp_exp_t>::max)() / 64 + 1);
      }
      return value.second;
   }
   BOOST_STATIC_CONSTEXPR number_type lowest()
   {
      return -(max)();
   }
   BOOST_STATIC_CONSTEXPR int digits   = static_cast<int>((Digits10 * 1000L) / 301L + ((Digits10 * 1000L) % 301L ? 2 : 1));
   BOOST_STATIC_CONSTEXPR int digits10 = Digits10;
   // Have to allow for a possible extra limb inside the gmp data structure:
   BOOST_STATIC_CONSTEXPR int  max_digits10 = Digits10 + 3 + ((GMP_LIMB_BITS * 301L) / 1000L);
   BOOST_STATIC_CONSTEXPR bool is_signed    = true;
   BOOST_STATIC_CONSTEXPR bool is_integer   = false;
   BOOST_STATIC_CONSTEXPR bool is_exact     = false;
   BOOST_STATIC_CONSTEXPR int  radix        = 2;
   static number_type          epsilon()
   {
      initializer.do_nothing();
      static std::pair<bool, number_type> value;
      if (!value.first)
      {
         value.first  = true;
         value.second = 1;
         mpf_div_2exp(value.second.backend().data(), value.second.backend().data(), std::numeric_limits<number_type>::digits - 1);
      }
      return value.second;
   }
   // What value should this be????
   static number_type round_error()
   {
      // returns epsilon/2
      initializer.do_nothing();
      static std::pair<bool, number_type> value;
      if (!value.first)
      {
         value.first  = true;
         value.second = 1;
      }
      return value.second;
   }
   BOOST_STATIC_CONSTEXPR long               min_exponent      = LONG_MIN;
   BOOST_STATIC_CONSTEXPR long               min_exponent10    = (LONG_MIN / 1000) * 301L;
   BOOST_STATIC_CONSTEXPR long               max_exponent      = LONG_MAX;
   BOOST_STATIC_CONSTEXPR long               max_exponent10    = (LONG_MAX / 1000) * 301L;
   BOOST_STATIC_CONSTEXPR bool               has_infinity      = false;
   BOOST_STATIC_CONSTEXPR bool               has_quiet_NaN     = false;
   BOOST_STATIC_CONSTEXPR bool               has_signaling_NaN = false;
   BOOST_STATIC_CONSTEXPR float_denorm_style has_denorm        = denorm_absent;
   BOOST_STATIC_CONSTEXPR bool               has_denorm_loss   = false;
   BOOST_STATIC_CONSTEXPR number_type        infinity() { return number_type(); }
   BOOST_STATIC_CONSTEXPR number_type        quiet_NaN() { return number_type(); }
   BOOST_STATIC_CONSTEXPR number_type        signaling_NaN() { return number_type(); }
   BOOST_STATIC_CONSTEXPR number_type        denorm_min() { return number_type(); }
   BOOST_STATIC_CONSTEXPR bool               is_iec559       = false;
   BOOST_STATIC_CONSTEXPR bool               is_bounded      = true;
   BOOST_STATIC_CONSTEXPR bool               is_modulo       = false;
   BOOST_STATIC_CONSTEXPR bool               traps           = true;
   BOOST_STATIC_CONSTEXPR bool               tinyness_before = false;
   BOOST_STATIC_CONSTEXPR float_round_style  round_style     = round_indeterminate;

 private:
   struct data_initializer
   {
      data_initializer()
      {
         std::numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<digits10> > >::epsilon();
         std::numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<digits10> > >::round_error();
         (std::numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<digits10> > >::min)();
         (std::numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<digits10> > >::max)();
      }
      void do_nothing() const {}
   };
   static const data_initializer initializer;
};

template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
const typename numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::data_initializer numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::initializer;

#ifndef BOOST_NO_INCLASS_MEMBER_INITIALIZATION

template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::digits;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::digits10;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::max_digits10;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::is_signed;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::is_integer;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::is_exact;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::radix;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST long numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::min_exponent;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST long numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::min_exponent10;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST long numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::max_exponent;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST long numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::max_exponent10;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::has_infinity;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::has_quiet_NaN;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::has_signaling_NaN;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST float_denorm_style numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::has_denorm;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::has_denorm_loss;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::is_iec559;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::is_bounded;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::is_modulo;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::traps;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::tinyness_before;
template <unsigned Digits10, boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST float_round_style numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<Digits10>, ExpressionTemplates> >::round_style;

#endif

template <boost::multiprecision::expression_template_option ExpressionTemplates>
class numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >
{
   typedef boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> number_type;

 public:
   BOOST_STATIC_CONSTEXPR bool is_specialized = false;
   static number_type(min)() { return number_type(); }
   static number_type(max)() { return number_type(); }
   static number_type                        lowest() { return number_type(); }
   BOOST_STATIC_CONSTEXPR int                digits       = 0;
   BOOST_STATIC_CONSTEXPR int                digits10     = 0;
   BOOST_STATIC_CONSTEXPR int                max_digits10 = 0;
   BOOST_STATIC_CONSTEXPR bool               is_signed    = false;
   BOOST_STATIC_CONSTEXPR bool               is_integer   = false;
   BOOST_STATIC_CONSTEXPR bool               is_exact     = false;
   BOOST_STATIC_CONSTEXPR int                radix        = 0;
   static number_type                        epsilon() { return number_type(); }
   static number_type                        round_error() { return number_type(); }
   BOOST_STATIC_CONSTEXPR int                min_exponent      = 0;
   BOOST_STATIC_CONSTEXPR int                min_exponent10    = 0;
   BOOST_STATIC_CONSTEXPR int                max_exponent      = 0;
   BOOST_STATIC_CONSTEXPR int                max_exponent10    = 0;
   BOOST_STATIC_CONSTEXPR bool               has_infinity      = false;
   BOOST_STATIC_CONSTEXPR bool               has_quiet_NaN     = false;
   BOOST_STATIC_CONSTEXPR bool               has_signaling_NaN = false;
   BOOST_STATIC_CONSTEXPR float_denorm_style has_denorm        = denorm_absent;
   BOOST_STATIC_CONSTEXPR bool               has_denorm_loss   = false;
   static number_type                        infinity() { return number_type(); }
   static number_type                        quiet_NaN() { return number_type(); }
   static number_type                        signaling_NaN() { return number_type(); }
   static number_type                        denorm_min() { return number_type(); }
   BOOST_STATIC_CONSTEXPR bool               is_iec559       = false;
   BOOST_STATIC_CONSTEXPR bool               is_bounded      = false;
   BOOST_STATIC_CONSTEXPR bool               is_modulo       = false;
   BOOST_STATIC_CONSTEXPR bool               traps           = false;
   BOOST_STATIC_CONSTEXPR bool               tinyness_before = false;
   BOOST_STATIC_CONSTEXPR float_round_style  round_style     = round_indeterminate;
};

#ifndef BOOST_NO_INCLASS_MEMBER_INITIALIZATION

template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::digits;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::digits10;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::max_digits10;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::is_signed;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::is_integer;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::is_exact;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::radix;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::min_exponent;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::min_exponent10;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::max_exponent;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::max_exponent10;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::has_infinity;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::has_quiet_NaN;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::has_signaling_NaN;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST float_denorm_style numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::has_denorm;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::has_denorm_loss;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::is_iec559;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::is_bounded;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::is_modulo;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::traps;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::tinyness_before;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST float_round_style numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_float<0>, ExpressionTemplates> >::round_style;

#endif

template <boost::multiprecision::expression_template_option ExpressionTemplates>
class numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >
{
   typedef boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> number_type;

 public:
   BOOST_STATIC_CONSTEXPR bool is_specialized = true;
   //
   // Largest and smallest numbers are bounded only by available memory, set
   // to zero:
   //
   static number_type(min)()
   {
      return number_type();
   }
   static number_type(max)()
   {
      return number_type();
   }
   static number_type                        lowest() { return (min)(); }
   BOOST_STATIC_CONSTEXPR int                digits       = INT_MAX;
   BOOST_STATIC_CONSTEXPR int                digits10     = (INT_MAX / 1000) * 301L;
   BOOST_STATIC_CONSTEXPR int                max_digits10 = digits10 + 3;
   BOOST_STATIC_CONSTEXPR bool               is_signed    = true;
   BOOST_STATIC_CONSTEXPR bool               is_integer   = true;
   BOOST_STATIC_CONSTEXPR bool               is_exact     = true;
   BOOST_STATIC_CONSTEXPR int                radix        = 2;
   static number_type                        epsilon() { return number_type(); }
   static number_type                        round_error() { return number_type(); }
   BOOST_STATIC_CONSTEXPR int                min_exponent      = 0;
   BOOST_STATIC_CONSTEXPR int                min_exponent10    = 0;
   BOOST_STATIC_CONSTEXPR int                max_exponent      = 0;
   BOOST_STATIC_CONSTEXPR int                max_exponent10    = 0;
   BOOST_STATIC_CONSTEXPR bool               has_infinity      = false;
   BOOST_STATIC_CONSTEXPR bool               has_quiet_NaN     = false;
   BOOST_STATIC_CONSTEXPR bool               has_signaling_NaN = false;
   BOOST_STATIC_CONSTEXPR float_denorm_style has_denorm        = denorm_absent;
   BOOST_STATIC_CONSTEXPR bool               has_denorm_loss   = false;
   static number_type                        infinity() { return number_type(); }
   static number_type                        quiet_NaN() { return number_type(); }
   static number_type                        signaling_NaN() { return number_type(); }
   static number_type                        denorm_min() { return number_type(); }
   BOOST_STATIC_CONSTEXPR bool               is_iec559       = false;
   BOOST_STATIC_CONSTEXPR bool               is_bounded      = false;
   BOOST_STATIC_CONSTEXPR bool               is_modulo       = false;
   BOOST_STATIC_CONSTEXPR bool               traps           = false;
   BOOST_STATIC_CONSTEXPR bool               tinyness_before = false;
   BOOST_STATIC_CONSTEXPR float_round_style  round_style     = round_toward_zero;
};

#ifndef BOOST_NO_INCLASS_MEMBER_INITIALIZATION

template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::digits;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::digits10;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::max_digits10;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::is_signed;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::is_integer;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::is_exact;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::radix;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::min_exponent;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::min_exponent10;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::max_exponent;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::max_exponent10;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::has_infinity;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::has_quiet_NaN;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::has_signaling_NaN;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST float_denorm_style numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::has_denorm;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::has_denorm_loss;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::is_iec559;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::is_bounded;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::is_modulo;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::traps;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::tinyness_before;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST float_round_style numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_int, ExpressionTemplates> >::round_style;

#endif

template <boost::multiprecision::expression_template_option ExpressionTemplates>
class numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >
{
   typedef boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> number_type;

 public:
   BOOST_STATIC_CONSTEXPR bool is_specialized = true;
   //
   // Largest and smallest numbers are bounded only by available memory, set
   // to zero:
   //
   static number_type(min)()
   {
      return number_type();
   }
   static number_type(max)()
   {
      return number_type();
   }
   static number_type lowest() { return (min)(); }
   // Digits are unbounded, use zero for now:
   BOOST_STATIC_CONSTEXPR int                digits       = INT_MAX;
   BOOST_STATIC_CONSTEXPR int                digits10     = (INT_MAX / 1000) * 301L;
   BOOST_STATIC_CONSTEXPR int                max_digits10 = digits10 + 3;
   BOOST_STATIC_CONSTEXPR bool               is_signed    = true;
   BOOST_STATIC_CONSTEXPR bool               is_integer   = false;
   BOOST_STATIC_CONSTEXPR bool               is_exact     = true;
   BOOST_STATIC_CONSTEXPR int                radix        = 2;
   static number_type                        epsilon() { return number_type(); }
   static number_type                        round_error() { return number_type(); }
   BOOST_STATIC_CONSTEXPR int                min_exponent      = 0;
   BOOST_STATIC_CONSTEXPR int                min_exponent10    = 0;
   BOOST_STATIC_CONSTEXPR int                max_exponent      = 0;
   BOOST_STATIC_CONSTEXPR int                max_exponent10    = 0;
   BOOST_STATIC_CONSTEXPR bool               has_infinity      = false;
   BOOST_STATIC_CONSTEXPR bool               has_quiet_NaN     = false;
   BOOST_STATIC_CONSTEXPR bool               has_signaling_NaN = false;
   BOOST_STATIC_CONSTEXPR float_denorm_style has_denorm        = denorm_absent;
   BOOST_STATIC_CONSTEXPR bool               has_denorm_loss   = false;
   static number_type                        infinity() { return number_type(); }
   static number_type                        quiet_NaN() { return number_type(); }
   static number_type                        signaling_NaN() { return number_type(); }
   static number_type                        denorm_min() { return number_type(); }
   BOOST_STATIC_CONSTEXPR bool               is_iec559       = false;
   BOOST_STATIC_CONSTEXPR bool               is_bounded      = false;
   BOOST_STATIC_CONSTEXPR bool               is_modulo       = false;
   BOOST_STATIC_CONSTEXPR bool               traps           = false;
   BOOST_STATIC_CONSTEXPR bool               tinyness_before = false;
   BOOST_STATIC_CONSTEXPR float_round_style  round_style     = round_toward_zero;
};

#ifndef BOOST_NO_INCLASS_MEMBER_INITIALIZATION

template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::digits;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::digits10;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::max_digits10;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::is_signed;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::is_integer;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::is_exact;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::radix;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::min_exponent;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::min_exponent10;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::max_exponent;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST int numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::max_exponent10;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::has_infinity;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::has_quiet_NaN;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::has_signaling_NaN;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST float_denorm_style numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::has_denorm;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::has_denorm_loss;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::is_iec559;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::is_bounded;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::is_modulo;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::traps;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST bool numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::tinyness_before;
template <boost::multiprecision::expression_template_option ExpressionTemplates>
BOOST_CONSTEXPR_OR_CONST float_round_style numeric_limits<boost::multiprecision::number<boost::multiprecision::gmp_rational, ExpressionTemplates> >::round_style;

#endif

#ifdef BOOST_MSVC
#pragma warning(pop)
#endif

} // namespace std

#endif
