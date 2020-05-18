///////////////////////////////////////////////////////////////
//  Copyright 2020 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MP_ADD_UNSIGNED_ADDC_32_HPP
#define BOOST_MP_ADD_UNSIGNED_ADDC_32_HPP

#include <boost/multiprecision/cpp_int/intel_intrinsics.hpp>

namespace boost { namespace multiprecision { namespace backends {

template <class CppInt1, class CppInt2, class CppInt3>
inline BOOST_MP_CXX14_CONSTEXPR void add_unsigned_constexpr(CppInt1& result, const CppInt2& a, const CppInt3& b) BOOST_MP_NOEXCEPT_IF(is_non_throwing_cpp_int<CppInt1>::value)
{
   using ::boost::multiprecision::std_constexpr::swap;

   // Nothing fancy, just let uintmax_t take the strain:
   double_limb_type carry = 0;
   unsigned         m(0), x(0);
   unsigned         as = a.size();
   unsigned         bs = b.size();
   minmax(as, bs, m, x);
   if (x == 1)
   {
      bool s = a.sign();
      result = static_cast<double_limb_type>(*a.limbs()) + static_cast<double_limb_type>(*b.limbs());
      result.sign(s);
      return;
   }
   result.resize(x, x);
   typename CppInt2::const_limb_pointer pa     = a.limbs();
   typename CppInt3::const_limb_pointer pb     = b.limbs();
   typename CppInt1::limb_pointer       pr     = result.limbs();
   typename CppInt1::limb_pointer       pr_end = pr + m;

   if (as < bs)
      swap(pa, pb);

   // First where a and b overlap:
   while (pr != pr_end)
   {
      carry += static_cast<double_limb_type>(*pa) + static_cast<double_limb_type>(*pb);
#ifdef __MSVC_RUNTIME_CHECKS
      *pr = static_cast<limb_type>(carry & ~static_cast<limb_type>(0));
#else
      *pr = static_cast<limb_type>(carry);
#endif
      carry >>= CppInt1::limb_bits;
      ++pr, ++pa, ++pb;
   }
   pr_end += x - m;
   // Now where only a has digits:
   while (pr != pr_end)
   {
      if (!carry)
      {
         if (pa != pr)
            std_constexpr::copy(pa, pa + (pr_end - pr), pr);
         break;
      }
      carry += static_cast<double_limb_type>(*pa);
#ifdef __MSVC_RUNTIME_CHECKS
      *pr = static_cast<limb_type>(carry & ~static_cast<limb_type>(0));
#else
      *pr = static_cast<limb_type>(carry);
#endif
      carry >>= CppInt1::limb_bits;
      ++pr, ++pa;
   }
   if (carry)
   {
      // We overflowed, need to add one more limb:
      result.resize(x + 1, x + 1);
      if (result.size() > x)
         result.limbs()[x] = static_cast<limb_type>(carry);
   }
   result.normalize();
   result.sign(a.sign());
}
//
// Core subtraction routine for all non-trivial cpp_int's:
//
template <class CppInt1, class CppInt2, class CppInt3>
inline BOOST_MP_CXX14_CONSTEXPR void subtract_unsigned_constexpr(CppInt1& result, const CppInt2& a, const CppInt3& b) BOOST_MP_NOEXCEPT_IF(is_non_throwing_cpp_int<CppInt1>::value)
{
   using ::boost::multiprecision::std_constexpr::swap;

   // Nothing fancy, just let uintmax_t take the strain:
   double_limb_type borrow = 0;
   unsigned         m(0), x(0);
   minmax(a.size(), b.size(), m, x);
   //
   // special cases for small limb counts:
   //
   if (x == 1)
   {
      bool      s  = a.sign();
      limb_type al = *a.limbs();
      limb_type bl = *b.limbs();
      if (bl > al)
      {
         ::boost::multiprecision::std_constexpr::swap(al, bl);
         s = !s;
      }
      result = al - bl;
      result.sign(s);
      return;
   }
   // This isn't used till later, but comparison has to occur before we resize the result,
   // as that may also resize a or b if this is an inplace operation:
   int c = a.compare_unsigned(b);
   // Set up the result vector:
   result.resize(x, x);
   // Now that a, b, and result are stable, get pointers to their limbs:
   typename CppInt2::const_limb_pointer pa      = a.limbs();
   typename CppInt3::const_limb_pointer pb      = b.limbs();
   typename CppInt1::limb_pointer       pr      = result.limbs();
   bool                                 swapped = false;
   if (c < 0)
   {
      swap(pa, pb);
      swapped = true;
   }
   else if (c == 0)
   {
      result = static_cast<limb_type>(0);
      return;
   }

   unsigned i = 0;
   // First where a and b overlap:
   while (i < m)
   {
      borrow = static_cast<double_limb_type>(pa[i]) - static_cast<double_limb_type>(pb[i]) - borrow;
      pr[i]  = static_cast<limb_type>(borrow);
      borrow = (borrow >> CppInt1::limb_bits) & 1u;
      ++i;
   }
   // Now where only a has digits, only as long as we've borrowed:
   while (borrow && (i < x))
   {
      borrow = static_cast<double_limb_type>(pa[i]) - borrow;
      pr[i]  = static_cast<limb_type>(borrow);
      borrow = (borrow >> CppInt1::limb_bits) & 1u;
      ++i;
   }
   // Any remaining digits are the same as those in pa:
   if ((x != i) && (pa != pr))
      std_constexpr::copy(pa + i, pa + x, pr + i);
   BOOST_ASSERT(0 == borrow);

   //
   // We may have lost digits, if so update limb usage count:
   //
   result.normalize();
   result.sign(a.sign());
   if (swapped)
      result.negate();
}


#ifdef BOOST_MP_HAS_IMMINTRIN_H
//
// This is the key addition routine where all the argument types are non-trivial cpp_int's:
//
template <class CppInt1, class CppInt2, class CppInt3>
inline BOOST_MP_CXX14_CONSTEXPR void add_unsigned(CppInt1& result, const CppInt2& a, const CppInt3& b) BOOST_MP_NOEXCEPT_IF(is_non_throwing_cpp_int<CppInt1>::value)
{
#ifndef BOOST_MP_NO_CONSTEXPR_DETECTION
   if (BOOST_MP_IS_CONST_EVALUATED(a.size()))
   {
      add_unsigned_constexpr(result, a, b);
   }
   else
#endif
   {
      using std::swap;

      // Nothing fancy, just let uintmax_t take the strain:
      unsigned m(0), x(0);
      unsigned as = a.size();
      unsigned bs = b.size();
      minmax(as, bs, m, x);
      if (x == 1)
      {
         bool s = a.sign();
         result = static_cast<double_limb_type>(*a.limbs()) + static_cast<double_limb_type>(*b.limbs());
         result.sign(s);
         return;
      }
      result.resize(x, x);
      typename CppInt2::const_limb_pointer pa = a.limbs();
      typename CppInt3::const_limb_pointer pb = b.limbs();
      typename CppInt1::limb_pointer       pr = result.limbs();
      //typename CppInt1::limb_pointer       pr_end = pr + m;

      if (as < bs)
         swap(pa, pb);
      // First where a and b overlap:
      unsigned      i = 0;
      unsigned char carry = 0;
#if defined(BOOST_MSVC) && !defined(BOOST_HAS_INT128) && defined(_M_X64)
      for (; i + 8 <= m; i += 8)
      {
         carry = _addcarry_u64(carry, *(unsigned long long*)(pa + i + 0), *(unsigned long long*)(pb + i + 0), (unsigned long long*)(pr + i));
         carry = _addcarry_u64(carry, *(unsigned long long*)(pa + i + 2), *(unsigned long long*)(pb + i + 2), (unsigned long long*)(pr + i + 2));
         carry = _addcarry_u64(carry, *(unsigned long long*)(pa + i + 4), *(unsigned long long*)(pb + i + 4), (unsigned long long*)(pr + i + 4));
         carry = _addcarry_u64(carry, *(unsigned long long*)(pa + i + 6), *(unsigned long long*)(pb + i + 6), (unsigned long long*)(pr + i + 6));
      }
#else
      for (; i + 4 <= m; i += 4)
      {
         carry = ::boost::multiprecision::detail::addcarry_limb(carry, pa[i + 0], pb[i + 0], pr + i);
         carry = ::boost::multiprecision::detail::addcarry_limb(carry, pa[i + 1], pb[i + 1], pr + i + 1);
         carry = ::boost::multiprecision::detail::addcarry_limb(carry, pa[i + 2], pb[i + 2], pr + i + 2);
         carry = ::boost::multiprecision::detail::addcarry_limb(carry, pa[i + 3], pb[i + 3], pr + i + 3);
      }
#endif
      for (; i < m; ++i)
         carry = ::boost::multiprecision::detail::addcarry_limb(carry, pa[i], pb[i], pr + i);
      for (; i < x && carry; ++i)
         carry = ::boost::multiprecision::detail::addcarry_limb(carry, pa[i], 0, pr + i);
      if (i == x && carry)
      {
         // We overflowed, need to add one more limb:
         result.resize(x + 1, x + 1);
         if (result.size() > x)
            result.limbs()[x] = carry;
      }
      else
         std::copy(pa + i, pa + x, pr + i);
      result.normalize();
      result.sign(a.sign());
   }
}

template <class CppInt1, class CppInt2, class CppInt3>
inline BOOST_MP_CXX14_CONSTEXPR void subtract_unsigned(CppInt1& result, const CppInt2& a, const CppInt3& b) BOOST_MP_NOEXCEPT_IF(is_non_throwing_cpp_int<CppInt1>::value)
{
#ifndef BOOST_MP_NO_CONSTEXPR_DETECTION
   if (1 || BOOST_MP_IS_CONST_EVALUATED(a.size()))
   {
      subtract_unsigned_constexpr(result, a, b);
   }
   else
#endif
   {
      using std::swap;

      // Nothing fancy, just let uintmax_t take the strain:
      unsigned         m(0), x(0);
      minmax(a.size(), b.size(), m, x);
      //
      // special cases for small limb counts:
      //
      if (x == 1)
      {
         bool      s = a.sign();
         limb_type al = *a.limbs();
         limb_type bl = *b.limbs();
         if (bl > al)
         {
            ::boost::multiprecision::std_constexpr::swap(al, bl);
            s = !s;
         }
         result = al - bl;
         result.sign(s);
         return;
      }
      // This isn't used till later, but comparison has to occur before we resize the result,
      // as that may also resize a or b if this is an inplace operation:
      int c = a.compare_unsigned(b);
      // Set up the result vector:
      result.resize(x, x);
      // Now that a, b, and result are stable, get pointers to their limbs:
      typename CppInt2::const_limb_pointer pa = a.limbs();
      typename CppInt3::const_limb_pointer pb = b.limbs();
      typename CppInt1::limb_pointer       pr = result.limbs();
      bool                                 swapped = false;
      if (c < 0)
      {
         swap(pa, pb);
         swapped = true;
      }
      else if (c == 0)
      {
         result = static_cast<limb_type>(0);
         return;
      }

      unsigned i = 0;
      unsigned char borrow = 0;
      // First where a and b overlap:
      for(; i + 4 <= m; i += 4)
      {
         borrow = boost::multiprecision::detail::subborrow_limb(borrow, pa[i], pb[i], pr + i);
         borrow = boost::multiprecision::detail::subborrow_limb(borrow, pa[i + 1], pb[i + 1], pr + i + 1);
         borrow = boost::multiprecision::detail::subborrow_limb(borrow, pa[i + 2], pb[i + 2], pr + i + 2);
         borrow = boost::multiprecision::detail::subborrow_limb(borrow, pa[i + 3], pb[i + 3], pr + i + 3);
      }
      for (; i < m; ++i)
         borrow = boost::multiprecision::detail::subborrow_limb(borrow, pa[i], pb[i], pr + i);

      // Now where only a has digits, only as long as we've borrowed:
      while (borrow && (i < x))
      {
         borrow = boost::multiprecision::detail::subborrow_limb(borrow, pa[i], 0, pr + i);
         ++i;
      }
      // Any remaining digits are the same as those in pa:
      if ((x != i) && (pa != pr))
         std_constexpr::copy(pa + i, pa + x, pr + i);
      BOOST_ASSERT(0 == borrow);

      //
      // We may have lost digits, if so update limb usage count:
      //
      result.normalize();
      result.sign(a.sign());
      if (swapped)
         result.negate();
   }  // constepxr.
}

#else

template <class CppInt1, class CppInt2, class CppInt3>
inline BOOST_MP_CXX14_CONSTEXPR void add_unsigned(CppInt1& result, const CppInt2& a, const CppInt3& b) BOOST_MP_NOEXCEPT_IF(is_non_throwing_cpp_int<CppInt1>::value)
{
   add_unsigned_constexpr(result, a, b);
}

template <class CppInt1, class CppInt2, class CppInt3>
inline BOOST_MP_CXX14_CONSTEXPR void subtract_unsigned(CppInt1& result, const CppInt2& a, const CppInt3& b) BOOST_MP_NOEXCEPT_IF(is_non_throwing_cpp_int<CppInt1>::value)
{
   subtract_unsigned_constexpr(result, a, b);
}

#endif

} } }  // namespaces


#endif


