#ifndef CRYPTO3_MP_MONTGOMERY_INT_MISC_HPP
#define CRYPTO3_MP_MONTGOMERY_INT_MISC_HPP

#include <boost/multiprecision/detail/bitscan.hpp> // lsb etc
#include <boost/integer/common_factor_rt.hpp> // gcd/lcm
#include <boost/functional/hash_fwd.hpp>


#ifdef BOOST_MSVC
#pragma warning(push)
#pragma warning(disable:4702)
#pragma warning(disable:4127) // conditional expression is constant
#pragma warning(disable:4146) // unary minus operator applied to unsigned type, result still unsigned
#endif


namespace boost {
    namespace multiprecision {
        namespace backends {


            template<
                    class R, unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            inline typename enable_if_c<is_integral<R>::value && !is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value, void>

            ::type eval_convert_to(R *result, const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                    Allocator1> &backend)

            BOOST_MP_NOEXCEPT_IF((is_non_throwing_cpp_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value)) {
                typedef mpl::int_<Checked1> checked_type;
                check_in_range<R>(backend, checked_type());

                *result = static_cast<R>(backend.limbs()[0]);
                unsigned shift = montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1>::limb_bits;
                for (unsigned i = 1; (i < backend.size()) && (shift < static_cast<unsigned>(std::numeric_limits<R>::digits)); ++i) {
                    *result += static_cast<R>(backend.limbs()[i]) << shift;
                    shift += montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1>::limb_bits;
                }
                if (backend.sign()) {
                    check_is_negative(boost::is_signed<R>());

                    *result = negate_integer(*result, boost::is_signed<R>());
                }
            }

            template<typename R, unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            inline typename enable_if_c<is_floating_point<R>::value && !is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value, void>

            ::type eval_convert_to(R *result, const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                    Allocator1, ParamsBackend1> &backend)

            BOOST_MP_NOEXCEPT_IF(is_arithmetic<R>::value) {
                typename montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                        Allocator1, ParamsBackend1>::const_limb_pointer p = backend.limbs();
                unsigned shift = montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1>::limb_bits;
                *result = static_cast<R>(*p);
                for (unsigned i = 1; i < backend.size(); ++i) {
                    *result += static_cast<R>(std::ldexp(static_cast<long double>(p[i]), shift));
                    shift += montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1>::limb_bits;
                }
                if (backend.sign()) {
                    *result = -*result;
                }
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            BOOST_MP_FORCEINLINE typename enable_if_c<!is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value, bool>

            ::type eval_is_zero(const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                    Allocator1, ParamsBackend1> &val) BOOST_NOEXCEPT {
                return (val.size() == 1) && (val.limbs()[0] == 0);
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            BOOST_MP_FORCEINLINE typename enable_if_c<!is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value, int>

            ::type eval_get_sign(const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                    Allocator1, ParamsBackend1> &val) BOOST_NOEXCEPT {
                return eval_is_zero(val) ? 0 : val.sign() ? -1 : 1;
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            BOOST_MP_FORCEINLINE typename enable_if_c<!is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value>

            ::type eval_abs(montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &result,
                            const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &val)

            BOOST_MP_NOEXCEPT_IF((is_non_throwing_cpp_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value)) {
                result = val;
                result.sign(false);
            }

//
// Get the location of the least-significant-bit:
//
            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            inline typename enable_if_c<!is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value, unsigned>

            ::type eval_lsb(const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &a) {
                using default_ops::eval_get_sign;
                if (eval_get_sign(a) == 0) {
                    BOOST_THROW_EXCEPTION(std::range_error("No bits were set in the operand."));
                }
                if (a.sign()) {
                    BOOST_THROW_EXCEPTION(std::range_error(
                            "Testing individual bits in negative values is not supported - results are undefined."));
                }

                //
                // Find the index of the least significant limb that is non-zero:
                //
                unsigned index = 0;
                while (!a.limbs()[index] && (index < a.size())) {
                    ++index;
                }
                //
                // Find the index of the least significant bit within that limb:
                //
                unsigned result = boost::multiprecision::detail::find_lsb(a.limbs()[index]);

                return result +
                       index * montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1>::limb_bits;
            }

//
// Get the location of the most-significant-bit:
//
            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            inline typename enable_if_c<!is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value, unsigned>

            ::type eval_msb_imp(const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &a) {
                //
                // Find the index of the most significant bit that is non-zero:
                //
                return (a.size() - 1) *
                       montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1>::limb_bits +
                       boost::multiprecision::detail::find_msb(a.limbs()[a.size() - 1]);
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            inline typename enable_if_c<!is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value, unsigned>

            ::type eval_msb(const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &a) {
                using default_ops::eval_get_sign;
                if (eval_get_sign(a) == 0) {
                    BOOST_THROW_EXCEPTION(std::range_error("No bits were set in the operand."));
                }
                if (a.sign()) {
                    BOOST_THROW_EXCEPTION(std::range_error(
                            "Testing individual bits in negative values is not supported - results are undefined."));
                }
                return eval_msb_imp(a);
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            inline typename enable_if_c<!is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value, bool>

            ::type eval_bit_test(const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &val,
                                 unsigned index) BOOST_NOEXCEPT {
                unsigned offset =
                        index / montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1>::limb_bits;
                unsigned shift =
                        index % montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1>::limb_bits;
                limb_type mask = shift ? limb_type(1u) << shift : limb_type(1u);
                if (offset >= val.size()) {
                    return false;
                }
                return val.limbs()[offset] & mask ? true : false;
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            inline typename enable_if_c<!is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value>

            ::type eval_bit_set(montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &val,
                                unsigned index) {
                unsigned offset =
                        index / montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1>::limb_bits;
                unsigned shift =
                        index % montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1>::limb_bits;
                limb_type mask = shift ? limb_type(1u) << shift : limb_type(1u);
                if (offset >= val.size()) {
                    unsigned os = val.size();
                    val.resize(offset + 1, offset + 1);
                    if (offset >= val.size()) {
                        return;
                    }  // fixed precision overflow
                    for (unsigned i = os; i <= offset; ++i) {
                        val.limbs()[i] = 0;
                    }
                }
                val.limbs()[offset] |= mask;
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            inline typename enable_if_c<!is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value>

            ::type eval_bit_unset(montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &val,
                                  unsigned index) BOOST_NOEXCEPT {
                unsigned offset =
                        index / montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1>::limb_bits;
                unsigned shift =
                        index % montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1>::limb_bits;
                limb_type mask = shift ? limb_type(1u) << shift : limb_type(1u);
                if (offset >= val.size()) {
                    return;
                }
                val.limbs()[offset] &= ~mask;
                val.normalize();
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            inline typename enable_if_c<!is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value>

            ::type eval_bit_flip(montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &val,
                                 unsigned index) {
                unsigned offset =
                        index / montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1>::limb_bits;
                unsigned shift =
                        index % montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1>::limb_bits;
                limb_type mask = shift ? limb_type(1u) << shift : limb_type(1u);
                if (offset >= val.size()) {
                    unsigned os = val.size();
                    val.resize(offset + 1, offset + 1);
                    if (offset >= val.size()) {
                        return;
                    }  // fixed precision overflow
                    for (unsigned i = os; i <= offset; ++i) {
                        val.limbs()[i] = 0;
                    }
                }
                val.limbs()[offset] ^= mask;
                val.normalize();
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            inline typename enable_if_c<!is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value>

            ::type eval_qr(const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &x,
                           const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &y,
                           montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &q,
                           montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &r)

            BOOST_MP_NOEXCEPT_IF((is_non_throwing_cpp_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value)) {
                divide_unsigned_helper(&q, x, y, r);
                q.sign(x.sign() != y.sign());
                r.sign(x.sign());
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            inline typename enable_if_c<!is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value>

            ::type eval_qr(const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &x,
                           limb_type y, montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &q,
                           montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &r)

            BOOST_MP_NOEXCEPT_IF((is_non_throwing_cpp_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value)) {
                divide_unsigned_helper(&q, x, y, r);
                q.sign(x.sign());
                r.sign(x.sign());
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1, class U>
            inline typename enable_if_c<is_integral<U>::value>::type eval_qr(
                    const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &x, U y,
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &q,
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &r)

            BOOST_MP_NOEXCEPT_IF((is_non_throwing_cpp_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value)) {
                using default_ops::eval_qr;
                montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> t;
                t = y;
                eval_qr(x, t, q, r);
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1, class Integer>
            inline typename enable_if_c<is_unsigned<Integer>::value && !is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value, Integer>

            ::type eval_integer_modulus(
                    const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &x, Integer val) {
                if ((sizeof(Integer) <= sizeof(limb_type)) || (val <= (std::numeric_limits<limb_type>::max)())) {
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> d;
                    divide_unsigned_helper(static_cast<montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                            Allocator1> *>(0), x, static_cast<limb_type>(val), d);
                    return d.limbs()[0];
                } else {
                    return default_ops::eval_integer_modulus(x, val);
                }
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1, class Integer>
            BOOST_MP_FORCEINLINE typename enable_if_c<is_signed<Integer>::value && !is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value, Integer>

            ::type eval_integer_modulus(
                    const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &x, Integer val) {
                return eval_integer_modulus(x, boost::multiprecision::detail::unsigned_abs(val));
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            inline typename enable_if_c<!is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value>

            ::type eval_gcd(montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &result,
                            const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &a,
                            limb_type v) {
                using default_ops::eval_lsb;
                using default_ops::eval_is_zero;
                using default_ops::eval_get_sign;

                int shift;

                montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> u(a);

                int s = eval_get_sign(u);

                /* GCD(0,x) := x */
                if (s < 0) {
                    u.negate();
                } else if (s == 0) {
                    result = v;
                    return;
                }
                if (v == 0) {
                    result = u;
                    return;
                }

                /* Let shift := lg K, where K is the greatest power of 2
                dividing both u and v. */

                unsigned us = eval_lsb(u);
                unsigned vs = boost::multiprecision::detail::find_lsb(v);
                shift = (std::min)(us, vs);
                eval_right_shift(u, us);
                if (vs) {
                    v >>= vs;
                }

                do {
                    /* Now u and v are both odd, so diff(u, v) is even.
                    Let u = min(u, v), v = diff(u, v)/2. */
                    if (u.size() <= 2) {
                        if (u.size() == 1) {
                            v = integer_gcd_reduce(*u.limbs(), v);
                        } else {
                            double_limb_type i;
                            i = u.limbs()[0] |
                                (static_cast<double_limb_type>(u.limbs()[1]) << sizeof(limb_type) * CHAR_BIT);
                            v = static_cast<limb_type>(integer_gcd_reduce(i, static_cast<double_limb_type>(v)));
                        }
                        break;
                    }
                    eval_subtract(u, v);
                    us = eval_lsb(u);
                    eval_right_shift(u, us);
                } while (true);

                result = v;
                eval_left_shift(result, shift);
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1, class Integer>
            inline typename enable_if_c<is_unsigned<Integer>::value && (sizeof(Integer) <= sizeof(limb_type)) &&
                                        !is_trivial_montgomery_int<
                                                montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                                                        Allocator1> >::value>

            ::type eval_gcd(montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &result,
                            const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &a,
                            const Integer &v) {
                eval_gcd(result, a, static_cast<limb_type>(v));
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1, class Integer>
            inline typename enable_if_c<is_signed<Integer>::value && (sizeof(Integer) <= sizeof(limb_type)) &&
                                        !is_trivial_montgomery_int<
                                                montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                                                        Allocator1> >::value>

            ::type eval_gcd(montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &result,
                            const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &a,
                            const Integer &v) {
                eval_gcd(result, a, static_cast<limb_type>(v < 0 ? -v : v));
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            inline typename enable_if_c<!is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value>

            ::type eval_gcd(montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &result,
                            const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &a,
                            const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &b) {
                using default_ops::eval_lsb;
                using default_ops::eval_is_zero;
                using default_ops::eval_get_sign;

                if (a.size() == 1) {
                    eval_gcd(result, b, *a.limbs());
                    return;
                }
                if (b.size() == 1) {
                    eval_gcd(result, a, *b.limbs());
                    return;
                }

                int shift;

                montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> u(a), v(b);

                int s = eval_get_sign(u);

                /* GCD(0,x) := x */
                if (s < 0) {
                    u.negate();
                } else if (s == 0) {
                    result = v;
                    return;
                }
                s = eval_get_sign(v);
                if (s < 0) {
                    v.negate();
                } else if (s == 0) {
                    result = u;
                    return;
                }

                /* Let shift := lg K, where K is the greatest power of 2
                dividing both u and v. */

                unsigned us = eval_lsb(u);
                unsigned vs = eval_lsb(v);
                shift = (std::min)(us, vs);
                eval_right_shift(u, us);
                eval_right_shift(v, vs);

                do {
                    /* Now u and v are both odd, so diff(u, v) is even.
                    Let u = min(u, v), v = diff(u, v)/2. */
                    s = u.compare(v);
                    if (s > 0) {
                        u.swap(v);
                    }
                    if (s == 0) {
                        break;
                    }
                    if (v.size() <= 2) {
                        if (v.size() == 1) {
                            u = integer_gcd_reduce(*v.limbs(), *u.limbs());
                        } else {
                            double_limb_type i, j;
                            i = v.limbs()[0] |
                                (static_cast<double_limb_type>(v.limbs()[1]) << sizeof(limb_type) * CHAR_BIT);
                            j = (u.size() == 1) ? *u.limbs() : u.limbs()[0] |
                                                               (static_cast<double_limb_type>(u.limbs()[1])
                                                                       << sizeof(limb_type) * CHAR_BIT);
                            u = integer_gcd_reduce(i, j);
                        }
                        break;
                    }
                    eval_subtract(v, u);
                    vs = eval_lsb(v);
                    eval_right_shift(v, vs);
                } while (true);

                result = u;
                eval_left_shift(result, shift);
            }

//
// Now again for trivial backends:
//
            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            BOOST_MP_FORCEINLINE typename enable_if_c<is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value>

            ::type eval_gcd(montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &result,
                            const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &a,
                            const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                                    Allocator1> &b) BOOST_NOEXCEPT {
                *result.limbs() = boost::integer::gcd(*a.limbs(), *b.limbs());
            }

// This one is only enabled for unchecked cpp_int's, for checked int's we need the checking in the default version:
            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            BOOST_MP_FORCEINLINE typename enable_if_c<is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value &&
                                                      (Checked1 == unchecked)>

            ::type eval_lcm(montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &result,
                            const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &a,
                            const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &b)

            BOOST_MP_NOEXCEPT_IF((is_non_throwing_cpp_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value)) {
                *result.

                        limbs() = boost::integer::lcm(*a.limbs(), *b.limbs());

                result.

                        normalize(); // result may overflow the specified number of bits
            }

            template<
                    class R, unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            inline typename enable_if_c<is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value &&
                                        is_signed_number<montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                                                Allocator1> >::value && boost::is_convertible<
                    typename montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                            Allocator1>::local_limb_type, R>::value>

            ::type eval_convert_to(R *result, const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                    Allocator1> &val) {
                typedef typename common_type<R, typename montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                        Allocator1>::local_limb_type>::type common_type;
                if (std::numeric_limits<R>::is_specialized && (static_cast<common_type>(*val.limbs()) >
                                                               static_cast<common_type>((std::numeric_limits<
                                                                       R>::max)()))) {
                    if (val.isneg()) {
                        if (static_cast<common_type>(*val.limbs()) >
                            -static_cast<common_type>((std::numeric_limits<R>::min)())) {
                            conversion_overflow(typename montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                                    Allocator1>::checked_type());
                        }
                        *result = (std::numeric_limits<R>::min)();
                    } else {
                        conversion_overflow(typename montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                                Allocator1>::checked_type());
                        *result = (std::numeric_limits<R>::max)();
                    }
                } else {
                    *result = static_cast<R>(*val.limbs());
                    if (val.isneg()) {
                        check_is_negative(mpl::bool_<is_signed_number<R>::value ||
                                                     (number_category<R>::value == number_kind_floating_point)>());
                        *result = negate_integer(*result, mpl::bool_<is_signed_number<R>::value ||
                                                                     (number_category<R>::value ==
                                                                      number_kind_floating_point)>());
                    }
                }
            }

            template<
                    class R, unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            inline typename enable_if_c<is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value &&
                                        is_unsigned_number<
                                                montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                                                        Allocator1> >::value && boost::is_convertible<
                    typename montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                            Allocator1>::local_limb_type, R>::value>

            ::type eval_convert_to(R *result, const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                    Allocator1> &val) {
                typedef typename common_type<R, typename montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                        Allocator1>::local_limb_type>::type common_type;
                if (std::numeric_limits<R>::is_specialized && (static_cast<common_type>(*val.limbs()) >
                                                               static_cast<common_type>((std::numeric_limits<
                                                                       R>::max)()))) {
                    conversion_overflow(typename montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                            Allocator1>::checked_type());
                    *result = (std::numeric_limits<R>::max)();
                } else {
                    *result = static_cast<R>(*val.limbs());
                }
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            inline typename enable_if_c<is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value, unsigned>

            ::type eval_lsb(const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &a) {
                using default_ops::eval_get_sign;
                if (eval_get_sign(a) == 0) {
                    BOOST_THROW_EXCEPTION(std::range_error("No bits were set in the operand."));
                }
                if (a.sign()) {
                    BOOST_THROW_EXCEPTION(std::range_error(
                            "Testing individual bits in negative values is not supported - results are undefined."));
                }
                //
                // Find the index of the least significant bit within that limb:
                //
                return boost::multiprecision::detail::find_lsb(*a.limbs());
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            inline typename enable_if_c<is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value, unsigned>

            ::type eval_msb_imp(const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &a) {
                //
                // Find the index of the least significant bit within that limb:
                //
                return boost::multiprecision::detail::find_msb(*a.limbs());
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            inline typename enable_if_c<is_trivial_montgomery_int<
                    montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> >::value, unsigned>

            ::type eval_msb(const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1, Allocator1, ParamsBackend1> &a) {
                using default_ops::eval_get_sign;
                if (eval_get_sign(a) == 0) {
                    BOOST_THROW_EXCEPTION(std::range_error("No bits were set in the operand."));
                }
                if (a.sign()) {
                    BOOST_THROW_EXCEPTION(std::range_error(
                            "Testing individual bits in negative values is not supported - results are undefined."));
                }
                return eval_msb_imp(a);
            }

            template<unsigned MinBits1, unsigned MaxBits1, cpp_integer_type SignType1, cpp_int_check_type Checked1,
                    class Allocator1, typename ParamsBackend1>
            inline std::size_t hash_value(const montgomery_int_backend<MinBits1, MaxBits1, SignType1, Checked1,
                    Allocator1, ParamsBackend1> &val) BOOST_NOEXCEPT {
                std::size_t result = 0;
                for (unsigned i = 0; i < val.size(); ++i) {
                    boost::hash_combine(result, val.limbs()[i]);
                }
                boost::hash_combine(result, val.sign());
                return result;
            }

#ifdef BOOST_MSVC
#pragma warning(pop)
#endif

        }
    }
} // namespaces

#endif
