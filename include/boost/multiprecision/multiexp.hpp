//---------------------------------------------------------------------------//
// Copyright (c) 2018-2020 Pavel Kharitonov <ipavrus@nil.foundation>
//
// Distributed under the Boost Software License, Version 1.0
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt
//---------------------------------------------------------------------------//

#ifndef BOOST_MULTIPRECISION_MULTIEXP_HPP
#define BOOST_MULTIPRECISION_MULTIEXP_HPP

#include <boost/multiprecision/modular/modular_adaptor.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>

namespace boost {
namespace multiprecision {

template <typename Backend>
inline int eval_multiexp(const Backend& vec_start, const Backend& scalar_start,
                         const Backend& num_groups, const Backend& bucket_size, 
                         const Backend& n, const Backend& workers_amount)
{ 
    using default_ops::eval_divide;
    using default_ops::eval_ceil;
    using default_ops::eval_multiply;
    using default_ops::eval_subtract;
    using default_ops::eval_add;

    Backend chunk_len;

    eval_divide(chunk_len, n, num_groups);
    eval_ceil(chunk_len, chunk_len);

    //part_res(num_groups);
    typename Backend::const_limb_pointer part_res;

    Backend start, end, one;

    //parallel for j
    for (size_t j = 0; j < num_groups.size(); ++j)
    {
        eval_multiply(start, j, chunk_len);
        eval_add(end, start, chunk_len);
        end = std::min(end, n);
        eval_subtract(end, one);
        
        part_res[j]  = eval_multiexp_subgroup(vec_start, scalar_start, start, end, workers_amount, bucket_size);
    }

    return ResultAggregation(part_res);
}

template <typename Backend>
inline Backend eval_multiexp_subgroup(const Backend& vec_start, const Backend& scalar_start,
                                      const Backend& start, const Backend& end, 
                                      const Backend& workers_amount, const Backend& bucket_size)
{
   using default_ops::eval_add;
   using default_ops::eval_logb;
   using default_ops::eval_divide;
   using default_ops::eval_ceil;

   typename Backend::const_limb_pointer vec_start_pointer = vec_start.limbs();
   typename Backend::const_limb_pointer scalar_start_pointer = scalar_start.limbs();

   Backend res, L, b, c;
   eval_add(res, scalar_start_pointer, start);
   eval_logb(L, *(res));
   eval_ceil(L, L);
   eval_divide(b, L, bucket_size);
   eval_ceil(b, b);
   eval_divide(c, b, workers_amount);
   eval_ceil(c, c);

   //part_sum(workers_amount * c)
   typename Backend::const_limb_pointer part_sum;

   //do parallel for j
   for (size_t j = 0; j < workers_amount.size(); ++j)
   {
      for (size_t k = 0; k <= c.size() - 1; ++k) {

         size_t bucket_start = j * bucket_size * c + k * bucket_size;
         modular_adaptor<Backend> b(2, std::numeric_limits<double>::infinity()), result;
         modular_adaptor<Backend> e(bucket_size, std::numeric_limits<double>::infinity());
         eval_pow(result, b, e);

         //buckets(result.m_base)
         typename Backend::const_limb_pointer buckets;

         for (size_t i = start; i <= end; ++i)
         {
            size_t idx = get_bits(*(scalar_start_pointer + i), bucket_start, bucket_size, 2);
            if (idx > 0)
            {
               buckets[idx - 1] = buckets[idx - 1] + *(vec_start + i);
            }
         }

         Backend acc = std::numeric_limits<Backend>::infinity();
      
         for (size_t i = 0; i <= bucket_size; ++i)
         {
            eval_add(acc, acc, buckets[i]);
            eval_add(part_sum[j], part_sum[j], acc);
         }
      }
   }
   return part_sum;
}


template <typename Backend>
inline Backend get_bits(const Backend& scalar_start, const Backend& start, const Backend& end, const Backend& repr)
{
   typename Backend::const_limb_pointer scalar_start_pointer = scalar_start.limbs();
   Backend res = 0, e = scalar_start;

   for (size_t i = start; i <= end; ++i)
   {
        modular_adaptor<Backend> b(repr, std::numeric_limits<double>::infinity()), result;
        modular_adaptor<Backend> e(i - start, std::numeric_limits<double>::infinity());
        eval_pow(result, b, e);
        if (scalar_start_pointer & 1) {
            res = res + i * result.m_base;
        }
        scalar_start >>= 1; 
   }
   return res;
}


template <typename Backend>
inline int result_aggregation(const Backend& r, const Backend& L)
{
   //part_res(L)
   typename Backend::const_limb_pointer part_res;

   for (size_t i = 0; i <= L.size(); ++i)
   {
      part_res[i] = sum_par(r, i);
   }

   Backend S = std::numeric_limits<double>::infinity();

   for (size_t i = 0; i <= L - 1; ++i)
   {  
      eval_multiplay(S, S, 2);
      eval_add(S, S, part_res[i]);
   }

   return S;
}

template <typename Backend>
inline Backend sum_par(const Backend& vec_start, const Backend& n)
{

   size_t h = std::ceil(n / std::log(n));

   //part_res(h)
   typename Backend::const_limb_pointer part_res;

   //do parallel for i
   for (size_t i = 0; i <= h - 1; ++i)
   {
      for (size_t j = 0; j <= std::log(n) - 1; ++j)
      {
         part_res[i] = part_res[i] + *(vec_start + (i * std::ceil(std::log(n)) + j));
      }
   }

   size_t parallel_boundary = std::ceil(std::log(n));
   size_t m                 = std::ceil(n / std::log(n));

   while (m > parallel_boundary)
   {
      h = std::ceil(std::log((m)));

      //do parallel for i
      for (size_t i = 0; i <= (m / h) - 1; ++i)
      {
         size_t d = h - 1;
         if (i == (m / h) - 1)
         {
            d = m - 1 - i * h;
         }

         for (size_t j = 1; j <= d; ++j)
         {
            part_res[i * h] = part_res[i * h] + part_res[i * h + j];
         }
         part_res[i] = part_res[i * h];
      }

      m = std::ceil((m / h));
   }

   for (size_t i = 1; i <= m - 1; ++i)
   {
      part_res[0] = part_res[0] + part_res[i];
   }

   return part_res[0];
}


template <typename Backend, expression_template_option ExpressionTemplates>
inline typename std::enable_if<number_category<Backend>::value == number_kind_integer, int>::type multiexp(
    const number<Backend, ExpressionTemplates>& vec_start, const number<Backend, ExpressionTemplates>& num_groups, 
    const number<Backend, ExpressionTemplates>& bucket_size, const number<Backend, ExpressionTemplates>& n,
    const number<Backend, ExpressionTemplates>& workers_amount)
{
   return eval_multiexp(vec_start.backend(), vec_start.backend(), num_groups.backend(), bucket_size.backend(), n.backend());
}

template <typename Backend, expression_template_option ExpressionTemplates>
inline typename std::enable_if<number_category<Backend>::value == number_kind_integer, int>::type multiexp(
    const number<Backend, ExpressionTemplates>& vec_start, const number<Backend, ExpressionTemplates>& scalar_start, 
    const number<Backend, ExpressionTemplates>& num_groups, const number<Backend, ExpressionTemplates>& bucket_size,
    const number<Backend, ExpressionTemplates>& n, const number<Backend, ExpressionTemplates>& workers_amount)
{
   return eval_multiexp(vec_start.backend(), scalar_start.backend(), num_groups.backend(), bucket_size.backend(), n.backend(), workers_amount.backend());
}

}
} // namespace boost::multiprecision

#endif //BOOST_MULTIPRECISION_MULTIEXP_HPP