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


template <typename Backend, expression_template_option ExpressionTemplates>
inline Backend sum_par(typename std::vector<Backend>::const_iterator vec_start, const Backend& n)
{
   using default_ops::eval_divide;
   using default_ops::eval_add;
   using default_ops::eval_subtract;
   using default_ops::eval_multiply;

   Backend h, res;

   eval_divide(h, n,  std::log(n));
   eval_ceil(h, h);

   std::vector<number<Backend, ExpressionTemplates>> part_res(h, std::numeric_limits<double>::infinity());

   Backend i = static_cast<limb_type>(0);
   Backend j = static_cast<limb_type>(0);
   //do parallel for i
   while (i != h) 
   {
      while(j != std::log(n))
      {
         eval_ceil(res, std::log(n));
         eval_multiply(res, i, res);
         eval_add(res, j, res);
         eval_add(res, vec_start, res);
         eval_add(part_res[i].backend(), part_res[i].backend(), *(res.backend()));
         ++j;
      }
      ++i;
   }

   Backend parallel_boundary, m;
   eval_ceil(parallel_boundary, std::log(n));
   eval_divide(m, n, std::log(n));
   eval_ceil(m, m);

   while (m != parallel_boundary)
   {
      eval_ceil(h, std::log((m)));

      //do parallel for i
      eval_divide(res, m, h);
      while (i != res)
      {
         Backend d, tmp;
         eval_subtract(d, h, 1);
         eval_divide(tmp, m, h);
         eval_subtract(tmp, tmp, 1);
         if (i == tmp)
         {
            eval_multiply(tmp, i, h);
            eval_subtract(d, m, 1);
            eval_subtract(d, d, tmp);
         }
         while(j != d)
         {
            part_res[i * h].backend() = part_res[i * h].backend() + part_res[i * h + j].backend();
            ++j;
         }
         eval_multiply(tmp, i, h);
         part_res[i].backend() = part_res[tmp].backend();
         ++i;
      }
      eval_divide(m, m, h);
      eval_ceil(m, m);
   }

   i = static_cast<limb_type>(1);
   while(i != m)
   {
      eval_add(part_res[0].backend(), part_res[0].backend(), part_res[i].backend());
   }

   return part_res[0];
}

template <typename Backend, expression_template_option ExpressionTemplates>
inline int result_aggregation(const Backend& r, const Backend& L)
{
   using default_ops::eval_add;
   using default_ops::eval_multiply;

   std::vector<number<Backend, ExpressionTemplates>> part_res(L, std::numeric_limits<double>::infinity());

   Backend i = static_cast<limb_type>(0);

   while(i != L) 
   {
    part_res[i].backend() = sum_par(r, i);
   }

   Backend S = std::numeric_limits<double>::infinity();

   for (size_t i = 0; i <= L - 1; ++i)
   {  
      eval_multiply(S, S, 2);
      eval_add(S, S, part_res[i].backend());
   }

   return S;
}

template <typename Backend>
inline Backend get_bits(const Backend& scalar_start, const Backend& start, const Backend& end, const Backend& repr)
{
   using default_ops::eval_add;

   Backend res = static_cast<limb_type>(0), i = start;
   Backend scalar = scalar_start;
   scalar >>=start;

   while (i <= end)
   {
        if (scalar & 1) {
          modular_adaptor<Backend> b(repr, std::numeric_limits<double>::infinity()), result;
          modular_adaptor<Backend> e(i - start, std::numeric_limits<double>::infinity());
          eval_pow(result, b, e);
          eval_add(res, res, result);
        }
        scalar >>= 1; 
        ++i;
   }
   return res;
}

template <typename Backend, expression_template_option ExpressionTemplates>
inline Backend eval_multiexp_subgroup(typename std::vector<std::pair<number<Backend, ExpressionTemplates>, number<Backend, ExpressionTemplates>>>::const_iterator & vec_start,
                                      typename std::vector<number<Backend, ExpressionTemplates>>::const_iterator & scalar_start,
                                      const Backend& start, const Backend& end, 
                                      const Backend& workers_amount, const Backend& bucket_size)
{
   using default_ops::eval_add;
   using default_ops::eval_logb;
   using default_ops::eval_divide;
   using default_ops::eval_ceil;
   using default_ops::eval_multiply;

   Backend res, L, b, c;
   eval_add(res, scalar_start, start);
   eval_logb(L, *(res.backend()));
   eval_ceil(L, L);
   eval_divide(b, L, bucket_size);
   eval_ceil(b, b);
   eval_divide(c, b, workers_amount);
   eval_ceil(c, c);

   eval_multiply(res, workers_amount, c);
   std::vector<number<Backend, ExpressionTemplates>> part_sum(res, std::numeric_limits<double>::infinity());

   //do parallel for j
   Backend j = static_cast<limb_type>(0);
   while (j != workers_amount) 
   {
      Backend k = static_cast<limb_type>(0);
      while(k != c)
      {  
         Backend bucket_start, res;
         eval_multiply(bucket_start, bucket_size, j);
         eval_multiply(bucket_start, bucket_start, c);
         eval_multiply(res, k, bucket_size);
         eval_add(bucket_start, res);
         modular_adaptor<Backend> b(2, std::numeric_limits<double>::infinity()), result;
         modular_adaptor<Backend> e(bucket_size, std::numeric_limits<double>::infinity());
         eval_pow(result, b, e);

         //buckets(result.m_base)
         std::vector<number<Backend, ExpressionTemplates>> buckets(result.m_base, std::numeric_limits<double>::infinity());

         Backend i = start, tmp;
         while (i != end)
         {
            size_t idx = get_bits(*(scalar_start[i].backend()), bucket_start, bucket_size, 2);
            if (idx > 0)
            {
              eval_add(tmp, vec_start, i);
              buckets[idx - 1].backend() = buckets[idx - 1].backend() + (*(tmp.backend());
            }
            ++i;
         }

         Backend acc = std::numeric_limits<Backend>::infinity();
        
         Backend i = static_cast<limb_type>(0);
         while(i != bucket_size())
         {
            eval_add(acc, acc, buckets[i].backend());
            eval_add(part_sum[j].backend(), part_sum[j].backend(), acc);
            ++i;
         }
         ++k;
      }
      ++j;
   }
   return part_sum;
}

template <typename Backend, expression_template_option ExpressionTemplates>
inline int eval_multiexp(typename std::vector<std::pair<number<Backend, ExpressionTemplates>, number<Backend, ExpressionTemplates>>>::const_iterator & vec_start,
                         typename std::vector<number<Backend, ExpressionTemplates>>::const_iterator & scalar_start,
                         const Backend& num_groups, const Backend& bucket_size, 
                         const Backend& n, const Backend& workers_amount)
{ 
    using default_ops::eval_divide;
    using default_ops::eval_ceil;
    using default_ops::eval_multiply;
    using default_ops::eval_subtract;
    using default_ops::eval_add;

    Backend res, chunk_len;

    eval_divide(res, n, num_groups);
    eval_ceil(chunk_len, res);

    std::vector<number<Backend, ExpressionTemplates>> part_res(num_groups);

    Backend start, end;
    Backend one = static_cast<limb_type>(1);
    Backend j = static_cast<limb_type>(0);

    //parallel for j
    while (j != num_groups)
    {
        eval_multiply(start, j, chunk_len);
        eval_add(end, start, chunk_len);
        end = std::min(end, n);
        eval_subtract(end, one);
        
        part_res[j]  = eval_multiexp_subgroup(vec_start, scalar_start, start, end, workers_amount, bucket_size);
        ++j;
    }

    return ResultAggregation(part_res, bucket_size);
}

template <typename Backend, expression_template_option ExpressionTemplates>
inline typename std::enable_if<number_category<Backend>::value == number_kind_integer, int>::type multiexp(
    typename std::vector<std::pair<number<Backend, ExpressionTemplates>, number<Backend, ExpressionTemplates>>>::const_iterator &vec_start,
    typename std::vector<number<Backend, ExpressionTemplates>>::const_iterator &scalar_start,
    const number<Backend, ExpressionTemplates>& num_groups, const number<Backend, ExpressionTemplates>& bucket_size,
    const number<Backend, ExpressionTemplates>& n, const number<Backend, ExpressionTemplates>& workers_amount)
{
   return eval_multiexp(vec_start, scalar_start, num_groups.backend(), bucket_size.backend(), n.backend(), workers_amount.backend());
}

}
} // namespace boost::multiprecision

#endif //BOOST_MULTIPRECISION_MULTIEXP_HPP