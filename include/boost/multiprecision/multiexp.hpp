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

namespace boost {
namespace multiprecision {

template <typename Backend>
inline Backend eval_multiexp(typename std::vector<Backend>::const_iterator vec_start,
                 		  	  typename std::vector<Backend>::const_iterator scalar_start,
                 		  	  const Backend& num_groups, const Backend& bucket_size, 
                 		  	  const Backend& n, const Backend& workers_amount)
{	
    size_t chunk_len = std::ceil(n  / num_groups);

    std::vector<Backend> part_res(num_groups);

    //do parallel for j
    for (size_t j = 0; j < num_groups; ++j)
    {
    	size_t start = j * chunk_len;
        size_t end   = std::min(start + chunk_len - 1, n - 1);
        part_res[j]  = eval_multiexp_subgroup(vec_start, scalar_start, start, end, workers_amount, bucket_size);
    }

    return ResultAggregation(part_res);
}

template <typename Backend>
inline Backend eval_multiexp_subgroup(typename std::vector<Backend>::const_iterator vec_start,
                 		  			   typename std::vector<Backend>::const_iterator scalar_start,
                    	  			   const Backend& start, const Backend& end, 
                          			   const Backend& workers_amount, const Backend& bucket_size)
{
   size_t L = std::log2(*(scalar_start + start));
   size_t b = std::ceil(L / bucket_size);
   size_t c = std::ceil(b / workers_amount);

   typename std::vector<Backend> part_sum(workers_amount * c, std::numeric_limits<double>::infinity());

   //do parallel for j
   for (size_t j = 0; j < workers_amount; ++j)
   {
      for (size_t k = 0; k <= c - 1; ++k) {

         size_t bucket_start = j * bucket_size * c + k * bucket_size;
         modular_adaptor<Backend> b(2, std::numeric_limits<double>::infinity()), result;
         modular_adaptor<Backend> e(bucket_size, std::numeric_limits<double>::infinity());
         eval_pow(result, b, e);
         typename std::vector<Backend> buckets(result.m_base, std::numeric_limits<double>::infinity());

         for (size_t i = start; i <= end; ++i)
         {
            size_t idx = get_bits(*(scalar_start + i), bucket_start, bucket_size, 2);
            if (idx > 0)
            {
               buckets[idx - 1] = buckets[idx - 1] + *(vec_start + i);
            }
         }

         size_t acc = std::numeric_limits<double>::infinity();
      
         for (size_t i = 0; i <= bucket_size; ++i)
         {
            acc         = acc + buckets[i];
            part_sum[j] = part_sum[j] + acc;
         }
      }
   }
   return part_sum;
}

template <typename Backend>
inline Backend get_bits(const Backend& scalar_start, const Backend& start, const Backend& end, const Backend& repr)
{
   Backend res = 0, e = scalar_start;

   for (size_t i = start; i <= end; ++i)
   {
        modular_adaptor<Backend> b(repr, std::numeric_limits<double>::infinity()), result;
        modular_adaptor<Backend> e(i - start, std::numeric_limits<double>::infinity());
        eval_pow(result, b, e);
        if (scalar_start & 1) {
      	    res = res + i * result.m_base;
      	}
      	scalar_start >>= 1;	
   }
   return res;
}

template <typename Backend>
inline Backend result_aggregation(typename std::vector<Backend> r, const Backend& L)
{
   typename std::vector<Backend> part_res(L, std::numeric_limits<double>::infinity());

   for (size_t i = 0; i <= L; ++i)
   {
      part_res[i] = sum_par(r, i);
   }

   Backend S = std::numeric_limits<double>::infinity();

   for (size_t i = 0; i <= L - 1; ++i)
   {
      S = S * 2;
      S = S + part_res[i];
   }

   return S;
}

template <typename Backend>
inline Backend sum_par(typename std::vector<Backend>::const_iterator vec_start, const Backend& n)
{
   size_t                        h = std::ceil(n / std::log(n));
   typename std::vector<Backend> part_res(h, std::numeric_limits<double>::infinity());

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