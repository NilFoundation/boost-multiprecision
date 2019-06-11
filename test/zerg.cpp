//
// Created by Zerg1996 on 2019-04-27.
//

//#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/montgomery_int/montgomery_int.hpp>
#include <boost/multiprecision/montgomery_int/montgomery_params.hpp>
#include <iostream>


int main()
{

    using namespace boost::multiprecision;
    using default_ops::eval_msb;

    cpp_int a(7);
    nil::crypto3::montgomery_params<cpp_int> tt(a);
    std::cout << "-----" << std::endl;
    //montgomery_int x(9299, tt);
    //montgomery_int x2(1315541, tt);
    montgomery_int x(3, tt);
    std::cout << "X1:" << x << std::endl;
    //montgomery_int x2(b, tt);
    montgomery_int x2(6, tt);
    std::cout << "X2:" << x2 << std::endl;

    x = x * x2;
    std::cout << "Mult:" << x << std::endl;
    return 0;
}
