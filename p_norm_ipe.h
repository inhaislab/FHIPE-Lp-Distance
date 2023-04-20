//
// Created by DongHyeon Yu on 22. 11. 28..
//

#ifndef LIB_FFHIPE_P_NORM_IPE_H
#define LIB_FFHIPE_P_NORM_IPE_H

#include <initializer_list>
#include "../libffhipe/ffhipe.h"


class p_norm_ipe {
public:

    /**
     * @param n dimension of input vector.
     * @param m max range of vector's elements.
     * @param p norm degree.
     */
    static void decrypt_value_test(ULONG n, ULONG m, ULONG p);
    static void encrypt(ULONG n, ULONG m, ULONG p);
    static void keygen(ULONG n, ULONG m, ULONG p);
    static void decrypt_performance_test(const std::initializer_list<ULONG>& n_arr, ULONG m, const std::initializer_list<ULONG>& p_arr, ULONG times);
    static void encrypt_keygen_performance_test(const std::initializer_list<ULONG>& n_arr, ULONG m, const std::initializer_list<ULONG>& p_arr, ULONG times);
};

#endif //LIB_FFHIPE_P_NORM_IPE_H
