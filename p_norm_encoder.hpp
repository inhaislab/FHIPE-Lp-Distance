//
// Created by DongHyeon Ryu on 22. 12. 5..
//

#ifndef LIB_FFHIPE_P_NORM_ENCODER_HPP
#define LIB_FFHIPE_P_NORM_ENCODER_HPP

#include <ffhipe.h>
#include <vector>
#include "encoder.hpp"

namespace ffhipe {
    template<typename T>
    class PNormEncoder : public ffhipe::CEncoder<T> {
    public:
        explicit PNormEncoder(ULONG _p): p(_p){};

        ULONG getN(ULONG n) override {
            return n;
        }
        /**
         *
         * @param x a normal state vector
         * @param encodedX
         * @return encodedX
         */
        VECTOR<MSG> & encodeX(const std::vector<T> &x, VECTOR<MSG> &encodedX) override {

            long n = x.size();
            long encoded_dimension = (p - 1) * n + 2;
            encodedX.SetLength(encoded_dimension);

            std::vector<T> _x = std::vector<T>(x);

            long coef;
            encodedX[0] = get_signed_binomial_coefficient(p, p);
            for (ULONG i = 0; i < p-1; i++)
            {
                coef = get_signed_binomial_coefficient(p, p - i - 1);
                for (ULONG j = 0; j < n; j++)
                {
                    encodedX[i * n + j + 1] = _x[j] * coef;
                    _x[j] *= x[j];
                }
            }
            for (ULONG j = 0; j < n; j++)
            {
                    encodedX[encoded_dimension - 1] += _x[j];
            }

            return encodedX;
        }

        /**
         *
         * @param y current state vector to be examined if it is abnormal or not.
         * @param encodedY
         * @return encodedY
         */
        VECTOR<MSG> & encodeY(const std::vector<T> &y, VECTOR<MSG> &encodedY) override {

            long n = y.size();
            long encoded_dimension = (p - 1) * n + 2;
            encodedY.SetLength(encoded_dimension);

            std::vector<T> _y = std::vector<T>(y);

            encodedY[encoded_dimension-1] = 1;
            for (ULONG i = p-1; i > 0; i--)
            {
                for (ULONG j = 0; j < n ; j++)
                {
                    encodedY[i * n - j] = _y[n - j - 1];
                    _y[n - j - 1] *= y[n - j - 1];
                }
            }
            for (ULONG j = 0; j < n; j++)
            {
                    encodedY[0] += _y[j];
            }

            return encodedY;
        }

    private:

        const ULONG p;
        std::vector<LONG> nC;

        /**
         * Time Complexity is O(k). However once calculate nCk the intermediate values are all stored.
         * Space complexity is O(k).
         */
        LONG get_signed_binomial_coefficient(const ULONG n, const ULONG k)
        {
            if(this->nC.empty()) {
                nC = std::vector<LONG>(n + 1);

                nC[0] = nC[n] = 1;

                for (ULONG i = 1; i < n; ++i) {
                    nC[i] = nC[i - 1] * (n + 1 - i) / i ;
                }
            }

            LONG sign = (k % 2 == 0 ? 1 : -1);

            return nC[k] * sign;
        }
    };
}

#endif //LIB_FFHIPE_P_NORM_ENCODER_HPP
