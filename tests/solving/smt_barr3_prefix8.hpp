/***************************************************************************
 *  Barr3 2-8-8-1 scaling fixture derived from the published 2-64-64-1 model.
 *
 *  Source checkpoint SHA-256:
 *  788fd56f21abb0cb5b83206b12d4d8d0ed718896c2140222a8fd129475f6c213
 *
 *  Extraction rule: first 8 first-layer neurons, top-left 8x8 block of the
 *  second layer, first 8 second-layer biases and first 8 output weights; the
 *  published output bias is retained. This prefix network is NOT the published
 *  Barrier-3 certificate. It is a deterministic Ariadne scaling fixture with
 *  the same tanh/affine computational structure and parameter distribution.
 *
 *  Concatenated float32 tensor bytes SHA-256:
 *  9431d092372b0bc57ffcf4b651d6bcd2f664c29d79988f7055880f79c4c53aec
 ***************************************************************************/

#ifndef ARIADNE_TEST_SMT_BARR3_PREFIX8_HPP
#define ARIADNE_TEST_SMT_BARR3_PREFIX8_HPP

#include <array>

namespace Ariadne::TestBarr3Prefix8 {

inline constexpr SizeType width = 8u;

inline constexpr std::array<double,16> w1 = {
    0x1.18260a0000000p-3, -0x1.c375000000000p-1, -0x1.670bf00000000p-3, -0x1.430f520000000p-1,
    -0x1.edda9c0000000p-2, 0x1.6fa2680000000p-3, 0x1.cf4fca0000000p-2, -0x1.eba0ce0000000p-3,
    -0x1.dcf21a0000000p-5, -0x1.a6de200000000p-2, -0x1.d3cd9a0000000p-8, 0x1.9455b60000000p-1,
    0x1.e2bb6c0000000p-2, -0x1.587c060000000p-1, 0x1.7a6e8a0000000p-5, -0x1.3759b20000000p-1,
};

inline constexpr std::array<double,8> b1 = {
    0x1.a88fe60000000p-2, -0x1.80a4640000000p-3, 0x1.35baa60000000p-1, 0x1.d18a6a0000000p-3,
    -0x1.67b1400000000p-1, 0x1.5af0700000000p-1, -0x1.e5c5860000000p-3, -0x1.874c320000000p-2,
};

inline constexpr std::array<double,64> w2 = {
    0x1.02bba80000000p-1, -0x1.9d552a0000000p-4, 0x1.92810c0000000p-3, 0x1.7438200000000p-9,
    0x1.c94c9c0000000p-3, -0x1.32a0220000000p-7, 0x1.91a3bc0000000p-2, 0x1.b7391e0000000p-2,
    0x1.a4b9a40000000p-5, 0x1.aff2be0000000p-2, 0x1.aef3620000000p-3, -0x1.93d2960000000p-4,
    0x1.a0e4340000000p-3, 0x1.61e43e0000000p-3, 0x1.734ec80000000p-2, -0x1.74d0920000000p-3,
    0x1.4fc9ae0000000p-2, -0x1.2f40ac0000000p-2, 0x1.e14ce80000000p-3, -0x1.e539200000000p-4,
    -0x1.fc7fd80000000p-3, -0x1.77d1e40000000p-4, -0x1.28329c0000000p-2, 0x1.95a3220000000p-4,
    -0x1.c7cf1c0000000p-4, -0x1.8381420000000p-3, 0x1.0d239e0000000p-3, 0x1.d8fce80000000p-4,
    -0x1.857f960000000p-3, 0x1.dcf6f40000000p-5, -0x1.19ec2c0000000p-3, 0x1.a2cfea0000000p-5,
    -0x1.4723340000000p-2, 0x1.4b93480000000p-3, 0x1.e417d80000000p-5, -0x1.dd582c0000000p-5,
    0x1.4f0e300000000p-2, 0x1.c6d3220000000p-4, 0x1.2155b00000000p-2, -0x1.530f360000000p-2,
    -0x1.7f09ce0000000p-4, -0x1.5be80c0000000p-2, -0x1.9f39000000000p-5, 0x1.2814fc0000000p-2,
    -0x1.a84bd40000000p-2, 0x1.455b4a0000000p-2, -0x1.360ad80000000p-11, 0x1.c7fe800000000p-5,
    -0x1.7b4c3e0000000p-3, 0x1.d697c40000000p-3, -0x1.4759bc0000000p-2, -0x1.24bd3c0000000p-2,
    0x1.180dc60000000p-3, 0x1.80d9de0000000p-3, -0x1.b6b0400000000p-4, -0x1.059d620000000p-2,
    0x1.05948a0000000p-2, 0x1.1ad3060000000p-2, 0x1.2e812e0000000p-5, 0x1.366f980000000p-6,
    -0x1.6eca8a0000000p-2, -0x1.0b41e80000000p-2, 0x1.8514100000000p-2, -0x1.32d6220000000p-2,
};

inline constexpr std::array<double,8> b2 = {
    0x1.e5c8840000000p-4, 0x1.502b9e0000000p-9, -0x1.7eae0e0000000p-5, 0x1.a7f2580000000p-10,
    0x1.ae365c0000000p-5, 0x1.92cab00000000p-5, 0x1.36d21e0000000p-4, -0x1.98063a0000000p-6,
};

inline constexpr std::array<double,8> w3 = {
    0x1.baad540000000p-2, -0x1.6e4a840000000p-2, 0x1.16cc080000000p-2, -0x1.2336780000000p-2,
    0x1.10a6e00000000p-2, -0x1.a7b7ea0000000p-2, 0x1.e3cfca0000000p-2, -0x1.29b4ce0000000p-1,
};

inline constexpr std::array<double,1> b3 = {
    0x1.d589e20000000p-5,
};

inline RealExpression constant(double value) {
    return RealExpression(Real(ExactDouble(value)));
}

inline RealExpression tanh_expression(RealExpression const& value) {
    // No primitive symbolic tanh node exists today; this is the exact real
    // identity tanh(z)=(exp(2z)-1)/(exp(2z)+1).
    RealExpression e=exp(2*value);
    return (e-1)/(e+1);
}

inline RealExpression network(RealExpression const& x0, RealExpression const& x1) {
    std::array<RealExpression,width> h1;
    for(SizeType i=0u; i!=width; ++i) {
        RealExpression z=constant(b1[i])
            + constant(w1[2u*i])*x0
            + constant(w1[2u*i+1u])*x1;
        h1[i]=tanh_expression(z);
    }
    std::array<RealExpression,width> h2;
    for(SizeType i=0u; i!=width; ++i) {
        RealExpression z=constant(b2[i]);
        for(SizeType j=0u; j!=width; ++j) {
            z=z+constant(w2[width*i+j])*h1[j];
        }
        h2[i]=tanh_expression(z);
    }
    RealExpression output=constant(b3[0]);
    for(SizeType i=0u; i!=width; ++i) {
        output=output+constant(w3[i])*h2[i];
    }
    return output;
}

} // namespace Ariadne::TestBarr3Prefix8

#endif
