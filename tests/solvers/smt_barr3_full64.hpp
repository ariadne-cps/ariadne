/***************************************************************************
 * Published Barr3 2-64-64-1 tanh network fixture.
 *
 * Original checkpoint SHA-256:
 * 788fd56f21abb0cb5b83206b12d4d8d0ed718896c2140222a8fd129475f6c213
 *
 * Internal raw tensor payload SHA-256:
 * d1c6646a9354d44092e23de07495c40d2a6575f2235aaa5be0d68b1e9b61bdbf
 *
 * Payload layout, all IEEE-754 binary32, little-endian:
 * W1[64][2], b1[64], W2[64][64], b2[64], W3[1][64], b3[1].
 ***************************************************************************/
#ifndef ARIADNE_TEST_SMT_BARR3_FULL64_HPP
#define ARIADNE_TEST_SMT_BARR3_FULL64_HPP

#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <vector>

namespace Ariadne::TestBarr3Full64 {

inline constexpr SizeType width=64u;
inline constexpr SizeType parameter_count=4417u;
inline constexpr SizeType w1_offset=0u;
inline constexpr SizeType b1_offset=128u;
inline constexpr SizeType w2_offset=192u;
inline constexpr SizeType b2_offset=4288u;
inline constexpr SizeType w3_offset=4352u;
inline constexpr SizeType b3_offset=4416u;

inline std::array<double,parameter_count> const& parameters() {
    static const std::array<double,parameter_count> values=[] {
        std::ifstream stream(ARIADNE_SMT_BARR3_DATA_PATH,std::ios::binary);
        if(not stream) {
            throw std::runtime_error("Cannot open Barr3 64x64 tensor payload");
        }
        std::vector<unsigned char> bytes(
            (std::istreambuf_iterator<char>(stream)),
            std::istreambuf_iterator<char>());
        if(bytes.size()!=parameter_count*4u) {
            throw std::runtime_error("Invalid Barr3 64x64 tensor payload size");
        }
        std::array<double,parameter_count> decoded{};
        for(SizeType i=0u;i!=parameter_count;++i) {
            SizeType const offset=4u*i;
            std::uint32_t bits=
                static_cast<std::uint32_t>(bytes[offset])
                | (static_cast<std::uint32_t>(bytes[offset+1u])<<8u)
                | (static_cast<std::uint32_t>(bytes[offset+2u])<<16u)
                | (static_cast<std::uint32_t>(bytes[offset+3u])<<24u);
            decoded[i]=static_cast<double>(std::bit_cast<float>(bits));
        }
        return decoded;
    }();
    return values;
}

inline RealExpression constant(double value) {
    return RealExpression(Real(ExactDouble(value)));
}

inline RealExpression tanh_expression(RealExpression const& value) {
    return tanh(value);
}

inline RealExpression network(RealExpression const& x0, RealExpression const& x1) {
    auto const& p=parameters();
    std::array<RealExpression,width> h1;
    for(SizeType i=0u;i!=width;++i) {
        RealExpression z=constant(p[b1_offset+i])
            + constant(p[w1_offset+2u*i])*x0
            + constant(p[w1_offset+2u*i+1u])*x1;
        h1[i]=tanh_expression(z);
    }
    std::array<RealExpression,width> h2;
    for(SizeType i=0u;i!=width;++i) {
        RealExpression z=constant(p[b2_offset+i]);
        for(SizeType j=0u;j!=width;++j) {
            z=z+constant(p[w2_offset+width*i+j])*h1[j];
        }
        h2[i]=tanh_expression(z);
    }
    RealExpression output=constant(p[b3_offset]);
    for(SizeType i=0u;i!=width;++i) {
        output=output+constant(p[w3_offset+i])*h2[i];
    }
    return output;
}


struct NetworkAndLie {
    RealExpression barrier;
    RealExpression lie;
    RealExpression db_dx;
    RealExpression db_dy;
};

inline NetworkAndLie network_and_lie(
    RealExpression const& x0,
    RealExpression const& x1)
{
    auto const& p=parameters();

    std::array<RealExpression,width> h1;
    std::array<RealExpression,width> dh1_dx;
    std::array<RealExpression,width> dh1_dy;
    for(SizeType i=0u;i!=width;++i) {
        RealExpression z=constant(p[b1_offset+i])
            + constant(p[w1_offset+2u*i])*x0
            + constant(p[w1_offset+2u*i+1u])*x1;
        h1[i]=tanh_expression(z);
        RealExpression factor=1-sqr(h1[i]);
        dh1_dx[i]=factor*constant(p[w1_offset+2u*i]);
        dh1_dy[i]=factor*constant(p[w1_offset+2u*i+1u]);
    }

    std::array<RealExpression,width> h2;
    std::array<RealExpression,width> dh2_dx;
    std::array<RealExpression,width> dh2_dy;
    for(SizeType i=0u;i!=width;++i) {
        RealExpression z=constant(p[b2_offset+i]);
        RealExpression dz_dx=RealExpression(0);
        RealExpression dz_dy=RealExpression(0);
        for(SizeType j=0u;j!=width;++j) {
            RealExpression weight=constant(p[w2_offset+width*i+j]);
            z=z+weight*h1[j];
            dz_dx=dz_dx+weight*dh1_dx[j];
            dz_dy=dz_dy+weight*dh1_dy[j];
        }
        h2[i]=tanh_expression(z);
        RealExpression factor=1-sqr(h2[i]);
        dh2_dx[i]=factor*dz_dx;
        dh2_dy[i]=factor*dz_dy;
    }

    RealExpression barrier=constant(p[b3_offset]);
    RealExpression db_dx=RealExpression(0);
    RealExpression db_dy=RealExpression(0);
    for(SizeType i=0u;i!=width;++i) {
        RealExpression weight=constant(p[w3_offset+i]);
        barrier=barrier+weight*h2[i];
        db_dx=db_dx+weight*dh2_dx[i];
        db_dy=db_dy+weight*dh2_dy[i];
    }

    // Published Barr3 dynamics:
    //   x_dot = y
    //   y_dot = -x-y+x^3/3
    //
    // Use an algebraically equivalent factored form for y_dot so validated
    // interval evaluation preserves the self-correlation in x^2 through sqr.
    RealExpression dx=x1;
    RealExpression dy=x0*(sqr(x0)/3-1)-x1;
    RealExpression lie=db_dx*dx+db_dy*dy;
    return {barrier,lie,db_dx,db_dy};
}

} // namespace Ariadne::TestBarr3Full64

#endif
