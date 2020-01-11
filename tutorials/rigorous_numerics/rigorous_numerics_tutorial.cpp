/***************************************************************************
 *            rigorous_numerics_tutorial.cpp
 *
 *  Copyright  20-17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <ariadne.hpp>

#define print(expr) { std::cout << #expr << ": " << (expr) << "\n"; }

using namespace Ariadne;

extern template Ariadne::Nat Ariadne::Error<Ariadne::FloatMP>::output_places;
extern template Ariadne::Nat Ariadne::Approximation<Ariadne::FloatMP>::output_places;

int main() {
    //! [numeric_demonstration]
    {
        auto r = 6 * atan(1/sqrt(Real(3)));
        print(r);
        FloatDPBounds xdp=r.get(double_precision);
        FloatMPBounds xmp=r.get(precision(128));
        ValidatedReal ymp=r.compute(Effort(128));
        FloatMPBall zmp=FloatMPBall(r.compute(Accuracy(128)).get(precision(128)));
        print(xdp); print(xmp); print(ymp); print(zmp);
        print(xmp.error()); print(xmp-xmp); print(xmp-xmp+1);
        print(zmp.error()); print(zmp-zmp); print(zmp-zmp+1);
        std::cout << std::endl;

        EffectiveNumber y=r;
        print(y);
    }
    //! [numeric_demonstration]


    //! [expression_demonstration]
    {
        RealVariable x("x");RealVariable y("y");
        RealConstant c("c",3.75_dy);
        RealExpression e = c * x * (1-x);
        print(x); print(c); print(e);
        Real x0=1/2_q; Real y0=-2/3_q;
        RealValuation v({x|x0,y|y0});
        print(v);
        auto x1=evaluate(e,v);
        print(x1); print(x1.get(double_precision))
        std::cout << std::endl;
    }
    //! [expression_demonstration]

    //! [linear_algebra_demonstration]
    {
        Matrix<FloatMPApproximation> A({{4,1,0},{1,4,1},{0,1,4}},precision(128));
        Vector<FloatMPApproximation> v({2.0,3,5},precision(128));
        print(inverse(A));
        FloatMPApproximation::set_output_places(30);
        print(inverse(A));
        print(solve(A,v));

        Matrix<FloatDPBounds> A_approx({{4,1,0},{1,4,1},{0,1,4}},double_precision);
        Vector<FloatDPBounds> v_approx({2.0,3,5},double_precision);
        print(A_approx);
        print(inverse(A_approx));
        print(solve(A_approx,v_approx));
        print(gs_solve(A_approx,v_approx));
        std::cout << std::endl;
    }
    //! [linear_algebra_demonstration]

    //! [function_demonstration]
    {
        Real a(1.875_dy), b(0.3_dec);
        auto id=EffectiveVectorMultivariateFunction::identity(EuclideanDomain(2));
        auto x=id[0]; auto y=id[1];
        auto h = EffectiveVectorMultivariateFunction{a-x*x-b*y,x};
//        Vector<FloatDPValue> v({0.5_x,1.0_x},double_precision);
          Vector<FloatDPValue> v({0.5,1.0},double_precision);
        print(v);
        print(h(v))
        print(evaluate(h,v));
        print(h.jacobian(v));
        print(jacobian(h,v));
        print(h.differential(v,3));
        std::cout << std::endl;

//        BoxDomainType dom({{0.0_x,1.0_x},{0.5_x,1.5_x}});
          BoxDomainType dom({{0,1},{0.5,1.5}});
        auto th = ValidatedVectorMultivariateTaylorFunctionModelDP(dom,h,ThresholdSweeper<FloatDP>(double_precision,1e-4));
        print(th);
        auto thh=compose(h,th);
        print(thh);
        print(evaluate(thh,v));
        print(h(h(v)));
        std::cout << std::endl;
    }
    //! [function_demonstration]

    //! [geometry_demonstration]
    {
        auto x=EffectiveScalarMultivariateFunction::coordinate(EuclideanDomain(2),0);
        auto y=EffectiveScalarMultivariateFunction::coordinate(EuclideanDomain(2),1);
        auto g = sqr(x)+4*sqr(y);
        auto h = EffectiveVectorMultivariateFunction{1+x+y*y,2+x-y};
        EffectiveConstraint c=(g<=1);
        ConstraintSet cs={c};
        RealBox bx={{-2,+2},{-2,+2}};
        ApproximateBoxType bbx={{-2,+4},{-2,+4}};
        BoundedConstraintSet bcs({{-2,+2},{-2,+2}},{g<=1});
        bcs=intersection(bx,cs);
        ConstrainedImageSet cis=image(bcs,h);
        print(cis);
        Figure fig(bbx,0,1); // Set up a figure for a space restricted to the bounding box, drawing first two coordinates
        fig << fill_colour(0.0,0.5,0.5) << bbx; // Use stream insertion to control graphics
        fig.set_fill_colour(0,1,1).draw(cis); // Use chained methods to control graphics
        fig.write("rigorous_numerics_tutorial"); //
        plot("rigorous_numerics_tutorial",PlanarProjectionMap(2,0,1),bbx,Colour(1,0.5,0.5),bbx,Colour(0,1,1),cis); // Plot in one command

    }
    //! [geometry_demonstration]

    return 0;
}

