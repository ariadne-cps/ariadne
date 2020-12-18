/***************************************************************************
 *            symbolic_usage.cpp
 *
 *  Copyright  2020  Pieter Collins
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

#include "ariadne.hpp"

using namespace Ariadne;


int main() {
    {
        //! [Constant_usage]
        RealConstant g("g",9.81_dec);
        //! [Constant_usage]
    }

    {
        //! [Variable_usage]
        IntegerVariable z("z");
        RealVariables x("x",3); // Defines variables x[i] with name xi for i=0,1,2
        //! [Variable_usage]
    }

    {
        //! [Expression_usage]
        RealVariable x("x"), v("v"); TimeVariable t;
        RealExpression a=-v+x*(1-x*x)+cos(t);
        Real a_val=evaluate(a,{t|0,x|1,v|0});
        //! [Expression_usage]
    }

    {
        //! [PrimedVariable_usage]
        RealVariable x("x"), v("v"); Real dt(0.1_dec);
        next(x)=x+v*dt;
        //! [PrimedVariable_usage]
    }

    {
        //! [DottedVariable_usage]
        RealVariable x("x"), v("v");
        dot(v)=-v+sin(x); // Define the differential equation dv/dt=-v+sin(x)
        //! [DottedVariable_usage]
    }

    {
        //! [Valuation_usage]
        RealVariable x("x"), y("y"); TimeVariable t;
        RealValuation val({x|2,y|3.5_dy,t|1/3_q});
        Real x_val=val[x];
        //! [Valuation_usage]
    }

    {
        //! [Space_usage]
        RealVariable x("x"), y("y"); TimeVariable t;
        RealSpace spc({x,y,t});
        RealValuation val({t|1/3_q,x|2,y|3.5_dy}); // Define values for the variables in arbitrary order
        RealVector vec=val[spc]; // Extract the values in the order \a [x,y,t]
        val=RealValuation(vec,spc); // Recover the original valuation
        //! [Space_usage]
    }
}
