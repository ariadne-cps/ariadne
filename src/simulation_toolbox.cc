/***************************************************************************
 *            simulation_toolbox.cc
 *
 *  Copyright  2008  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <iomanip>

#include "macros.h"
#include "logging.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "function.h"
#include "point.h"

#include "simulation_toolbox.h"


namespace Ariadne {

class NoCrossingException { };
class DegenerateCrossingException { };
class NonInvertibleFunctionException { };






SimulationToolbox::
SimulationToolbox()
{
}


tribool
SimulationToolbox::
active(const ScalarFunction& guard,
       const Point& point) const
{
    Float value=guard.evaluate(point);
    if(value<0) { return false; }
    else if(value>0) { return true; }
    else { return indeterminate; }
}

Point
SimulationToolbox::
reset_step(const VectorFunction& map,
           const Point& point) const
{
    return map.evaluate(point);
}



Point
SimulationToolbox::
integration_step(const VectorFunction& f,
                 const Point& pt,
                 const TimeType& h) const
{
    // Use 4th-Order Runge-Kutta method as default.
    Vector<RealType> v=pt;
    Vector<RealType> k1=f.evaluate(v);
    v+=h*k1;
    Vector<RealType> k2=f.evaluate(v);
    v=pt; v+=(h/2)*k2;
    Vector<RealType> k3=f.evaluate(v);
    v=pt; v+=(h/2)*k3;
    Vector<RealType> k4=f.evaluate(v);

    return pt+(h/6)*(k1+2*k2+2*k3+k4);
}




// Compute the crossing time using bisections
SimulationToolbox::TimeType
SimulationToolbox::
crossing_time(const ScalarFunction& g,
              const VectorFunction& f,
              const Point& pt,
              const TimeType& h) const
{
    if(g.evaluate(pt)>=0.0) { return 0.0; }

    if(g.evaluate(Vector<RealType>(pt+h*f.evaluate(pt)))<0.0) { throw NoCrossingException(); }
    if(g.evaluate(this->integration_step(f,pt,h))<0.0) { throw NoCrossingException(); }


    // Use bisection to find crossing
    TimeType h0=0.0; TimeType h1=h;
    while(h0!=h1) {
        TimeType hc=(h0+h1)/2;
        if(g.evaluate(this->integration_step(f,pt,hc))<0.0) {
            h0=hc;
        } else {
            h1=hc;
        }
    }

    return h1;

}

}  // namespace Ariadne
