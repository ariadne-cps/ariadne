/***************************************************************************
 *            floatmp.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

#include "numeric/floatmp.h"
#include "numeric/interval.h"

namespace Ariadne {
  namespace Numeric {

    Interval<FloatMP> 
    operator+(const FloatMP& x1, const FloatMP& x2) 
    {
      return Interval<FloatMP>(add_down(x1,x2),add_up(x1,x2));
    }
    
    Interval<FloatMP> 
    operator-(const FloatMP& x1, const FloatMP& x2) 
    {
      return Interval<FloatMP>(sub_down(x1,x2),sub_up(x1,x2));
    }

    Interval<FloatMP> 
    operator*(const FloatMP& x1, const FloatMP& x2) 
    {
      return Interval<FloatMP>(mul_down(x1,x2),mul_up(x1,x2));
    }

    Interval<FloatMP> 
    operator/(const FloatMP& x1, const FloatMP& x2) 
    {
      return Interval<FloatMP>(div_down(x1,x2),div_up(x1,x2));
    }

  } 
}
