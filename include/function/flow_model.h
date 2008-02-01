/***************************************************************************
 *            flow_model.h
 *
 *  Copyright  2007 Pieter Collins
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
 
/*! \file flow_model.h
 *  \brief Models for flows. 
 */

#ifndef ARIADNE_FLOW_MODEL_H
#define ARIADNE_FLOW_MODEL_H

#include <iosfwd>
#include <string>
#include <sstream>

#include "base/types.h"
#include "base/array.h"
#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"

#include "function/taylor_series.h"
#include "function/taylor_variable.h"

namespace Ariadne {
  
  namespace Function {
  
    template<class R> 
    class FlowModel 
    {
      typedef Numeric::Interval<R> I;
     public:
      FlowModel(const array< TaylorSeries<TaylorVariable<I> > >& atstv) 
        : _variables(atstv) { }
      TaylorDerivative<I> evaluate(const I& t) const {
        TaylorDerivative<I> result(this->_variables.size(),this->_variables[0][0].argument_size(),this->_variables[0][0].degree());
        for(uint i=0; i!=this->_variables.size(); ++i) {
          result[i]=this->_variables[i][0]; 
          for(uint j=1; j<=_variables[i].degree(); ++j) {
            result[i]+=this->_variables[i][j]*pow(t,j); } }
        return result; }
      const array< TaylorSeries< TaylorVariable<I> > >& variables() const {
        return this->_variables; }
     private:
      const array< TaylorSeries< TaylorVariable<I> > > _variables;
    };


    template<class R>
    std::ostream& operator<<(std::ostream& os, const FlowModel<R>& fm) {
      return os << fm.variables(); 
    }

  }
}

#endif /* ARIADNE_FLOW_MODEL_H */
