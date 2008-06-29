/***************************************************************************
 *            integrator_base.h
 *
 *  Copyright  2006-8  Alberto Casagrande, Pieter Collins
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
 
#ifndef ARIADNE_INTEGRATOR_BASE_H
#define ARIADNE_INTEGRATOR_BASE_H

#include "integrator_interface.h"

namespace Ariadne {

template<class ES, class TM=Rational>
class IntegratorBase
  : public IntegratorInterface<ES,TM>
{
  typedef typename ES::real_type R;
  ES evolution_step(const VectorField<R>& f, const ES& s,
                    const TM& t1, const TM& t2, 
                    const Box<R>& bb) const;
};


       

template<class ES, class TM>
ES
IntegratorBase<ES,TM>::
evolution_step(const VectorField<R>& f, 
               const ES& s,
               const TM& t1, 
               const TM& t2, 
               const Box<R>& bb) const
{
  ARIADNE_ASSERT(t2>=t1);
  ES es=this->integration_step(f,s,t1,bb); 
  if(t1==t2) { return es; }
  else { return this->reachability_step(f,es,TM(t2-t1),bb); }
}
          


} // namespace Ariadne

#endif /* ARIADNE_INTEGRATOR_BASE_H */

