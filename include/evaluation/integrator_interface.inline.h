/***************************************************************************
 *            integrator_interface.inline.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
#include "integrator_interface.h"

namespace Ariadne {

template<class BS>
BS
Evaluation::IntegratorInterface<BS>::
evolution_step(const System::VectorField<R>& f, 
               const BS& s,
               const Numeric::Rational& t1, 
               const Numeric::Rational& t2, 
               const Geometry::Box<R>& bb) const
{
  ARIADNE_ASSERT(t2>=t1);
  BS es=this->integration_step(f,s,t1,bb); 
  if(t1==t2) { return es; }
  else { return this->reachability_step(f,es,Numeric::Rational(t2-t1),bb); }
}
          

template<class R> inline
std::ostream&
Evaluation::operator<<(std::ostream& os, const IntegratorInterface<R>& i) 
{
  return i.write(os);
}

} // namespace Ariadne
