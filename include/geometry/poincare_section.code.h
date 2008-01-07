/***************************************************************************
 *            poincare_section.code.h
 *
 *  Copyright  2007  Pieter Collins
 *  pieter.collins@cwi.nl
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
 
#include "function/function_interface.h"


namespace Ariadne {
  
template<class R> inline
Geometry::PoincareSection<R>::~PoincareSection() 
{
  delete _inclusion_map_ptr;
  delete _projection_map_ptr;
  delete _constraint_ptr;
}

template<class R> inline
Geometry::PoincareSection<R>::PoincareSection(const Function::FunctionInterface<R>& im, 
                                              const Function::FunctionInterface<R>& pm, 
                                              const Function::FunctionInterface<R>& c) 
  :  _inclusion_map_ptr(im.clone()),
     _projection_map_ptr(pm.clone()),
     _constraint_ptr(c.clone())
{
  assert(c.result_size()==1u);
  assert(im.argument_size()+1u==c.argument_size());
  assert(im.result_size()==c.argument_size());
  assert(pm.argument_size()==c.argument_size());
  assert(pm.result_size()+1u==c.argument_size());
}

template<class R> inline
Geometry::PoincareSection<R>::PoincareSection(const Geometry::PoincareSection<R>& sec)
  :  _inclusion_map_ptr(sec._inclusion_map_ptr->clone()),
     _projection_map_ptr(sec._projection_map_ptr->clone()),
     _constraint_ptr(sec._constraint_ptr->clone())
{
}

template<class R> inline
Geometry::PoincareSection<R>*
Geometry::PoincareSection<R>::clone() const
{
  return new PoincareSection<R>(*this);
}

template<class R> inline
dimension_type
Geometry::PoincareSection<R>::dimension() const
{
  return this->_inclusion_map_ptr->argument_size();
}

template<class R> inline
smoothness_type
Geometry::PoincareSection<R>::smoothness() const
{
  return std::min(std::min(this->_projection_map_ptr->smoothness(),this->_inclusion_map_ptr->smoothness()),
                  this->_constraint_ptr->smoothness());
}

template<class R> inline
const Function::FunctionInterface<R>&
Geometry::PoincareSection<R>::inclusion_map() const
{
  return *this->_inclusion_map_ptr;
}

template<class R> inline
const Function::FunctionInterface<R>&
Geometry::PoincareSection<R>::projection_map() const
{
  return *this->_projection_map_ptr;
}

template<class R> inline
const Function::FunctionInterface<R>&
Geometry::PoincareSection<R>::crossing_condition() const
{
  return *this->_constraint_ptr;
}

template<class R> inline
std::ostream&
Geometry::PoincareSection<R>::write(std::ostream& os) const
{
  return os << "PoincareSection( inclusion_map=" << this->inclusion_map()
            << ", projection_map=" << this->projection_map() 
            << ", constraint=" << this->crossing_condition() << " )";
}

}
