/***************************************************************************
 *            zonotope.template.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
 

namespace Ariadne {
    
namespace Geometry { template<class BS> ListSet<BS> subdivide(const BS&, const typename BS::real_type&); }

template<class R> template<class XX> 
Geometry::Zonotope<R>::Zonotope(dimension_type d, size_type m, const XX* ptr)
  : _centre(d,ptr), _generators(d,m,ptr+d), _error(d)
{
}

template<class R> template<class XX> 
Geometry::Zonotope<R>::Zonotope(const Zonotope<XX>& z)
  : _centre(z.centre()), _generators(z.generators()), _error(z.error())
{
}

template<class BS>
Geometry::ListSet<BS> 
Geometry::subdivide(const BS& bs, const typename BS::real_type& r)
{
  typedef typename BS::real_type R;
  ListSet<BS> result;
  ListSet<BS> working(bs);
  while(!working.size()==0) {
    BS set=working.pop();
    if(set.radius()<r) {
      result.adjoin(set);
    } else {
      working.adjoin(split(set));
    }
  }
  return result;
}





} // namespace Ariadne

