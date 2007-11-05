/***************************************************************************
 *            parallelotope.inline.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 *  along with this program; if not, write to bouthe Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
namespace Ariadne {

    
    
template<class XC,class XG> inline 
Geometry::Parallelotope<XC,XG>::Parallelotope(dimension_type n)
  : Zonotope<XC,XG>(n,n) 
{ 
}


template<class XC,class XG> template<class RC, class RG> inline 
Geometry::Parallelotope<XC,XG>::Parallelotope(const Point<RC>& c, const LinearAlgebra::Matrix<RG>& G)
  : Zonotope<XC,XG>(c,G)
{
  if (G.number_of_rows()!=G.number_of_columns()) {
    ARIADNE_THROW(InvalidGenerators,"Parallelotope::Parallelotope(Vector c, Matrix G)"," G="<<G<<" is not a square matrix");
  }
}


template<class XC,class XG> template<class RR> inline 
Geometry::Parallelotope<XC,XG>::Parallelotope(const Rectangle<RR>& r)
  : Zonotope<XC,XG>(r) 
{ 
}


template<class XC, class XG> inline 
Geometry::Parallelotope<XC,XG>::Parallelotope(const std::string& str)
  : Zonotope<XC,XG>(str) 
{ 
  if (this->dimension()!=this->number_of_generators()) {
    ARIADNE_THROW(InvalidGenerators,"Parallelotope::Parallelotope(string str)"," str=\""<<str<<"\" does not describle a parallelotope");
  }
}


template<class XC, class XG> template<class RC, class RG> inline 
Geometry::Parallelotope<XC,XG>::Parallelotope(const Zonotope<RC,RG>& z)
  : Zonotope<XC,XG>(z) 
{ 
  if (z.dimension()!=z.number_of_generators()) {
    ARIADNE_THROW(InvalidGenerators,"Parallelotope::Parallelotope(Zonotope z)"," z="<<z<<" is not a parallelotope");
  }
}

template<class XC, class XG> template<class RC, class RG> inline 
Geometry::Parallelotope<XC,XG>::Parallelotope(const Parallelotope<RC,RG>& original)
  : Zonotope<XC,XG>(original) 
{ 
}


template<class XC, class XG> template<class RR> inline 
Geometry::Parallelotope<XC,XG>& 
Geometry::Parallelotope<XC,XG>::operator=(const Rectangle<RR>& r)
{
  Zonotope<XC,XG>& z=*this; z=r; return *this;
}




template<class XC, class XG> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const Parallelotope<XC,XG>& p) 
{
  return p.write(os);
}


template<class XC, class XG> inline
std::istream& 
Geometry::operator>>(std::ostream& is, Parallelotope<XC,XG>& p) 
{
  return p.read(is);
}




} // namespace Ariadne
