/***************************************************************************
 *            polytope.inline.h
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it   Pieter.Collins@cwi.nl
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
  namespace Geometry {

  
    template<class R> template<class Rl> inline
    Polytope<R>::Polytope(const Polyhedron<Rl>& p)
      : _dimension(p.dimension()), _number_of_vertices(), _data()
    {   
      (*this)=Geometry::polytope(Polyhedron<R>(p));
    }
    
    template<class R> template<class Rl> inline
    Polytope<R>::Polytope(const Polytope<Rl>& p)
      : _dimension(p.dimension()), _number_of_vertices(p.number_of_vertices()), _data(p.data())
    { 
    }
    

    
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const Polytope<R>& p)
    {
      return p.write(os); 
    }
    
    
    template<class R> inline
    std::istream& operator>>(std::istream& os, Polytope<R>& p) 
    {
      return p.read(os); 
    }
    
    
    
    
    template<class R>
    class PolytopeVerticesIterator
      : public boost::iterator_facade<PolytopeVerticesIterator<R>,
                                      Point<R>,
                                      boost::forward_traversal_tag,
                                      Point<R> const&,
                                      Point<R> const*
                                     >
    {
     public:
      PolytopeVerticesIterator(const Polytope<R>& plytp, const size_type& j) 
        : _p(&plytp), _j(j) { }
      bool equal(const PolytopeVerticesIterator<R>& other) const {
        return this->_j==other._j && this->_p ==other._p; }
      void increment() {
        ++this->_j; }
      const Point<R>& dereference() const {
        this->_v=this->_p->vertex(this->_j); return this->_v; }
     private:
      const Polytope<R>* _p; size_type _j; mutable Point<R> _v;
    };
    


   
    

  }
}
