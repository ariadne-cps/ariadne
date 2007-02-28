/***************************************************************************
 *            point_list.inline.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it   Pieter.Collins@cwi.nl
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

    template<class R>
    class PointListIterator
    {
     public:
      PointListIterator(const LinearAlgebra::Matrix<R>& A, size_type j) : _mx(A), _j(j) { }
      Point<R> operator*() const { return Point<R>(column(_mx,_j)); }
      void operator++() { ++_j; }
      bool operator==(const PointListIterator& other) { 
        return (this->_j==other._j) && (&(this->_mx)==&(other._mx)); }
      bool operator!=(const PointListIterator& other) { 
        return !(*this==other); }
     private:
      const LinearAlgebra::Matrix<R>& _mx;
      size_type _j;
    };
    
    template<class R>
    class PointListReference
    {
     public:
      PointListReference(PointList<R>& pl, size_type j) : _pl(pl), _j(j) { }
      void operator=(const Point<R>& pt) {
        for(dimension_type i=0; i!=_pl.dimension(); ++i) { _pl._pts(i,_j)=pt[i]; } }
      R& operator[] (size_type i) { return _pl._pts(i,_j); }
      operator Point<R>() const { return _pl.get(_j); }
     private:
      PointList<R>& _pl;
      size_type _j;
    };
    
    template<class R> inline
    PointList<R>::PointList(dimension_type d)
      : _size(0), _pts(d+1,1) 
    { 
    }

    template<class R> inline
    PointList<R>::PointList(dimension_type d, size_type n)
      : _size(n), _pts(d+1,n) 
    { 
    }

    template<class R> inline
    PointList<R>::PointList(const LinearAlgebra::Matrix<R>& G)
      : _size(G.number_of_columns()), _pts(G) 
    { 
    }
      
    template<class R> template<class Rl> inline
    PointList<R>::PointList(const PointList<Rl>& ptl)
      : _size(ptl.size()), _pts(ptl.generators()) 
    { 
    }
      
    template<class R> inline
    PointList<R>::PointList(const PointList<R>& ptl)
      : _size(ptl._size), _pts(ptl.generators()) 
    { 
    }

    template<class R> inline
    PointList<R>& 
    PointList<R>::operator=(const PointList<R>& ptl)
    {
      if(this!=&ptl) { 
        this->_size=ptl.size(); 
        this->_pts=ptl.generators(); 
      } 
      return *this; 
    }
      
    template<class R> inline
    dimension_type 
    PointList<R>::dimension() const
    {
      return _pts.number_of_rows()-1; 
    }

    template<class R> inline
    size_type 
    PointList<R>::size() const
    {
      return _size; 
    }

    template<class R> inline
    size_type 
    PointList<R>::capacity() const
    {
      return _pts.number_of_columns(); 
    }


    template<class R> inline
    Point<R> 
    PointList<R>::operator[] (const size_type j) const
    {
      Point<R> result(this->dimension());
      for(size_type i=0; i!=this->dimension(); ++i) {
        result[i]=this->_pts(i,j);
      }
      return result;
    }
    
    template<class R> inline
    PointListReference<R>  
    PointList<R>::operator[](const size_type j)
    { 
      return PointListReference<R>(*this,j);
    }


    template<class R> inline
    void 
    PointList<R>::pop_back()
    {
      --this->_size; 
    }
      
    template<class R> inline
    void 
    PointList<R>::adjoin(const Point<R>& pt)
    {
      this->push_back(pt); 
    }
    template<class R> inline
    void 
    PointList<R>::clear()
    {
      this->_size=0; 
    }

      
    template<class R> inline
    typename PointList<R>::const_iterator 
    PointList<R>::begin() const
    {
      return const_iterator(this->_pts,0); 
    }

    template<class R> inline
    typename PointList<R>::const_iterator 
    PointList<R>::end() const
    {
      return const_iterator(this->_pts,this->_size); 
    }
      


    template<class R> inline
    void 
    PointList<R>::_set(size_type j, dimension_type i,const R& x)
    {
      _pts(i,j)=x; 
    }
      
    template<class R> inline
    Point<R> 
    PointList<R>::_get(size_type j) const
    {
      return Point<R>(_pts.column(j)); 
    }

    template<class R> inline
    R 
    PointList<R>::_get(size_type j, dimension_type i) const 
    {
      return _pts(i,j); 
    }
    


    template<class R> inline
    std::ostream& 
    operator<<(std::ostream& os, const PointList<R>& pts)
    {
      return pts.write(os);
    }
    
    
  }
}
