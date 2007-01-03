/***************************************************************************
 *            grid.inline.h
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
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
 *  Foundation, Inc., 59 Templece Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 

namespace Ariadne {
  namespace Geometry {

    template<class R> inline
    Grid<R>::~Grid() 
    {
    }
    
    
    template<class R> inline
    index_type Grid<R>::subdivision_index(dimension_type d, const real_type& x) const 
    {
      index_type n=subdivision_interval(d,x);
      if(subdivision_coordinate(d,n) == x) { return n; }
      throw std::runtime_error("Value is not a grid coordinate");
    }
    
    
    template<class R> inline
    index_type 
    Grid<R>::subdivision_lower_index(dimension_type d, const real_type& x) const
    {
      return subdivision_interval(d,x);
    }
    
    
    template<class R> inline
    index_type 
    Grid<R>::subdivision_upper_index(dimension_type d, const real_type& x) const 
    {
      index_type n=subdivision_interval(d,x);
      return subdivision_coordinate(d,n) == x ? n : n+1;
    }



   
    template<class R> inline
    grid_type
    IrregularGrid<R>::type() const 
    {
      return IRREGULAR;
    }
    
    
    template<class R> inline
    dimension_type 
    IrregularGrid<R>::dimension() const 
    { 
      return this->_subdivision_coordinates.size();
    }
    
    
    template<class R> inline
    R 
    IrregularGrid<R>::subdivision_coordinate(dimension_type d, index_type n) const 
    {
      check_coordinate(*this,d,__PRETTY_FUNCTION__);
      if(!(0<=n && uint(n)<_subdivision_coordinates[d].size())) {
        std::cerr << "d=" << d << ", n=" << n << ", size=" << _subdivision_coordinates[d].size() << std::endl;
        throw std::runtime_error("index does not lie in range of finite grid");
      }
      return _subdivision_coordinates[d][n];
    }
    
    
    template<class R> inline
    index_type 
    IrregularGrid<R>::subdivision_interval(dimension_type d, const real_type& x) const 
    {
      typename std::vector<R>::const_iterator pos;
      check_coordinate(*this,d,__PRETTY_FUNCTION__);
      if(x<_subdivision_coordinates[d].front() || x>_subdivision_coordinates[d].back()) {
        std::cerr << "d.front()=" << _subdivision_coordinates[d].front()
                  <<  " d.back()=" << _subdivision_coordinates[d].back()
                  <<  " x="<< x <<std::endl<<std::flush;
        throw std::runtime_error("point does not lie in extent of finite grid");
      }
      pos = std::upper_bound(_subdivision_coordinates[d].begin(),
                             _subdivision_coordinates[d].end(), x);
      return (pos - _subdivision_coordinates[d].begin()) - 1;
    }

      /*! \brief Tests whether the grid contains the given lattice rectangle 
       * within its bounds. */
    template<class R> inline 
    bool 
    IrregularGrid<R>::bounds_enclose(const Rectangle<R>& r) const 
    {
      return subset(r,Rectangle<R>(this->bounds())); 
    }
            
      
    template<class R> inline 
    IndexArray 
    IrregularGrid<R>::lower() const 
    { 
      return IndexArray(dimension(),0);
    }

    template<class R> inline 
    IndexArray 
    IrregularGrid<R>::upper() const 
    { 
      IndexArray result(this->dimension());
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        result[i]=_subdivision_coordinates[i].size()-1;
      }
      return result;
    }
      
    template<class R> inline 
    Combinatoric::LatticeBlock 
    IrregularGrid<R>::block() const 
    { 
      return Combinatoric::LatticeBlock(lower(),upper()); 
    }

    template<class R> inline 
    SizeArray 
    IrregularGrid<R>::sizes() const 
    { 
      return block().sizes(); 
    }

    template<class R> inline 
    size_type 
    IrregularGrid<R>::capacity() const 
    { 
      return block().size(); 
    }

    template<class R> inline 
    size_type 
    IrregularGrid<R>::size(dimension_type i) const 
    { 
      return _subdivision_coordinates[i].size()-1; 
    }
    
    
    
    
    
  template<class R> inline 
  grid_type 
  RegularGrid<R>::type() const 
  {
    return REGULAR;
  }
  
    
  template<class R> inline
  RegularGrid<R>::RegularGrid(const array<R>& sl)
    : _subdivision_lengths(sl) 
  {
  }

  template<class R> inline
  RegularGrid<R>::RegularGrid(const dimension_type& n, const R& l)
    : _subdivision_lengths(n,l) 
  {
  }


  template<class R> inline
  dimension_type 
  RegularGrid<R>::dimension() const 
  { 
    return _subdivision_lengths.size(); 
  }


  template<class R> inline
  R
  RegularGrid<R>::subdivision_coordinate(dimension_type d, index_type n) const
  { 
    return mul_approx(_subdivision_lengths[d] , n); 
  }


  template<class R> inline
  R
  RegularGrid<R>::subdivision_length(dimension_type d) const 
  { 
    return _subdivision_lengths[d]; 
  }

  template<class R> inline
  index_type
  RegularGrid<R>::subdivision_interval(dimension_type d, const real_type& x) const 
  {
    index_type result = int_down<index_type>(div_down(x,_subdivision_lengths[d]));
    return result;
  }

  template<class R> inline
  bool 
  RegularGrid<R>::bounds_enclose(const Rectangle<R>& r) const 
  {
    return true; 
  }
     
  template<class R> inline
  bool
  RegularGrid<R>::bounds_superset(const Rectangle<R>& r) const
  { 
    return true; 
  }





     
  template<class R> inline
  FiniteGrid<R>::FiniteGrid(const Grid<R>& g, const Combinatoric::LatticeBlock& b) 
    : _grid_ptr(&g), _grid_type(g.type()), _bounds(b)
  { 
  }
      
  template<class R> inline
  FiniteGrid<R>::FiniteGrid(const Grid<R>& g, const Rectangle<R>& bb) 
    : _grid_ptr(&g), _grid_type(g.type()), 
      _bounds(over_approximation(bb,g).lattice_set())
  { 
  }
      
  template<class R> inline
  const Grid<R>& 
  FiniteGrid<R>::grid() const
  { 
    return *_grid_ptr; 
  }
      
  template<class R> inline
  const Combinatoric::LatticeBlock& 
  FiniteGrid<R>::bounds() const 
  {
    return _bounds; 
  }
     

  template<class R> inline
  R 
  FiniteGrid<R>::subdivision_coordinate(dimension_type d, index_type n) const
  {
    return _grid_ptr->subdivision_coordinate(d,n);
  }
  
  template<class R> inline
  index_type 
  FiniteGrid<R>::subdivision_interval(dimension_type d, const real_type& x) const 
  {
    return _grid_ptr->subdivision_interval(d,x);
  }
      
  template<class R> inline
  index_type 
  FiniteGrid<R>::subdivision_lower_index(dimension_type d, const real_type& x) const 
  {
    return _grid_ptr->subdivision_lower_index(d,x);
  }

  template<class R> inline
  index_type 
  FiniteGrid<R>::subdivision_upper_index(dimension_type d, const real_type& x) const 
  {
    return _grid_ptr->subdivision_upper_index(d,x);
  }





    template<class R> inline
    std::istream& operator>>(std::istream& is, Grid<R>& g)
    {
      return g.read(is);
    }

    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const Grid<R>& g)
    {
      return g.write(os);
    }

    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const FiniteGrid<R>& fg)
    {
      return fg.write(os);
    }

  }
}
