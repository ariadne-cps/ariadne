/***************************************************************************
 *            grid.code.h
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

#include "grid.h"

#include <ostream>
#include <algorithm>

#include "../exceptions.h"
#include "../base/stlio.h"

#include "../numeric/arithmetic.h"

#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"

#include "../geometry/grid_set.h" // for GridBlock<R>

namespace Ariadne {
  namespace Geometry {


    template<class R> 
    Grid<R>::~Grid()
    {
    }
  
  

    template<class R> 
    Grid<R>::Grid(const dimension_type& d, const R& l)
      : _subdivision_coordinates(d,array<R>(1,static_cast<R>(0))),
        _subdivision_lengths(d,l),
        _centre_positions(d,0)
    {
      this->create();
    }
  

    template<class R> 
    Grid<R>::Grid(const LinearAlgebra::Vector<R>& v)
      : _subdivision_coordinates(v.size(),array<R>(1,static_cast<R>(0))),
        _subdivision_lengths(v.data()),
        _centre_positions(v.size(),0)
    {
      this->create();
    }
  

    template<class R> 
    Grid<R>::Grid(const Point<R>& pt, const LinearAlgebra::Vector<R>& v)
      : _subdivision_coordinates(v.size(),array<R>(1)),
        _subdivision_lengths(v.data()),
        _centre_positions(v.size(),0)
    {
      if(pt.dimension() != v.size()) {
        throw IncompatibleDimensions(__PRETTY_FUNCTION__);
      }
      for(size_type d=0; d!=pt.dimension(); ++d) {
        this->_subdivision_coordinates[d][0]=v[d];
      }
      this->create();
    }
  


    template<class R>
    Grid<R>::Grid(const array< array<R> >& sc, const array<R>& sl, const array<size_type>& cp) 
      : _subdivision_coordinates(sc),
        _subdivision_lengths(sl),
        _centre_positions(cp)
    {
      this->create();
    }

     
    template<class R>
    Grid<R>::Grid(const dimension_type& d, const size_type* nsc, const R** sc, const R* sl, const size_type* cp) 
      : _subdivision_coordinates(d),
        _subdivision_lengths(sl,sl+d),
        _centre_positions(cp,cp+d)
    {
      for(size_type i=0; i!=d; ++i) {
        this->_subdivision_coordinates[i]=array<R>(sc[i],sc[i]+nsc[i]);
      }
      this->create();
    }

     

    template<class R>
    Grid<R>::Grid(const ListSet< Rectangle<R> >& ls)
      : _subdivision_coordinates(ls.dimension()),
        _subdivision_lengths(ls.dimension(),static_cast<R>(1)),
        _centre_positions(ls.dimension(),0u)
    {
      dimension_type d=ls.dimension();
      R one=static_cast<R>(1);

      std::vector<R>* rectangle_bounds=new std::vector<R>[d];

      for(typename ListSet< Rectangle<R> >::const_iterator riter=ls.begin(); riter!=ls.end(); ++riter) {
        for(dimension_type i=0; i!=d; ++i) {
          rectangle_bounds[i].push_back(riter->lower_bound(i));
          rectangle_bounds[i].push_back(riter->upper_bound(i));
        }
      }

      for(dimension_type i=0; i!=d; ++i) {
        std::vector<R>& bounds=rectangle_bounds[i];
        bounds.push_back(static_cast<R>(0));
        std::sort(bounds.begin(),bounds.end());
        typename std::vector<R>::const_iterator bounds_end=std::unique(bounds.begin(),bounds.end());
        bounds.resize(bounds_end-bounds.begin());
        R bounds_floor=floor(bounds.front());
        R bounds_ceil=ceil(bounds.back());
        if(bounds_floor==bounds.front()) {
          bounds_floor=floor(sub_up(bounds_floor,one));
        }
        if(bounds_ceil==bounds.back()) {
          bounds_ceil=ceil(add_down(bounds_ceil,one));
        }

        array<R>& coordinates=this->_subdivision_coordinates[i];
        coordinates.resize(bounds.size()+2);
        coordinates[0]=bounds_floor;
        std::copy(bounds.begin(),bounds.end(),coordinates.begin()+1);
        coordinates[bounds.size()+1u]=bounds_ceil;
        
        this->_centre_positions[i] = std::find(coordinates.begin(),coordinates.end(),static_cast<R>(0))-coordinates.begin();
      }
      delete[] rectangle_bounds;
      this->create();
    }

    template<class R>
    Grid<R>::Grid(const Grid<R>& g1, Grid<R>& g2)
      : _subdivision_coordinates(g1.dimension()),
        _subdivision_lengths(g1._subdivision_lengths),
        _centre_positions(g1.dimension(),0)
    {
      if(g1._subdivision_lengths!=g2._subdivision_lengths) {
        throw std::runtime_error("Cannot merge two grids with different subdivision lengths");
      }
      for(dimension_type d=0; d!=this->dimension(); ++d) {
        array<R>& sc(this->_subdivision_coordinates[d]);
        const array<R>& sc1(g1._subdivision_coordinates[d]);
        const array<R>& sc2(g2._subdivision_coordinates[d]);
        sc.resize(sc1.size()+sc2.size());
        std::merge(sc1.begin(),sc1.end(),sc2.begin(),sc2.end(),sc.begin());
        index_type n=std::unique(sc.begin(),sc.end())-sc.begin();
        sc.resize(n);
        R x=sc1[g1._centre_positions[d]];
        this->_centre_positions[d]=std::find(sc.begin(),sc.end(),x)-sc.begin();
      }

      this->create();
    }



    template<class R>
    void
    Grid<R>::create() 
    {
      dimension_type d=this->dimension();
      this->_centre_pointers.resize(d);
      this->_lower_origin.resize(d);
      this->_upper_origin.resize(d);
      this->_lower_bounds.resize(d);
      this->_upper_bounds.resize(d);
      this->_lower_indices.resize(d);
      this->_upper_indices.resize(d);
      
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        this->_centre_pointers[i]=this->_subdivision_coordinates[i].begin()+_centre_positions[i];
        this->_lower_indices[i]=0u-this->_centre_positions[i];
        this->_upper_indices[i]=this->_subdivision_coordinates[i].size()-this->_centre_positions[i]-1;
        this->_lower_bounds[i]=this->_subdivision_coordinates[i][0];
        this->_upper_bounds[i]=this->_subdivision_coordinates[i][this->_subdivision_coordinates[i].size()-1];
        this->_lower_origin[i]=sub_approx(this->_lower_bounds[i],mul_approx(this->_subdivision_lengths[i],this->_lower_indices[i]));
        this->_upper_origin[i]=sub_approx(this->_upper_bounds[i],mul_approx(this->_subdivision_lengths[i],this->_upper_indices[i]));
      }
      return;
    }



    template<class R>
    dimension_type
    Grid<R>::dimension() const
    {
      return this->_subdivision_coordinates.size();
    }


    template<class R>
    Combinatoric::LatticeBlock
    Grid<R>::lattice_block() const 
    {
      return Combinatoric::LatticeBlock(this->_lower_indices,this->_upper_indices);
    }


    template<class R>
    R
    Grid<R>::subdivision_coordinate(dimension_type d, index_type n) const 
    {
      if(n > this->_upper_indices[d]) {
        return add_approx(this->_upper_origin[d],mul_approx(this->_subdivision_lengths[d],n));
      } else if(n < this->_lower_indices[d]) {
        return add_approx(this->_lower_origin[d],mul_approx(this->_subdivision_lengths[d],n));
      } else {
        return this->_centre_pointers[d][n];
      }
    }


    template<class R> 
    index_type 
    Grid<R>::subdivision_index(dimension_type d, const real_type& x) const 
    {
      index_type n=subdivision_lower_index(d,x);
      if(subdivision_coordinate(d,n) == x) { 
        return n; 
      } else {
        throw std::runtime_error("Value is not a grid coordinate");
      }
    }
    

    template<class R> 
    index_type 
    Grid<R>::subdivision_lower_index(dimension_type d, const real_type& x) const 
    {
      typename array<R>::const_iterator pos;
      if(x > this->_upper_bounds[d]) {
        return int_down<index_type>(div_down(sub_down(x,this->_upper_origin[d]),this->_subdivision_lengths[d]));
      } else if(x < this->_lower_bounds[d]) {
        return int_down<index_type>(div_down(sub_down(x,this->_lower_origin[d]),this->_subdivision_lengths[d]));
      } else {
        pos = std::upper_bound(_subdivision_coordinates[d].begin(),
                               _subdivision_coordinates[d].end(), x);
        return (pos - _centre_pointers[d])-1;
      }
    }


    template<class R> 
    index_type 
    Grid<R>::subdivision_upper_index(dimension_type d, const real_type& x) const 
    {
      typename array<R>::const_iterator pos;
      if(x > this->_upper_bounds[d]) {
        return int_up<index_type>(div_up(sub_up(x,this->_upper_origin[d]),this->_subdivision_lengths[d]));
      } else if(x < this->_lower_bounds[d]) {
        return int_up<index_type>(div_up(sub_up(x,this->_lower_origin[d]),this->_subdivision_lengths[d]));
      } else {
        pos = std::lower_bound(_subdivision_coordinates[d].begin(),
                               _subdivision_coordinates[d].end(), x);
        return (pos - _centre_pointers[d])+1;
      }
    }


    template<class R> 
    bool 
    Grid<R>::operator==(const Grid<R>& g) const
    {
      if(this==&g) { 
        return true; 
      } else {
        return this->_subdivision_coordinates==g._subdivision_coordinates
          && this->_subdivision_lengths==g._subdivision_lengths
          && this->_lower_indices==g._lower_indices;
      }
    }
  

    template<class R> 
    bool 
    Grid<R>::operator!=(const Grid<R>& g) const
    {
      return !(*this==g);
    }

                   
    template<class R>
    IndexArray 
    Grid<R>::index(const Point<R>& s) const
    {
      IndexArray res(s.dimension());
      for(size_type i=0; i!=res.size(); ++i) {
        res[i]=subdivision_index(i,s[i]);
      }
      return res;
    }


    template<class R>
    IndexArray  
    Grid<R>::lower_index(const Rectangle<R>& r) const {
      IndexArray res(r.dimension());
      for(size_type i=0; i!=res.size(); ++i) {
        res[i]=subdivision_lower_index(i,r.lower_bound(i));
      }
      return res;
    }


    template<class R>
    IndexArray  
    Grid<R>::upper_index(const Rectangle<R>& r) const {
      IndexArray res(r.dimension());
      for(size_type i=0; i!=res.size(); ++i) {
        res[i]=subdivision_upper_index(i,r.upper_bound(i));
      }
      return res;
    }


    template<class R>
    Point<R> 
    Grid<R>::point(const IndexArray& a) const
    {
      Point<R> res(a.size());
      for(size_type i=0; i!=res.dimension(); ++i) {
        res[i]=subdivision_coordinate(i,a[i]);
      }
      return res;
    }


    template<class R>
    Rectangle<R> 
    Grid<R>::rectangle(const Combinatoric::LatticeBlock& lb) const
    {
      Rectangle<R> res(lb.dimension());
      for(size_type i=0; i!=res.dimension(); ++i) {
        res.set_lower_bound(i,this->subdivision_coordinate(i,lb.lower_bound(i)));
        res.set_upper_bound(i,this->subdivision_coordinate(i,lb.upper_bound(i)));
      }
      return res;
    }


    template<class R>
    std::ostream&
    Grid<R>::write(std::ostream& os) const
    {
      os << "Grid( subdivision_coordinates=" << this->_subdivision_coordinates
         << ", subdivision_lengths=" << this->_subdivision_lengths
         << ", centre_positions=" << this->_centre_positions
         << ", origin=" << this->point(IndexArray(this->dimension(),0)) << " )";
      return os;


      
      os << "Grid( subdivisions=";
      for(dimension_type d=0; d!=this->dimension(); ++d) {
        if(d==0) { os << "[ "; } else { os << ", "; }
        os << "[ " << this->_subdivision_lengths[d] << "; " << this->_subdivision_coordinates[d][0];
        for(size_type i=1; i!=this->_subdivision_coordinates[d].size(); ++i) {
          os << ", " << this->_subdivision_coordinates[d][i];
        }
        os << "; " << this->_subdivision_lengths[d] << "]";
        os << " ]";
      }
      os << ", origin=" << this->point(IndexArray(this->dimension(),0)) << " )";
      return os;
    }


    template<class R>
    std::istream&
    Grid<R>::read(std::istream& is) 
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }






    template<class R>
    FiniteGrid<R>::~FiniteGrid() 
    {
      if(this->_own_ptr) {
        delete this->_grid_ptr;
      }
    }


    template<class R> inline
    FiniteGrid<R>::FiniteGrid(const Grid<R>& g, const Combinatoric::LatticeBlock& b) 
      : _own_ptr(0), _grid_ptr(&g), _lattice_block(b)
    { 
    }
      
    template<class R> inline
    FiniteGrid<R>::FiniteGrid(const Grid<R>& g, const Rectangle<R>& bb) 
      : _own_ptr(0), _grid_ptr(&g), _lattice_block(over_approximation(bb,g).lattice_set())
    { 
    }
      
    template<class R>
    FiniteGrid<R>::FiniteGrid(const Rectangle<R>& bb, const size_type& s)
      : _own_ptr(1), _grid_ptr(0), _lattice_block(bb.dimension())
    {
      dimension_type d=bb.dimension();
      array<R> subdivision_lengths(bb.dimension());
      for(dimension_type i=0; i!=bb.dimension(); ++i) {
        subdivision_lengths[i]=div_approx(bb[i].length(),s);
      }
      this->_grid_ptr=new Grid<R>(LinearAlgebra::Vector<R>(d,subdivision_lengths.begin()));
      GridBlock<R> bounding_box=over_approximation(bb,this->grid());
      this->_lattice_block=bounding_box.lattice_set();
    }


    template<class R>
    dimension_type 
    FiniteGrid<R>::dimension() const 
    {
      return this->_grid_ptr->dimension();
    }


    template<class R>
    Rectangle<R>
    FiniteGrid<R>::extent() const
    {
      return this->_grid_ptr->rectangle(this->_lattice_block);
    }

    
    template<class R>
    std::ostream&
    FiniteGrid<R>::write(std::ostream& os) const
    {
      return os << "FiniteGrid(grid=" << this->grid() << ", lattice_block=" << this->lattice_block() << ")";
    }


  }
}
