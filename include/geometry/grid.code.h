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

#include "../base/stlio.h"

#include "../numeric/arithmetic.h"

#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"

#include "../geometry/grid_set.h" // for GridBlock<R>

namespace Ariadne {
  namespace Geometry {

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
    IrregularGrid<R>*
    IrregularGrid<R>::clone() const
    {
      std::cerr << "WARNING: Cloning FiniteGrid<R> causes memory leak" << std::endl;
      return new IrregularGrid<R>(*this);
    }

    template<class R>
    IrregularGrid<R>::~IrregularGrid()
    {
    }

    template<class R>
    IrregularGrid<R>::IrregularGrid(const Rectangle<R>& r, size_type n)
      : _subdivision_coordinates(r.dimension())
    {
      for(dimension_type i=0; i!=r.dimension(); ++i) {
        R lower(r.lower_bound(i));
        R upper(r.upper_bound(i));
        R step=div_approx(sub_approx(upper,lower),(long uint)n);
        _subdivision_coordinates[i].push_back(lower);
        for(size_type j=1; j!=n; ++j) {
          _subdivision_coordinates[i].push_back(add_approx(lower,mul_approx((long uint)j,step)));
        }
        _subdivision_coordinates[i].push_back(upper);
      }
      create();
    }

    template<class R>
    IrregularGrid<R>::IrregularGrid(const Rectangle<R>& r, SizeArray sz)
      : _subdivision_coordinates(r.dimension())
    {
      for(dimension_type i=0; i!=r.dimension(); ++i) {
        R lower(r.lower_bound(i));
        R upper(r.upper_bound(i));
        R step=div_approx(sub_approx(upper,lower),(long uint)sz[i]);
        _subdivision_coordinates[i].push_back(lower);
        for(size_type j=1; j!=sz[i]; ++j) {
          _subdivision_coordinates[i].push_back(add_approx(lower,mul_approx((long uint)j,step)));
        }
        _subdivision_coordinates[i].push_back(upper);
      }
      create();
    }

    template<class R>
    IrregularGrid<R>::IrregularGrid(const array< std::vector<R> >& sp)
      : _subdivision_coordinates(sp)
    {
      create();
    }

    template<class R>
    IrregularGrid<R>::IrregularGrid(const ListSet<R,Rectangle>& ls)
      : _subdivision_coordinates(ls.dimension())
    {
      for(typename ListSet<R,Rectangle>::const_iterator riter=ls.begin(); riter!=ls.end(); ++riter) {
        for(dimension_type n=0; n!=ls.dimension(); ++n) {
          _subdivision_coordinates[n].push_back(riter->lower_bound(n));
          _subdivision_coordinates[n].push_back(riter->upper_bound(n));
        }
      }
      create();
    }

    template<class R>
    IrregularGrid<R>::IrregularGrid(const IrregularGrid<R>& g1, IrregularGrid<R>& g2)
      : _subdivision_coordinates(g1.dimension())
    {
      for(dimension_type d=0; d!=dimension(); ++d) {
        std::vector<R>& sc(_subdivision_coordinates[d]);
        const std::vector<R>& sc1(g1._subdivision_coordinates[d]);
        const std::vector<R>& sc2(g2._subdivision_coordinates[d]);
        sc.resize(sc1.size()+sc2.size());
        std::merge(sc1.begin(),sc1.end(),sc2.begin(),sc2.end(),sc.begin());
      }
      create();
    }

    template<class R>
    void
    IrregularGrid<R>::create()
    {
       for(dimension_type i=0; i!=dimension(); ++i) {
        std::vector<R>& pos=_subdivision_coordinates[i];
        std::sort(pos.begin(),pos.end());
        typename std::vector<R>::iterator newend=std::unique(pos.begin(),pos.end());
        pos.resize(std::distance(pos.begin(),newend));
      }
    }

    template<class R>
    GridBlock<R>
    IrregularGrid<R>::extent() const
    { 
      return GridBlock<R>(*this,this->lattice_block());
    }
    
    template<class R>
    array< std::vector<index_type> >
    IrregularGrid<R>::index_translation(const IrregularGrid<R>& from, const IrregularGrid<R>& to)
    {
      check_equal_dimensions(from,to,"IrregularGrid<R>::index_translation(IrregularGrid<R>,IrregularGrid<R>)");
      array< std::vector<index_type> > result(from.dimension());
      for(dimension_type d=0; d!=from.dimension(); ++d) {
        for(size_type n=0; n!=from.size(d); ++n) {
          index_type i=to.subdivision_index(d,from.subdivision_coordinate(d,n));
          result[d].push_back(i);
        }
      }
      return result;
    }


    template<class R>
    RegularGrid<R>*
    RegularGrid<R>::clone() const
    {
      std::cerr << "WARNING: Cloning InfiniteGrid<R> causes memory leak" << std::endl;
      return new RegularGrid<R>(*this);
    }

    template<class R>
    inline bool 
    Grid<R>::operator==(const Grid<R>& g) const 
    { 
       return this==&g;
    }
            
    template<class R>
    inline bool 
    Grid<R>::operator!=(const Grid<R>& g) const 
    { 
       return !(*this==g);
    }

    // FIXME: Use dynamic_cast<> 
    template<class R>
    bool 
    IrregularGrid<R>::operator==(const Grid<R>& g) const 
    {
      if (g.type()!=IRREGULAR) {
         return false;
      }
      return this->_subdivision_coordinates==
            ((const IrregularGrid<R>&)(g))._subdivision_coordinates; 
    }
      
    template<class R>
    inline bool 
    IrregularGrid<R>::operator!=(const Grid<R>& g) const 
    { 
      return !(*this==g);
    }

    // FIXME: Use dynamic_cast<> 
    template<class R>
    bool 
    RegularGrid<R>::operator==(const Grid<R>& g) const 
    { 
      if (g.type()!=REGULAR) {
        return false;
      }
      return this->_subdivision_lengths==
          ((const RegularGrid<R>&)g)._subdivision_lengths; 
    }
            
    template<class R>
    bool 
    RegularGrid<R>::operator!=(const Grid<R>& g) const 
    { 
      return !(*this==g);
    }

    template<class R>
    FiniteGrid<R>::FiniteGrid(const Rectangle<R>& bb, const size_type& s)
      : _grid_ptr(0), _grid_type(REGULAR), _lattice_block(bb.dimension())
    {
      array<R> subdivision_lengths(bb.dimension());
      for(dimension_type i=0; i!=bb.dimension(); ++i) {
        subdivision_lengths[i]=div_approx(bb[i].length(),s);
      }
      this->_grid_ptr=new RegularGrid<R>(subdivision_lengths);
      GridBlock<R> bounding_box=over_approximation(bb,*_grid_ptr);
      this->_lattice_block=bounding_box.lattice_set();
    }

    template<class R>
    Rectangle<R>
    FiniteGrid<R>::extent() const
    {
      dimension_type dim=this->dimension();
      const Combinatoric::LatticeBlock &block=this->_lattice_block;
       
      switch(this->_grid_type) {
        case REGULAR:
        {
          const RegularGrid<R> *this_grid=((RegularGrid<R> *)(this->_grid_ptr));
          
          Point<R> l(dim),u(dim);
          
          for (dimension_type i=0; i< dim; i++) {
            l[i]=mul_approx(this_grid->subdivision_length(i),block.lower_bound(i));
            u[i]=mul_approx(this_grid->subdivision_length(i),block.upper_bound(i));
          }
          
          return Rectangle<R>(l,u);
        }
        case IRREGULAR:
        {
          const IrregularGrid<R> *this_grid=((IrregularGrid<R> *)(this->_grid_ptr));
          return (Rectangle<R>)(this_grid->extent());
        }
        default:
          throw std::runtime_error("FiniteGrid<R>::bounding_box(): not implemented for this grid type.");
      }
           
      return Rectangle<R>(0);
    }

    template<class R>
    dimension_type 
    FiniteGrid<R>::dimension() const 
    {
      return this->_grid_ptr->dimension();
    }
    

    
    template<class R>
    std::ostream&
    IrregularGrid<R>::write(std::ostream& os) const
    {
      return os << "IrregularGrid(" << this->_subdivision_coordinates << ")";
    }

    
    template<class R>
    std::ostream&
    RegularGrid<R>::write(std::ostream& os) const
    {
        return os << "RegularGrid( subdivision_lengths=" 
                  << this->_subdivision_lengths << " )\n";
    }

    
    template<class R>
    std::istream&
    IrregularGrid<R>::read(std::istream& is)
    {
      char c;
      char grid_string[]="IrregularGrid(";

      try {
        int i=0; 
        do {
          is >> c;
          if(c != grid_string[i]) {
             throw std::ios_base::failure("Ariadne::Geometry::IrregularGrid<R>::read: The input object is not an IrregularGrid" );
          }
          i++;
        } while (grid_string[i]!='\0');
        
        is >> this->_subdivision_coordinates >> c;

        if (c != ')') {
             throw std::ios_base::failure("Ariadne::Geometry::IrregularGrid<R>::read: The input IrregularGrid has not a proper format" );
        }

      } catch(...) {
        throw; 
      }

      return is; 
    }

    
    template<class R>
    std::istream&
    RegularGrid<R>::read(std::istream& is)
    {
      char c;
      char grid_string[]="RegularGrid( subdivision_lengths=";

      try {
        int i=0; 
        do {
          is >> c;
          if(c != grid_string[i]) {
             throw std::ios_base::failure("Ariadne::Geometry::RegularGrid<R>::read: The input object is not an RegularGrid" );
          }
          i++;
        } while (grid_string[i]!='\0');
        
        is >> this->_subdivision_lengths >> c;

        if (c != ')') {
             throw std::ios_base::failure("Ariadne::Geometry::RegularGrid<R>::read: The input RegularGrid has not a proper format" );
        }

      }
      catch(...) {
        throw; 
      }

      return is; 
    }

    
    
    template<class R>
    std::ostream&
    FiniteGrid<R>::write(std::ostream& os) const
    {
      return os << "FiniteGrid(grid=" << this->grid() << ", lattice_block=" << this->lattice_block() << ")";
    }


  }
}
