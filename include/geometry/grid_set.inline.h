/***************************************************************************
 *            grid_set.h
 *
 *  10 January 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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
    const Grid<R>& 
    GridCell<R>::grid() const 
    { 
      return *this->_grid_ptr; 
    }

    template<class R> inline
    dimension_type 
    GridCell<R>::dimension() const 
    {
      return this->_lattice_set.dimension(); 
    }

    template<class R> inline
    const Combinatoric::LatticeCell& 
    GridCell<R>::lattice_set() const 
    {
      return this->_lattice_set; 
    }

    template<class R> inline
    tribool 
    GridCell<R>::bounded() const { 
      return true; 
    }

    template<class R> inline
    Rectangle<R> 
    GridCell<R>::bounding_box() const 
    {
      return *this; 
    }




    template<class R> inline
    GridBlock<R>::GridBlock(const Grid<R>* gptr, const Combinatoric::LatticeBlock& lc)
      : _grid_ptr(gptr), _lattice_set(lc)
    {
    }


    template<class R> inline
    const Grid<R>& 
    GridBlock<R>::grid() const 
    {
      return *this->_grid_ptr; 
    }


    template<class R> inline
    dimension_type 
    GridBlock<R>::dimension() const 
    {
      return this->_lattice_set.dimension(); 
    }


    template<class R> inline
    const Combinatoric::LatticeBlock& GridBlock<R>::lattice_set() const 
    {
      return this->_lattice_set; 
    }


    template<class R> inline
    tribool 
    GridBlock<R>::empty() const 
    {
      return this->_lattice_set.empty(); 
    }


    template<class R> inline
    tribool 
    GridBlock<R>::bounded() const 
    {
      return true; 
    }


    template<class R> inline
    Rectangle<R> 
    GridBlock<R>::bounding_box() const 
    {
      return *this; 
    }



    template<class R> inline
    typename GridBlock<R>::const_iterator 
    GridBlock<R>::begin() const 
    {
      return const_iterator(*this->_grid_ptr,_lattice_set.begin()); 
    }


    template<class R> inline
    typename GridBlock<R>::const_iterator 
    GridBlock<R>::end() const 
    {
      return const_iterator(*this->_grid_ptr,_lattice_set.end()); 
    }



    template<class R> inline
    const Grid<R>& 
    GridCellListSet<R>::grid() const 
    {
      return *this->_grid_ptr; 
    }

    template<class R> inline
    dimension_type 
    GridCellListSet<R>::dimension() const 
    {
      return this->_lattice_set.dimension(); 
    }

    template<class R> inline
    tribool 
    GridCellListSet<R>::empty() const 
    {
      return this->_lattice_set.empty(); 
    }

    template<class R> inline
    tribool 
    GridCellListSet<R>::bounded() const 
    {
      return true; 
    }

    template<class R> inline
    size_type 
    GridCellListSet<R>::size() const 
    {
      return _lattice_set.size(); 
    }

    template<class R> inline
    const Combinatoric::LatticeCellListSet& 
    GridCellListSet<R>::lattice_set() const 
    {
      return _lattice_set; 
    }

    template<class R> inline
    GridCell<R> 
    GridCellListSet<R>::operator[] (const size_type i) const 
    {
      return GridCell<R>(grid(),_lattice_set[i]); 
    }


    template<class R> inline
    typename GridCellListSet<R>::const_iterator 
    GridCellListSet<R>::begin() const 
    {
      return const_iterator(*this->_grid_ptr,_lattice_set.begin()); 
    }


    template<class R> inline
    typename GridCellListSet<R>::const_iterator 
    GridCellListSet<R>::end() const 
    {
      return const_iterator(*this->_grid_ptr,_lattice_set.end()); 
    }


    template<class R> inline
    void
    GridCellListSet<R>::unique_sort()
    {
      this->_lattice_set.unique_sort();
    }


    template<class R> inline
    void 
    GridCellListSet<R>::adjoin(const GridCell<R>& c) 
    {
      _lattice_set.adjoin(c.lattice_set()); 
    }


    template<class R> inline
    void 
    GridCellListSet<R>::adjoin(const GridBlock<R>& bl) 
    {
      _lattice_set.adjoin(bl.lattice_set()); 
    }


    template<class R> inline
    void 
    GridCellListSet<R>::adjoin(const GridCellListSet<R>& cls) 
    {
      _lattice_set.adjoin(cls.lattice_set()); 
    }


    template<class R> template<class SetInterface> inline
    void 
    GridCellListSet<R>::adjoin_over_approximation(const SetInterface& s) 
    {
        this->adjoin(over_approximation(s,this->grid()));
    }

    
    template<class R> template<class SetInterface> inline
    void 
    GridCellListSet<R>::adjoin_under_approximation(const SetInterface& s) 
    {
        this->adjoin(under_approximation(s,this->grid()));
    }






    template<class R> inline
    bool 
    GridMaskSet<R>::operator==(const GridMaskSet<R>& gms) 
    {
      throw Deprecated(__PRETTY_FUNCTION__);
    }


    template<class R> inline
    bool 
    GridMaskSet<R>::operator!=(const GridMaskSet<R>& gms) 
    {
      return !(*this==gms);
    }

    template<class R> inline
    const Grid<R>& 
    GridMaskSet<R>::grid() const 
    {
      return *this->_grid_ptr; 
    }

    template<class R> inline
    GridBlock<R> 
    GridMaskSet<R>::bounds() const 
    {
      return GridBlock<R>(*this->_grid_ptr,_lattice_set.block()); 
    }

    template<class R> inline
    const Combinatoric::LatticeMaskSet& 
    GridMaskSet<R>::lattice_set() const 
    {
      return _lattice_set; 
    }

    template<class R> inline
    dimension_type 
    GridMaskSet<R>::dimension() const 
    {
      return _lattice_set.dimension(); 
    }

    template<class R> inline
    size_type 
    GridMaskSet<R>::capacity() const 
    {
      return _lattice_set.capacity(); 
    }
      
    template<class R> inline
    const Combinatoric::LatticeBlock& 
    GridMaskSet<R>::block() const 
    {
      return _lattice_set.block(); 
    }
  
    //template<class R> inline const SizeArray& GridMaskSet<R>::sizes() const { return _lattice_set.sizes(); } 

    template<class R> inline
    const BooleanArray& 
    GridMaskSet<R>::mask() const 
    {
      return _lattice_set.mask(); 
    }

    template<class R> inline
    tribool 
    GridMaskSet<R>::empty() const 
    {
      return _lattice_set.empty(); 
    }
      
    template<class R> inline
    tribool 
    GridMaskSet<R>::bounded() const 
    {
      return _lattice_set.bounded(); 
    }
      
    template<class R> inline
    size_type 
    GridMaskSet<R>::size() const 
    {
      return _lattice_set.size(); 
    }
      
    template<class R> inline
    GridCell<R> 
    GridMaskSet<R>::operator[](size_type i) const 
    {
      return GridCell<R>(*this->_grid_ptr,_lattice_set[i]); 
    }

    template<class R> inline
    typename GridMaskSet<R>::const_iterator 
    GridMaskSet<R>::begin() const 
    {
      return const_iterator(*this->_grid_ptr,this->_lattice_set.begin()); 
    }
    template<class R> inline
    typename GridMaskSet<R>::const_iterator 
    GridMaskSet<R>::end() const 
    {
      return const_iterator(*this->_grid_ptr,this->_lattice_set.end()); 
    }


    template<class R> inline
    void 
    GridMaskSet<R>::remove(const GridCell<R>& gc) 
    {
      ARIADNE_CHECK_SAME_GRID(*this,gc,"void GridMaskSet::remove(GridCell gc)");
      this->_lattice_set.remove(gc._lattice_set);
    }

    template<class R> inline
    void 
    GridMaskSet<R>::remove(const GridCellListSet<R>& gcls) 
    {
      ARIADNE_CHECK_SAME_GRID(*this,gcls,"void GridMaskSet::remove(GridCellListSet gcls)");
      this->_lattice_set.remove(gcls._lattice_set);
    }

    template<class R> inline
    void 
    GridMaskSet<R>::remove(const GridMaskSet<R>& gms) 
    {
      ARIADNE_CHECK_SAME_GRID(*this,gms,"void GridMaskSet::remove(GridMaskSet gms)");
      this->_lattice_set.remove(gms._lattice_set);
    }

    template<class R> inline
    void
    GridMaskSet<R>::adjoin_unbounded_cell() 
    {
      this->_lattice_set.adjoin_unbounded_cell();
    }
    
    template<class R> inline
    void 
    GridMaskSet<R>::adjoin(const GridCell<R>& gc) 
    {
      ARIADNE_CHECK_SAME_GRID(*this,gc,"void GridMaskSet::adjoin(GridCell gc)");
      this->_lattice_set.adjoin(gc._lattice_set);
    }

    template<class R> inline
    void 
    GridMaskSet<R>::adjoin(const GridBlock<R>& gb) 
    {
      ARIADNE_CHECK_SAME_GRID(*this,gb,"void GridMaskSet::adjoin(GridBlock gb)");
      this->_lattice_set.adjoin(gb._lattice_set);
    }

    template<class R> inline
    void 
    GridMaskSet<R>::adjoin(const GridCellListSet<R>& gcls) 
    {
      ARIADNE_CHECK_SAME_GRID(*this,gcls,"void GridMaskSet::adjoin(GridCellListSet gcls)");
      this->_lattice_set.adjoin(gcls._lattice_set);
    }

    template<class R> template<class S> inline
    void 
    GridMaskSet<R>::adjoin_over_approximation(const S& s) 
    {
      this->adjoin(over_approximation(s,FiniteGrid<R>(this->grid(),this->block())));
    }

    template<class R> template<class S> inline
    void 
    GridMaskSet<R>::adjoin_under_approximation(const S& s) 
    {
      this->adjoin(under_approximation(s,FiniteGrid<R>(this->grid(),this->block())));
    }

    template<class R> template<class S> inline
    void 
    GridMaskSet<R>::restrict_over_approximation(const S& s) 
    {
      this->restrict(over_approximation(s,FiniteGrid<R>(this->grid(),this->block())));
    }


    template<class R> inline
    void 
    GridMaskSet<R>::adjoin(const GridMaskSet<R>& gms) 
    {
      ARIADNE_CHECK_SAME_GRID(*this,gms,"void GridMaskSet::adjoin(GridMaskSet gms)");
      this->_lattice_set.adjoin(gms._lattice_set);
    }

    template<class R> inline
    void 
    GridMaskSet<R>::restrict(const GridCellListSet<R>& gcls) 
    {
      ARIADNE_CHECK_SAME_GRID(*this,gcls,"void GridMaskSet::restrict(GridCellListSet gcls)");
      this->_lattice_set.restrict(gcls._lattice_set);
    }

    template<class R> inline
    void 
    GridMaskSet<R>::restrict(const GridMaskSet<R>& gms) 
    {
      ARIADNE_CHECK_SAME_GRID(*this,gms,"void GridMaskSet::restrict(GridMaskSet gms)");
      this->_lattice_set.restrict(gms._lattice_set);
    }


    template<class R> inline
    GridMaskSet<R> 
    GridMaskSet<R>::neighbourhood() const 
    {
        return GridMaskSet(this->grid(),this->_lattice_set.neighbourhood());
    }


    template<class R> inline
    GridMaskSet<R> 
    GridMaskSet<R>::adjoining() const 
    {
        return GridMaskSet(this->grid(),this->_lattice_set.adjoining());
    }

    
    template<class R> template<class BS> inline
    GridMaskSet<R>::operator ListSet<BS> () const 
    {
      ListSet<BS> result(this->dimension());
      Rectangle<R> r(this->dimension());
      BS bs(this->dimension());
      for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        r=*iter;
        bs=BS(r);
        result.push_back(bs);
      }
      return result;
    }

    
    template<class R, class BS> inline
    GridMaskSet<R> 
    over_approximation(const ListSet<BS>& ls, const FiniteGrid<R>& g) 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(ls,g,"over_approximation(ListSet<R,BS>,FiniteGrid<R>)");

      GridMaskSet<R> result(g);
      for(typename ListSet<BS>::const_iterator iter=ls.begin(); iter!=ls.end(); ++iter) {
       result.adjoin(over_approximation(*iter,g.grid()));
      }
        
      return result;
    }
    
    
    
    
    
    
    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const GridCell<R>& gc) {
      return gc.write(os);
    }
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const GridBlock<R>& gb) {
      return gb.write(os);
    }
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const GridCellListSet<R>& gcls) {
      return gcls.write(os);
    }
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const GridMaskSet<R>& gms) {
      return gms.write(os);
    }
    
    template<class Base, class Value>
    class GridSetIterator 
      : public boost::iterator_adaptor<GridSetIterator<Base,Value>,Base,Value,boost::use_default,Value>
    { 
     public:
      typedef typename Value::real_type real_type;
      GridSetIterator(const Grid<real_type>& g, Base i) 
        : GridSetIterator::iterator_adaptor_(i), _grid(g) { }
     private:
      friend class boost::iterator_core_access;
      Value dereference() const { return Value(_grid,*this->base_reference()); }
      const Grid<real_type>& _grid;
    };




    template<class R>
    class GridCellListSetIterator {
      typedef GridCellListSetIterator Self;
     public:
      typedef std::forward_iterator_tag iterator_category;
      typedef GridCell<R> value_type;
      typedef GridCell<R> reference;
      typedef const GridCell<R>* pointer;
      typedef int difference_type;
     public:
      GridCellListSetIterator(const Grid<R>& g, Combinatoric::LatticeCellListSetIterator iter);
      GridCell<R> operator*() const { return this->dereference(); }
      Self& operator++() { this->increment(); return *this; }
      Self operator++(int) { Self tmp=*this; this->increment(); return tmp; }
      bool operator==(const Self& other) const { return this->equal(other); }
      bool operator!=(const Self& other) const { return !this->equal(other); }
     private:
      bool equal(const GridCellListSetIterator& other) const { return this->_iter==other._iter; }
      void increment() { ++_iter; }
      GridCell<R> dereference() const { return GridCell<R>(_grid_ref,*_iter); }
     private:
      Grid<R>& _grid_ref;
      Combinatoric::LatticeCellListSet::const_iterator _iter;
    };

    template<class R>
    class GridMaskSetIterator {
      typedef GridMaskSetIterator Self;
     public:
      typedef std::forward_iterator_tag iterator_category;
      typedef GridCell<R> value_type;
      typedef GridCell<R> reference;
      typedef const GridCell<R>* pointer;
      typedef int difference_type;
     public:
      GridMaskSetIterator(const Grid<R>& g, Combinatoric::LatticeMaskSetIterator iter);
      GridCell<R> operator*() const { return this->dereference(); }
      Self& operator++() { this->increment(); return *this; }
      Self operator++(int) { Self tmp=*this; this->increment(); return tmp; }
      bool operator==(const Self& other) const { return this->equal(other); }
      bool operator!=(const Self& other) const { return !this->equal(other); }
     private:
      bool equal(const GridMaskSetIterator& other) const { return this->_iter==other._iter; }
      void increment() { ++_iter; }
      GridCell<R> dereference() const { return GridCell<R>(_grid_ref,*_iter); }
     private:
      Grid<R>& _grid_ref;
      Combinatoric::LatticeMaskSet::const_iterator _iter;
    };

    
    
    
    
  }
}
