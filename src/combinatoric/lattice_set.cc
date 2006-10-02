/***************************************************************************
 *            lattice_set.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <algorithm>

#include "base/stlio.h"

#include "combinatoric/array_operations.h"
#include "combinatoric/lattice_set.h"

namespace Ariadne {
  namespace Combinatoric {
    
    template<typename Ptr>
    class array_less {
     public:
      array_less(size_type n) : _array_size(n) { }
      bool operator() (Ptr p1, Ptr p2) {
        Ptr p1_end=p1+_array_size;
        while(p1!=p1_end) {
          if(*p1<*p2) {
            return true;
          }
          else if(*p2<*p1) {
            return false;
          }
          ++p1;
          ++p2;
        }
        return false;
      }
     private:
      size_type _array_size;
    };
      
    template<typename Ptr>
    class array_eq {
     public:
      array_eq(size_type n) : _array_size(n) { }
      bool operator() (Ptr p1, Ptr p2) {
        Ptr p1_end=p1+_array_size;
        while(p1!=p1_end) {
          if(*p1!=*p2) {
            return false;
          }
          ++p1;
          ++p2;
        }
        return true;
      }
     private:
      size_type _array_size;
    };
      
    LatticeCell::LatticeCell(const std::string& s)
      : _lower()
    {
      std::stringstream ss(s);
      ss >> *this;
    }
     
    bool
    operator<(const LatticeCell& lc1, const LatticeCell& lc2) {
      assert(lc1.dimension()==lc2.dimension());
      for(dimension_type i=0; i!=lc1.dimension(); ++i) {
        if(lc1.lower_bound(i)<lc2.lower_bound(i)) { 
          return true;
        }
        else if(lc1.lower_bound(i)>lc2.lower_bound(i)) {
          return false;
        }
      }
      return false;
    }

    std::istream& 
    operator>>(std::istream& is, LatticeCell& lc)
    {
      char c;
      is >> c;
      if(c=='(') {
        /* Representation as a literal (l1,l2,...,ln) */
        std::vector< int > v;
        int i;
        c=',';
        while(c==',') {
          is >> i;
          v.push_back(i);
          c=' ';
          while( is && c==' ') {
            is >> c;
          }
        }
        if(is) {
          assert(c=')');
        }
        
        IndexArray l(v.size());
        for(size_type i=0; i!=v.size(); ++i) {
          l[i]=v[i];
        }
        lc=LatticeCell(l);
      }
      else {
        is.putback(c);
        /* representation as lower and upper corners */
        /* FIXME */
        // throw invalid_input("Not implemented");
      }
      return is;
    }
    
    
    
    LatticeBlock::LatticeBlock(const std::string& s)
      : _lower(), _upper()
    {
      std::stringstream ss(s);
      ss >> *this;
    }
     
    SizeArray
    LatticeBlock::sizes() const
    {
      SizeArray result(this->dimension());
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        result[i]=std::max(0,this->_upper[i]-this->_lower[i]);
      }
      return result;
    }

    SizeArray
    LatticeBlock::strides() const
    {
      SizeArray result(this->dimension()+1);
      result[0]=1;
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        result[i+1]=std::max(0,this->_upper[i]-this->_lower[i])*result[i];
      }
      return result;
    }

    size_type
    LatticeBlock::size() const
    {
      size_type result=1;
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        result*=std::max(0,this->_upper[i]-this->_lower[i]);
      }
      return result;
    }

    bool 
    LatticeBlock::empty() const
    {
      if(this->dimension()==0) {
        return true;
      }
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        if(this->lower_bound(i)>this->upper_bound(i)) {
          return true;
        }
      }
      return false;
    }
    
    bool 
    LatticeBlock::empty_interior() const
    {
      if(this->dimension()==0) {
        return true;
      }
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        if(this->lower_bound(i)>=this->upper_bound(i)) {
          return true;
        }
      }
      return false;
    }
    
    LatticeBlock
    LatticeBlock::neighbourhood() const
    {
      dimension_type n=this->dimension();
      IndexArray lower(n);
      IndexArray upper(n);
      for(dimension_type i=0; i!=n; ++i) {
        lower[i]=this->lower_bound(i)-1;
        upper[i]=this->upper_bound(i)+1;
      }
      return LatticeBlock(lower,upper);
    }

    std::istream& 
    operator>>(std::istream& is, LatticeBlock& r)
    {
      char c;
      is >> c;
      is.putback(c);
      if(c=='[') {
        /* Representation as a literal [a1,b1]x[a2,b2]x...x[an,bn] */
        std::vector< Interval<int> > v;
        Interval<int> i;
        c='x';
        while(c=='x') {
          is >> i;
          v.push_back(i);
          c=' ';
          while( is && c==' ') {
            is >> c;
          }
        }
        if(is) {
          is.putback(c);
        }
        
        IndexArray l(v.size());
        IndexArray u(v.size());
        for(size_type i=0; i!=v.size(); ++i) {
          l[i]=v[i].lower();
          u[i]=v[i].upper();
        }
        r=LatticeBlock(l,u);
      }
      else {
        /* representation as lower and upper corners */
        /* FIXME */
        // throw invalid_input("Not implemented");
      }
      return is;
    }

    
    LatticeBlock
    LatticeTransformation::operator() (const LatticeCell& c) const
    {
      dimension_type n=c.dimension();
      IndexArray lower(n);
      IndexArray upper(n);
      for(dimension_type i=0; i!=n; ++i) {
        lower[i]=_transformation[i][c.lower_bound(i)];
        upper[i]=_transformation[i][c.upper_bound(i)];
      }
      return LatticeBlock(lower,upper);
    }

    LatticeBlock
    LatticeTransformation::operator() (const LatticeBlock& r) const
    {
      dimension_type n=r.dimension();
      IndexArray lower(n);
      IndexArray upper(n);
      for(dimension_type i=0; i!=n; ++i) {
        lower[i]=_transformation[i][r.lower_bound(i)];
        upper[i]=_transformation[i][r.upper_bound(i)];
      }
      return LatticeBlock(lower,upper);
    }



    LatticeCellListSet::LatticeCellListSet(const LatticeCell& c) 
      : _list(c.dimension()) 
    {
      this->adjoin(c); 
    }
    
    LatticeCellListSet::LatticeCellListSet(const LatticeBlock& r) 
      : _list(r.dimension()) 
    {
      this->adjoin(r); 
    }
    
    LatticeCellListSet::LatticeCellListSet(const LatticeMaskSet& ms) 
      : _list(ms.dimension()) 
    {
      this->adjoin(ms); 
    }
    
    LatticeCellListSet::LatticeCellListSet(const LatticeBlockListSet& rls) 
      : _list(rls.dimension()) 
    {
      this->adjoin(rls); 
    }
    
    LatticeBlock
    LatticeCellListSet::bounding_block() const
    { 
      if(this->empty()) {
        return LatticeBlock(this->dimension());
      }

      LatticeCell c=(*this)[0];
      IndexArray lower(c.lower());
      IndexArray upper(c.upper());

      for(size_type i=1; i!=this->size(); ++i) {
        c=(*this)[i];
        assign_min(lower,c.lower());
        assign_max(upper,c.upper());
      }
      return LatticeBlock(lower,upper);
    }
    
    void
    LatticeCellListSet::unique_sort()
    {
      dimension_type n=this->dimension();
      
      std::vector<index_type*> pointers(this->size());
      index_type* ptr=&(*this->_list.begin())[0];
      for(size_type i=0; i!=this->size(); ++i) {
        pointers[i]=ptr;
        ptr+=n;
      }
      
      std::sort(pointers.begin(),pointers.end(),array_less<index_type*>(n));
      std::vector<index_type*>::iterator p=std::unique(pointers.begin(),pointers.end(), array_eq<index_type*>(n));
      pointers.erase(p,pointers.end());
      
      array_vector<index_type> sorted(n);
      sorted.resize(pointers.size());
      for(size_type i=0; i!=sorted.size(); ++i) {
        sorted[i]=range<index_type*>(pointers[i],pointers[i]+n);
      }
      
      this->_list.swap(sorted);
    }

    void 
    LatticeCellListSet::adjoin(const LatticeBlock& r) 
    { 
      assert(this->dimension() == r.dimension());
      for(LatticeBlock::const_iterator i=r.begin(); i!=r.end(); ++i) {
        this->adjoin(*i);
      }
    }

    void 
    LatticeCellListSet::adjoin(const LatticeCellListSet& cl) {
      for(LatticeCellListSet::const_iterator i=cl.begin(); i!=cl.end(); ++i) {
        this->adjoin(*i);
      }
    }

    void 
    LatticeCellListSet::adjoin(const LatticeBlockListSet& rl) {
      for(LatticeBlockListSet::const_iterator i=rl.begin(); i!=rl.end(); ++i) {
        this->adjoin(*i);
      }
    }

    void 
    LatticeCellListSet::adjoin(const LatticeMaskSet& ms) {
      assert(this->dimension() == ms.dimension());
      for(LatticeMaskSet::const_iterator i=ms.begin(); i!=ms.end(); ++i) {
        this->adjoin(*i);
      }
    }


    LatticeBlock
    LatticeBlockListSet::bounding_block() const
    { 
      if(this->empty()) {
        return LatticeBlock(this->dimension());
      }
      
      LatticeBlock r=(*this)[0];
      IndexArray lower(r.lower());
      IndexArray upper(r.upper());

      for(size_type i=1; i!=this->size(); ++i) {
        r=(*this)[i];
        assign_min(lower,r.lower());
        assign_max(upper,r.upper());
      }
      return LatticeBlock(lower,upper);
    }

    void
    LatticeBlockListSet::adjoin(const LatticeBlockListSet& rls)
    {
      for(LatticeBlockListSet::const_iterator rect_iter=rls.begin(); rect_iter!=rls.end(); ++rect_iter) {
        this->adjoin(*rect_iter);
      }
    }


    LatticeMaskSet::LatticeMaskSet(const LatticeBlock& bd, const LatticeCellListSet& cls)
      : _block(bd), _unbounded(false)
    {
      this->_compute_cached_attributes();
      this->adjoin(cls);
    }
    
    LatticeMaskSet::LatticeMaskSet(const LatticeBlock& bd, const LatticeBlockListSet& rls)
      : _block(bd), _unbounded(false)
    {
      this->_compute_cached_attributes();
      this->adjoin(rls);
    }
    
    LatticeMaskSet::LatticeMaskSet(const LatticeCellListSet& cls)
      : _block(cls.bounding_block()), _unbounded(false)
    {
      this->_compute_cached_attributes();
      this->adjoin(cls);
    }
    
    LatticeMaskSet::LatticeMaskSet(const LatticeBlockListSet& rls)
      : _block(rls.bounding_block()), _unbounded(false)
    {
      this->_compute_cached_attributes();
      this->adjoin(rls);
    }
    
    LatticeMaskSet::LatticeMaskSet(const LatticeMaskSet& ms)
      : _block(ms._block), _mask(ms._mask), _unbounded(false)
    {
      this->_compute_cached_attributes();
    }
    

    LatticeCell
    LatticeMaskSet::operator[](size_type i) const
    {
      const_iterator iter=this->begin();
      while(i!=0) {
        --i;
        ++iter;
      }
      return *iter;
    }

    void
    LatticeMaskSet::clear()
    {
      _mask=BooleanArray(_mask.size(),false);
    }



    /* Compute the index of a lattice cell in a grid. */
    size_type
    LatticeMaskSet::index(const LatticeCell& c) const
    {
      return this->_index(c.position());
    }
    
     /* Compute the index of a lattice cell in a grid. */
    size_type
    LatticeMaskSet::_index(const IndexArray& pos) const
    {
      size_type result=0;
      for(dimension_type i=0; i!=pos.size(); ++i) {
        result += size_type(pos[i]-this->_lower[i])*this->_strides[i];
      }
      return result;
    }

 
/*
    size_type
    LatticeMaskSet::index(const IndexArray& pos) const
    {
      index_type result=this->_origin_index;
      dimension_type n=pos.size();
      for(dimension_type i=0; i!=n; ++i) {
        result += pos[i]*this->_strides[i];
      }
      return size_type(result);
    }
*/

    LatticeCell
    LatticeMaskSet::cell(size_type index) const
    {
      LatticeCell result(this->dimension());
      IndexArray& array=result._lower;
      dimension_type n=this->dimension();
      for(dimension_type i=n-1; i!=0; --i) {
        array[i] = index/this->_strides[i]+this->_lower[i];
        index = index%this->_strides[i];
      }
      array[0]=index;
      return result;
    }


    void 
    LatticeMaskSet::adjoin(const LatticeBlock& b) 
    {
      LatticeBlock r=regular_intersection(this->block(),b);
      if(!subset(b,this->block())) {
        this->_unbounded=true;
      }
      
      dimension_type n=this->dimension();
      const IndexArray& rlower(r.lower());
      const IndexArray& rupper(r.upper());
      const IndexArray& glower(this->block().lower());
      const SizeArray gstrides(this->block().strides());

      if(n==1) {
        for(size_type i=rlower[0]-glower[0]; 
          i!=size_type(rupper[0]-glower[0]); ++i) {
          _mask[i]=true;
        }
        return;
      }

      if(n==2) {
        SizeArray rsizes=r.sizes();
        size_type index=this->_index(rlower);
        for(size_type loop_end=index+rsizes[1]*gstrides[1]; index!=loop_end; index+=gstrides[1]-rsizes[0]) {
          for(size_type inner_loop_end=index+rsizes[0]; index!=inner_loop_end; index+=1) {
            _mask[index]=true;
          }
        }

        return;
      }


      if(n==0) {
        return;
      }


      /* dim>2 */
      SizeArray rsizes=r.sizes();
      size_type index=this->_index(rlower);
      IndexArray rposition=rlower;

      while(rposition[n-1]!=rupper[n-1]) {
        _mask[index]=true;
        
        dimension_type d=0;
        rposition[d]+=1;
        index+=gstrides[d];
        while(rposition[d]==rupper[d] && (d+1u)!=n ) {
          rposition[d]=rlower[d];
          index-=gstrides[d]*rsizes[d];
          d+=1;
          rposition[d]+=1;
          index+=gstrides[d];
        }
      }
    }

    void 
    LatticeMaskSet::adjoin(const LatticeCellListSet& cl) {
      for(LatticeCellListSet::const_iterator i=cl.begin(); i!=cl.end(); ++i) {
        if(subset(*i,this->block())) {
          this->adjoin(*i);
        } else {
          this->_unbounded=true;
        }
      }
    }

    void 
    LatticeMaskSet::adjoin(const LatticeBlockListSet& rl) {
      for(LatticeBlockListSet::const_iterator i=rl.begin(); i!=rl.end(); ++i) {
        this->adjoin(*i);
      }
    }

    void 
    LatticeMaskSet::adjoin(const LatticeMaskSet& lm) 
    {
      if(this->_block==lm._block) {
        this->_mask |= lm._mask;
      }
      else {
        for(LatticeMaskSet::const_iterator iter=lm.begin(); iter!=lm.end(); ++iter) {
          this->adjoin(*iter);
        }
      }
    }

    void
    LatticeMaskSet::_compute_cached_attributes() 
    {
      if(_block.empty()) {
        IndexArray origin(_block.dimension(),0);
        _block=LatticeBlock(origin,origin);
      }
      _lower=_block.lower();
      _upper=_block.upper();
      _sizes=_block.sizes();
      _strides=_block.strides();
      _mask.resize(_strides[_block.dimension()]);
    }
    


    LatticeBlock 
    regular_intersection(const LatticeBlock& r1, const LatticeBlock& r2) 
    {
      LatticeBlock result(r1.dimension());
      for(dimension_type i=0; i!=r1.dimension(); ++i) {
        result.set_lower_bound(i,std::max(r1.lower_bound(i),r2.lower_bound(i)));
        result.set_upper_bound(i,std::min(r1.upper_bound(i),r2.upper_bound(i)));
      }
      return result;
    }

    LatticeMaskSet
    regular_intersection(const LatticeMaskSet& A, const LatticeMaskSet& B) 
    {
      assert(A.block()==B.block());
      return LatticeMaskSet(A.block(),A.mask() & B.mask(), !A.bounded() & !B.bounded());
    }

    LatticeMaskSet
    join(const LatticeMaskSet& A, const LatticeMaskSet& B) 
    {
      assert(A.block()==B.block());
      return LatticeMaskSet(A.block(),A.mask() | B.mask(), !A.bounded() | !B.bounded());
    }

    LatticeMaskSet
    difference(const LatticeMaskSet& A, const LatticeMaskSet& B) 
    {
      assert(A.block()==B.block());
      return LatticeMaskSet(A.block(),A.mask() - B.mask(), !A.bounded() - !B.bounded());
    }

    bool 
    disjoint(const LatticeBlock& r1, const LatticeBlock& r2) 
    {
      for(dimension_type i=0; i!=r1.dimension(); ++i) {
        if(r1.upper_bound(i)<r2.lower_bound(i)
            || r1.upper_bound(i)<r2.lower_bound(i))
        {
          return true;
        }
      }
      return false;
    }
    
    bool 
    disjoint(const LatticeBlock& r, const LatticeMaskSet& ms) 
    {
      return !interiors_intersect(r.neighbourhood(),ms);
    }
    
    bool 
    disjoint(const LatticeMaskSet& ms1, const LatticeMaskSet& ms2) 
    {
      return !interiors_intersect(ms1.neighbourhood(),ms2);
    }
    
    bool 
    interiors_intersect(const LatticeBlock& r1, const LatticeBlock& r2) 
    {
      for(dimension_type i=0; i!=r1.dimension(); ++i) {
        if(r1.upper_bound(i)<=r2.lower_bound(i)
            || r1.upper_bound(i)<=r2.lower_bound(i))
        {
          return false;
        }
      }
      return true;
    }
    
    bool 
    interiors_intersect(const LatticeBlock& r, const LatticeMaskSet& ms) 
    {
      LatticeBlock rstr=regular_intersection(r,ms.block());
      if(rstr.empty()) {
        return false;
      }
      for(LatticeBlock::const_iterator i=rstr.begin(); i!=rstr.end(); ++i) {
        if(subset(*i,ms)) {
          return true;
        }
      }
      return false;
    }
    
    bool 
    interiors_intersect(const LatticeMaskSet& a, const LatticeMaskSet& b) 
    {
      if(a.block()==b.block()) {
        BooleanArray::const_iterator aiter=a.mask().begin();
        BooleanArray::const_iterator biter=a.mask().begin();
        BooleanArray::const_iterator aend=a.mask().end();
        while(aiter!=aend) {
          if(*aiter & *biter) {
            return true;
          }
          ++aiter;
          ++biter;
        }
        return false;
      } else {
        for(LatticeMaskSet::const_iterator aiter=a.begin(); aiter!=a.end(); ++aiter) {
          if(interiors_intersect(*aiter,b)) {
            return true;
          }
        }
        return false;
      }
    }
    

    bool
    subset(const LatticeCell& lc, const LatticeBlock& lr)
    {
      return subset(LatticeBlock(lc),lr);
    }
     

    bool 
    subset(const LatticeCell& c, const LatticeMaskSet& ms) 
    {
      return ms.mask()[ms.index(c)];
    }
    
    bool 
    subset(const LatticeBlock& r1, const LatticeBlock& r2) 
    {
      if(r1.empty()) { 
        return true; 
      }
      for(dimension_type i=0; i!=r1.dimension(); ++i) {
        if(r1.lower_bound(i)<r2.lower_bound(i)
            || r1.upper_bound(i)>r2.upper_bound(i))
        {
          return false;
        }
      }
      return true;
    }

    bool 
    subset(const LatticeBlock& r, const LatticeMaskSet& ms) 
    {
      if(r.empty()) {
        return true;
      }
      if(!subset(r,ms.block())) {
        return false;
      }
      for(LatticeBlock::const_iterator i=r.begin(); i!=r.end(); ++i) {
        if(!subset(*i,ms)) {
          return false;
        }
      }
      return true;
    }
    
    bool 
    subset(const LatticeMaskSet& a, const LatticeMaskSet& b) {
      if(a.block() == b.block()) {
        return a.mask() <= b.mask();
      } else {
        for(LatticeMaskSet::const_iterator aiter=a.begin(); aiter!=a.end(); ++aiter) {
          if(!subset(*aiter,b)) {
            return false;
          }
        }
        return true;
      }
    }
      

    LatticeMaskSet LatticeMaskSet::neighbourhood() const {
      LatticeMaskSet result=*this;
      const IndexArray& lower_bound=result.block().lower();
      const IndexArray& upper_bound=result.block().upper();
      for(LatticeMaskSet::const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        LatticeCell cell=*iter;
        IndexArray lower(this->dimension());
        IndexArray upper(this->dimension());
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          lower[i]=std::max(cell.lower()[i]-1,lower_bound[i]);
          upper[i]=std::min(cell.upper()[i]+1,upper_bound[i]);
        }
        LatticeBlock block(lower,upper);
        result.adjoin(block);
      }
      return result;
    }

    LatticeMaskSet LatticeMaskSet::adjoining() const {
      LatticeMaskSet result=*this;
      const IndexArray& lower_bound=result.block().lower();
      const IndexArray& upper_bound=result.block().upper();
      for(LatticeMaskSet::const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        LatticeCell cell=*iter;
        IndexArray lower=(cell.lower());
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          if(lower[i]>lower_bound[i]) {
            lower[i]-=1;
            result.adjoin(LatticeCell(lower));
            lower[i]+=1;
          }
          if(lower[i]<upper_bound[i]-1) {
            lower[i]+=1;
            result.adjoin(LatticeCell(lower));
            lower[i]-=1;
          }
        }
      }
      return result;
    }

    std::ostream& 
    operator<<(std::ostream& os, const LatticeCell& c) 
    {
      if(c.dimension()==0) {
        return os<<"Empty";
      }
      os << c[0];
      for(dimension_type i=1; i!=c.dimension(); ++i) {
        os << "x" << c[i];
      }
      return os;
    }
        
    std::ostream& 
    operator<<(std::ostream& os, const LatticeBlock& r) 
    {
      if(r.empty() || r.dimension()==0) {
        return os<<"Empty";
      }
      os << r[0];
      for(dimension_type i=1; i!=r.dimension(); ++i) {
        os << "x" << r[i];
      }
      return os;
    }
    
    std::ostream& 
    operator<<(std::ostream& os, const LatticeMaskSet& ms) 
    {
      return os << "LatticeMaskSet(\n  block=" << ms.block() << "\n  mask=" << ms.mask() << "\n)\n";
    }
        
    std::ostream& 
    operator<<(std::ostream& os, const LatticeCellListSet& cls) 
    {
      return Utility::write_sequence(os,cls.begin(),cls.end(),'[',']',',');
    }
    
    std::ostream& 
    operator<<(std::ostream& os, const LatticeBlockListSet& rls) 
    {
      return Utility::write_sequence(os,rls.begin(),rls.end(),'[',']',',');
    }
    
    
    
    class chompfstream : public std::ofstream { };
    
    chompfstream& 
    operator<<(chompfstream& cfs, const LatticeMaskSet& ms) 
    {
      std::ostream& os=cfs;
      dimension_type d=ms.dimension();
      for(LatticeMaskSet::const_iterator iter=ms.begin(); iter!=ms.end(); ++iter) {
        LatticeCell c=*iter;
        for(size_type i=0; i!=d; ++i) {
          uint k=c.lower_bound(i);
          os << k << " ";
        }
        os << "\n";
      }
      return cfs;
    }
    
    
  } 
}
