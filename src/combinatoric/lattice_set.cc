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

#include "base/sequence_io.h"

#include "combinatoric/array_operations.h"
#include "combinatoric/lattice_set.h"

namespace Ariadne {

    
    template<class Ptr>
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
      
    template<class Ptr>
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
        result[i]=std::max((index_type)0,this->_upper[i]-this->_lower[i]);
      }
      return result;
    }

    SizeArray
    LatticeBlock::strides() const
    {
      SizeArray result(this->dimension()+1);
      result[0]=1;
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        result[i+1]=std::max((index_type)0,this->_upper[i]-this->_lower[i])*result[i];
      }
      return result;
    }

    size_type
    LatticeBlock::size() const
    {
      size_type result=1;
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        result*=std::max((index_type)0,this->_upper[i]-this->_lower[i]);
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
    LatticeTransformation::operator() (const LatticeCell& lc) const
    {
      dimension_type n=lc.dimension();
      IndexArray lower(n);
      IndexArray upper(n);
      for(dimension_type i=0; i!=n; ++i) {
        lower[i]=_transformation[i][lc.lower_bound(i)];
        upper[i]=_transformation[i][lc.upper_bound(i)];
      }
      return LatticeBlock(lower,upper);
    }

    LatticeBlock
    LatticeTransformation::operator() (const LatticeBlock& lb) const
    {
      dimension_type n=lb.dimension();
      IndexArray lower(n);
      IndexArray upper(n);
      for(dimension_type i=0; i!=n; ++i) {
        lower[i]=_transformation[i][lb.lower_bound(i)];
        upper[i]=_transformation[i][lb.upper_bound(i)];
      }
      return LatticeBlock(lower,upper);
    }



    LatticeCellListSet::LatticeCellListSet(const LatticeCell& lc) 
      : _list(lc.dimension()) 
    {
      this->adjoin(lc); 
    }
    
    LatticeCellListSet::LatticeCellListSet(const LatticeBlock& lb) 
      : _list(lb.dimension()) 
    {
      this->adjoin(lb); 
    }
    
    LatticeCellListSet::LatticeCellListSet(const LatticeMaskSet& ms) 
      : _list(ms.dimension()) 
    {
      this->adjoin(ms); 
    }
    
    LatticeBlock
    LatticeCellListSet::bounding_block() const
    { 
      if(this->empty()) {
        return LatticeBlock(this->dimension());
      }

      LatticeCell lc=(*this)[0];
      IndexArray lower(lc.lower_corner());
      IndexArray upper(lc.upper_corner());

      for(size_type i=1; i!=this->size(); ++i) {
        lc=(*this)[i];
        assign_min(lower,lc.lower_corner());
        assign_max(upper,lc.upper_corner());
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
        sorted[i]=range_array<index_type*>(pointers[i],pointers[i]+n);
      }
      
      this->_list.swap(sorted);
    }

    void 
    LatticeCellListSet::adjoin(const LatticeBlock& lb) 
    { 
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,lb,"void LatticeCellListSet::adjoin(LatticeBlock)");
      for(LatticeBlock::const_iterator i=lb.begin(); i!=lb.end(); ++i) {
        this->adjoin(*i);
      }
    }

    void 
    LatticeCellListSet::adjoin(const LatticeCellListSet& lcls) {
      for(LatticeCellListSet::const_iterator citer=lcls.begin(); citer!=lcls.end(); ++citer) {
        this->adjoin(*citer);
      }
    }

    void 
    LatticeCellListSet::adjoin(const LatticeMaskSet& lms) {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,lms,"void LatticeCellListSet::adjoin(LatticeMaskSet)");
      for(LatticeMaskSet::const_iterator citer=lms.begin(); citer!=lms.end(); ++citer) {
        this->adjoin(*citer);
      }
    }

    void 
    LatticeCellListSet::restrict(const LatticeBlock& lb) 
    { 
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,lb,"void LatticeCellListSet::restrict(LatticeBlock)");
      LatticeCellListSet result(this->dimension());
      for(LatticeCellListSet::const_iterator i=this->begin(); i!=this->end(); ++i) {
        if(subset(*i,lb)) {
          result.adjoin(*i);
        }
      }
      this->swap(result);
    }



    void 
    LatticeCellListSet::restrict(LatticeCellListSet& lcl) 
    { 
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,lcl,"void LatticeCellListSet::restrict(LatticeBlock)");
      lcl.unique_sort();
      LatticeCellListSet result(this->dimension());
      for(LatticeCellListSet::const_iterator i=this->begin(); i!=this->end(); ++i) {
        if(std::binary_search(lcl.begin(),lcl.end(),*i)) {
          result.adjoin(*i);
        }
      }
      this->swap(result);
    }

    void 
    LatticeCellListSet::remove(LatticeCellListSet& lcl) 
    { 
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,lcl,"void LatticeCellListSet::restrict(LatticeBlock)");
      lcl.unique_sort();
      LatticeCellListSet result(this->dimension());
      for(LatticeCellListSet::const_iterator i=this->begin(); i!=this->end(); ++i) {
        if(!std::binary_search(lcl.begin(),lcl.end(),*i)) {
          result.adjoin(*i);
        }
      }
      this->swap(result);
    }

    void 
    LatticeCellListSet::restrict(const LatticeMaskSet& lms) 
    { 
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,lms,"void LatticeCellListSet::restrict(LatticeBlock)");
      this->unique_sort();
      LatticeCellListSet result(this->dimension());
      for(LatticeCellListSet::const_iterator i=this->begin(); i!=this->end(); ++i) {
        if(subset(*i,lms)) {
          result.adjoin(*i);
        }
      }
      this->swap(result);
    }

    void 
    LatticeCellListSet::remove(const LatticeMaskSet& lms) 
    { 
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,lms,"void LatticeCellListSet::restrict(LatticeBlock)");
      this->unique_sort();
      LatticeCellListSet result(this->dimension());
      for(LatticeCellListSet::const_iterator i=this->begin(); i!=this->end(); ++i) {
        if(!subset(*i,lms)) {
          result.adjoin(*i);
        }
      }
      this->swap(result);
    }

    void
    LatticeCellListSet::clear()
    {
      _list.clear();
    }




    LatticeMaskSet::LatticeMaskSet(const LatticeBlock& lb, const LatticeCellListSet& lcls)
      : _block(lb), _unbounded(false)
    {
      this->_compute_cached_attributes();
      this->adjoin(lcls);
    }
    
    LatticeMaskSet::LatticeMaskSet(const LatticeCellListSet& lcls)
      : _block(lcls.bounding_block()), _unbounded(false)
    {
      this->_compute_cached_attributes();
      this->adjoin(lcls);
    }
    
    LatticeMaskSet::LatticeMaskSet(const LatticeMaskSet& lms)
      : _block(lms._block), _mask(lms._mask), _unbounded(false)
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
    LatticeMaskSet::index(const LatticeCell& lc) const
    {
      return this->_index(lc.position());
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
    LatticeMaskSet::remove(const LatticeCellListSet& lcls) 
    {
      for(LatticeCellListSet::const_iterator i=lcls.begin(); i!=lcls.end(); ++i) {
        if(subset(*i,this->block())) {
          this->remove(*i);
        }
      }
    }


    void 
    LatticeMaskSet::remove(const LatticeMaskSet& lms) 
    {
      for(LatticeMaskSet::const_iterator i=lms.begin(); i!=lms.end(); ++i) {
        if(subset(*i,this->block())) {
          this->remove(*i);
        }
      }
    }


    void 
    LatticeMaskSet::adjoin(const LatticeBlock& lb) 
    {
      if(lb.empty()) {
        return;
      }

      LatticeBlock r=regular_intersection(this->block(),lb);
      if(!subset(lb,this->block())) {
        this->_unbounded=true;
      }
      
      dimension_type n=this->dimension();
      const IndexArray& rlower(r.lower_corner());
      const IndexArray& rupper(r.upper_corner());
      const IndexArray& glower(this->block().lower_corner());
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
    LatticeMaskSet::adjoin(const LatticeCellListSet& lcls) {
      for(LatticeCellListSet::const_iterator i=lcls.begin(); i!=lcls.end(); ++i) {
        if(subset(*i,this->block())) {
          this->adjoin(*i);
        } else {
          this->_unbounded=true;
        }
      }
    }

    void 
    LatticeMaskSet::adjoin(const LatticeMaskSet& lms) 
    {
      if(this->_block==lms._block) {
        this->_mask |= lms._mask;
      }
      else {
        for(LatticeMaskSet::const_iterator iter=lms.begin(); iter!=lms.end(); ++iter) {
          this->adjoin(*iter);
        }
      }
    }

    
    void 
    LatticeMaskSet::restrict(const LatticeCellListSet& lcls) {
      LatticeMaskSet copy=*this;
      this->clear();
      for(LatticeCellListSet::const_iterator iter=lcls.begin(); iter!=lcls.end(); ++iter) {
        if(subset(*iter,copy)) {
          this->adjoin(*iter);
        }
      }
    }

    void 
    LatticeMaskSet::restrict(const LatticeMaskSet& lms) {
      if(this->_block==lms._block) {
        this->_unbounded &= lms._unbounded;
        this->_mask &= lms._mask;
      }
      else {
        LatticeMaskSet copy=*this;
        this->clear();
        for(LatticeMaskSet::const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
          if(subset(*iter,lms)) {
            this->adjoin(*iter);
          }
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
      _lower=_block.lower_corner();
      _upper=_block.upper_corner();
      _sizes=_block.sizes();
      _strides=_block.strides();
      _mask.resize(_strides[_block.dimension()]);
    }
    


    LatticeBlock 
    rectangular_hull(const LatticeBlock& lb1, const LatticeBlock& lb2) 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(lb1,lb2,"LatticeBlock::rectangular_hull(LatticeBlock lb1, LatticeBlock lb2)");
      LatticeBlock result(lb1.dimension());
      for(dimension_type i=0; i!=lb1.dimension(); ++i) {
        result.set_lower_bound(i,std::min(lb1.lower_bound(i),lb2.lower_bound(i)));
        result.set_upper_bound(i,std::max(lb1.upper_bound(i),lb2.upper_bound(i)));
      }
      return result;
    }

    LatticeBlock 
    regular_intersection(const LatticeBlock& lb1, const LatticeBlock& lb2) 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(lb1,lb2,"LatticeBlock::rectangular_hull(LatticeBlock lb1, LatticeBlock lb2)");
      LatticeBlock result(lb1.dimension());
      for(dimension_type i=0; i!=lb1.dimension(); ++i) {
        result.set_lower_bound(i,std::max(lb1.lower_bound(i),lb2.lower_bound(i)));
        result.set_upper_bound(i,std::min(lb1.upper_bound(i),lb2.upper_bound(i)));
      }
      return result;
    }

    LatticeCellListSet
    regular_intersection(const LatticeCellListSet& lcls, const LatticeMaskSet& lms) 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(lcls,lms,"LatticeCellListSet::rectangular_hull(LatticeCellListSet lcls, LatticeMaskSet lb2)");
      LatticeCellListSet result(lcls.dimension());
      for(LatticeCellListSet::const_iterator citer=lcls.begin(); citer!=lcls.end(); ++citer) {
        if(subset(*citer,lms)) {
          result.adjoin(*citer);
        }
      }
      return result;
    }

    LatticeCellListSet
    regular_intersection(const LatticeMaskSet& lms, const LatticeCellListSet& lcls) 
    {
      return regular_intersection(lcls,lms);
    }

    LatticeMaskSet
    regular_intersection(const LatticeBlock& lb, const LatticeMaskSet& lms) 
    {
      LatticeMaskSet res(lms.block()); 
      res.adjoin(lb);
      return regular_intersection(res,lms);
    }

    LatticeMaskSet
    regular_intersection(const LatticeMaskSet& lms, const LatticeBlock& lb) 
    {
      return regular_intersection(lb,lms);
    }

    LatticeMaskSet
    regular_intersection(const LatticeMaskSet& lms1, const LatticeMaskSet& lms2) 
    {
      if(lms1.block()==lms2.block()) { 
        return LatticeMaskSet(lms1.block(),lms1.mask() & lms2.mask(), !lms1.bounded() && !lms2.bounded());
      } else {
        ARIADNE_CHECK_BOUNDED(lms1,"LatticeMaskSet regular_intersection(LatticeMaskSet lms1, LatticeMaskSet lms2)");
        ARIADNE_CHECK_BOUNDED(lms2,"LatticeMaskSet regular_intersection(LatticeMaskSet lms1, LatticeMaskSet lms2)");
        LatticeMaskSet res(regular_intersection(lms1.block(),lms2.block()));
        res.adjoin(regular_intersection(LatticeCellListSet(lms1),lms2));
        return res;
      }
    }

    LatticeMaskSet
    join(const LatticeMaskSet& lms1, const LatticeMaskSet& lms2) 
    {
      if(lms1.block()==lms2.block()) { 
        return LatticeMaskSet(lms1.block(),lms1.mask() | lms2.mask(), !lms1.bounded() | !lms2.bounded());
      } else {
        ARIADNE_CHECK_BOUNDED(lms1,"LatticeMaskSet join(LatticeMaskSet lms1, LatticeMaskSet lms2)");
        ARIADNE_CHECK_BOUNDED(lms2,"LatticeMaskSet join(LatticeMaskSet lms1, LatticeMaskSet lms2)");
        LatticeMaskSet res(rectangular_hull(lms1.block(),lms2.block()));
        res.adjoin(LatticeCellListSet(lms1));
        res.adjoin(LatticeCellListSet(lms2));
        return res;
      }
    }

    LatticeMaskSet
    difference(const LatticeMaskSet& lms1, const LatticeMaskSet& lms2) 
    {
      if(lms1.block()==lms2.block()) { 
        return LatticeMaskSet(lms1.block(),lms1.mask() - lms2.mask(), !lms1.bounded() - !lms2.bounded());
      } else {
        ARIADNE_CHECK_BOUNDED(lms1,"LatticeMaskSet difference(LatticeMaskSet lms1, LatticeMaskSet lms2)");
        ARIADNE_CHECK_BOUNDED(lms2,"LatticeMaskSet difference(LatticeMaskSet lms1, LatticeMaskSet lms2)");
        LatticeMaskSet res(lms1.block());
        res.adjoin(difference(LatticeCellListSet(lms1),lms2));
        return res;
      }
    }

    LatticeCellListSet
    difference(const LatticeCellListSet& lcls, const LatticeMaskSet& lms) 
    {
      LatticeCellListSet result(lcls.dimension());
      for(LatticeCellListSet::const_iterator citer=lcls.begin(); citer!=lcls.end(); ++citer) {
        if(!subset(*citer,lms)) {
          result.adjoin(*citer);
        }
      }
      return result;
    }

    bool 
    disjoint(const LatticeBlock& lb1, const LatticeBlock& lb2) 
    {
      for(dimension_type i=0; i!=lb1.dimension(); ++i) {
        if(lb1.upper_bound(i)<lb2.lower_bound(i)
            || lb1.upper_bound(i)<lb2.lower_bound(i))
        {
          return true;
        }
      }
      return false;
    }
    
    bool 
    disjoint(const LatticeBlock& lb, const LatticeMaskSet& lms) 
    {
      return !overlap(lb.neighbourhood(),lms);
    }
    
    bool 
    disjoint(const LatticeMaskSet& lms1, const LatticeMaskSet& lms2) 
    {
      return !overlap(lms1.neighbourhood(),lms2);
    }
    
    bool 
    overlap(const LatticeBlock& lb1, const LatticeBlock& lb2) 
    {
      for(dimension_type i=0; i!=lb1.dimension(); ++i) {
        if(lb1.upper_bound(i)<=lb2.lower_bound(i)
            || lb1.upper_bound(i)<=lb2.lower_bound(i))
        {
          return false;
        }
      }
      return true;
    }
    

    bool 
    overlap(const LatticeBlock& lb, const LatticeMaskSet& lms) 
    {
      LatticeBlock rlb=regular_intersection(lb,lms.block());
      if(rlb.empty()) {
        return false;
      }
      for(LatticeBlock::const_iterator biter=rlb.begin(); biter!=rlb.end(); ++biter) {
        if(subset(*biter,lms)) {
          return true;
        }
      }
      return false;
    }
    
    bool 
    overlap(const LatticeCellListSet& lcls, const LatticeMaskSet& lms) 
    {
      for(LatticeCellListSet::const_iterator citer=lcls.begin(); citer!=lcls.end(); ++citer) {
        if(subset(*citer,lms)) {
          return true;
        }
      }
      return false;
    }
    
    bool
    overlap(const LatticeMaskSet& lms1, const LatticeMaskSet& lms2) 
    {
      if(lms1.block()==lms2.block()) {
        BooleanArray::const_iterator iter1=lms1.mask().begin();
        BooleanArray::const_iterator iter2=lms2.mask().begin();
        BooleanArray::const_iterator end1=lms1.mask().end();
        while(iter1!=end1) {
          if(*iter1 & *iter2) {
            return true;
          }
          ++iter1;
          ++iter2;
        }
        return false;
      } else {
        for(LatticeMaskSet::const_iterator iter1=lms1.begin(); iter1!=lms1.end(); ++iter1) {
          if(subset(*iter1,lms2)) {
            return true;
          }
        }
        return false;
      }
    }
    

    bool
    subset(const LatticeCell& lc, const LatticeBlock& lb)
    {
      return subset(LatticeBlock(lc),lb);
    }
     

    bool 
    subset(const LatticeBlock& lb1, const LatticeBlock& lb2) 
    {
      if(lb1.empty()) { 
        return true; 
      }
      for(dimension_type i=0; i!=lb1.dimension(); ++i) {
        if(lb1.lower_bound(i)<lb2.lower_bound(i)
            || lb1.upper_bound(i)>lb2.upper_bound(i))
        {
          return false;
        }
      }
      return true;
    }

    bool
    subset(const LatticeCell& lc, const LatticeCellListSet& lcls)
    {
      return std::find(lcls.begin(),lcls.end(),lc)!=lcls.end();
    }
     
    bool 
    subset(const LatticeCellListSet& lcls, const LatticeBlock& lb) 
    {
      for(LatticeCellListSet::const_iterator citer=lcls.begin(); citer!=lcls.end(); ++citer) {
        if(!subset(*citer,lb)) {
          return false;
        }
      }
      return true;
    }

    bool 
    subset(const LatticeMaskSet& lms, const LatticeBlock& lb) 
    {
      if(!lms.bounded()) { return false; }
      for(LatticeMaskSet::const_iterator citer=lms.begin(); citer!=lms.end(); ++citer) {
        if(!subset(*citer,lb)) {
          return false;
        }
      }
      return true;
    }

    bool 
    subset(const LatticeCell& lc, const LatticeMaskSet& lms) 
    {
      if(!subset(lc,lms.block())) {
        return !lms.bounded();
      }
      return lms.mask()[lms.index(lc)];
    }
    
    bool 
    subset(const LatticeBlock& lb, const LatticeMaskSet& lms) 
    {
      if(lb.empty()) {
        return true;
      }
      if(!subset(lb,lms.block())) {
        return !lms.bounded();
      }
      for(LatticeBlock::const_iterator biter=lb.begin(); biter!=lb.end(); ++biter) {
        if(!subset(*biter,lms)) {
          return false;
        }
      }
      return true;
    }
    
    bool 
    subset(const LatticeCellListSet& lc, const LatticeMaskSet& lms) 
    {
      for(LatticeCellListSet::const_iterator citer=lc.begin(); citer!=lc.end(); ++citer) {
        if(!subset(*citer,lms)) {
          return false;
        }
      }
      return true;
    }
    
    bool 
    subset(const LatticeMaskSet& lms1, const LatticeMaskSet& lms2) {
      //FIXME: possible bug if lms1 and lms2 unbounded
      if(!lms1.bounded() && lms2.bounded()) { 
        return false;
      }
      if(lms1.block() == lms2.block()) {
        return lms1.mask() <= lms2.mask();
      } else {
        for(LatticeMaskSet::const_iterator iter1=lms1.begin(); iter1!=lms1.end(); ++iter1) {
          if(!subset(*iter1,lms2)) {
            return false;
          }
        }
        return true;
      }
    }
      

    LatticeMaskSet LatticeMaskSet::neighbourhood() const {
      LatticeMaskSet result=*this;
      const IndexArray& lower_bound=result.block().lower_corner();
      const IndexArray& upper_bound=result.block().upper_corner();
      for(LatticeMaskSet::const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        LatticeCell cell=*iter;
        IndexArray lower(this->dimension());
        IndexArray upper(this->dimension());
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          lower[i]=std::max(cell.lower_bound(i)-1,lower_bound[i]);
          upper[i]=std::min(cell.upper_bound(i)+1,upper_bound[i]);
        }
        LatticeBlock block(lower,upper);
        result.adjoin(block);
      }
      return result;
    }

    LatticeMaskSet LatticeMaskSet::adjoining() const {
      LatticeMaskSet result=*this;
      const IndexArray& lower_bound=result.block().lower_corner();
      const IndexArray& upper_bound=result.block().upper_corner();
      for(LatticeMaskSet::const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        LatticeCell cell=*iter;
        IndexArray lower=(cell.lower_corner());
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

    std::istream& 
    operator>>(std::istream& is, LatticePoint& lpt)
    {
      std::vector<int> v;
      read_sequence(is, v, '(', ')');
      lpt=LatticePoint(v.size());
      for(dimension_type i=0; i!=lpt.dimension(); ++i) {
        lpt[i]=v[i];
      }
      return is;
    }
    
    std::istream& 
    operator>>(std::istream& is, LatticeCell& lc)
    {
      /* Representation as a literal (l1,l2,...,ln) */
      std::vector<int> v;
      read_sequence(is, v, '(', ')');
      IndexArray l(v.size());
      for(size_type i=0; i!=v.size(); ++i) {
        l[i]=v[i];
      }
      lc=LatticeCell(l);
      return is;
    }
    
    
    std::istream& 
    operator>>(std::istream& is, LatticeBlock& r)
    {
      char c;
      is >> c;
      is.putback(c);
      if(c=='[') {
        /* Representation as a literal [a1,b1]x[a2,b2]x...x[an,bn] */
        std::vector< int > v;
        c='x';
        while(c=='x') {
          char cl,cm,cr;
          int l,u;
          is >> cl >> l >> cm >> u >> cr;
          if(cl!='[' || (cm!=',' && cm!=';') || cr!=']') {
            ARIADNE_THROW(InvalidInput,"istream& LatticeBlock::read(istream&)","");
          }
          v.push_back(l); v.push_back(u);
          c=' ';
          while( is && c==' ') {
            is >> c;
          }
        }
        if(is) {
          is.putback(c);
        }
        
        IndexArray l(v.size()/2);
        IndexArray u(v.size()/2);
        for(size_type i=0; i!=v.size()/2; ++i) {
          l[i]=v[2*i];
          u[i]=v[2*i+1];
        }
        r=LatticeBlock(l,u);
      }
      else {
        /* representation as lower and upper corners */
        /* FIXME */
        ARIADNE_THROW(InvalidInput,"istream& LatticeBlock::read(istream&)","");
      }
      return is;
    }
    
    std::ostream& 
    operator<<(std::ostream& os, const LatticePoint& lpt) 
    {
      os << "(" << lpt[0];
      for(dimension_type i=1; i!=lpt.dimension(); ++i) {
        os << "," << lpt[i];
      }
      return os << ")";
    }
        
    std::ostream& 
    operator<<(std::ostream& os, const LatticeCell& lc) 
    {
      if(lc.dimension()==0) {
        return os<<"Empty";
      }
      os << lc[0];
      for(dimension_type i=1; i!=lc.dimension(); ++i) {
        os << "x" << lc[i];
      }
      return os;
    }
        
    std::ostream& 
    operator<<(std::ostream& os, const LatticeBlock& lb) 
    {
      if(lb.empty() || lb.dimension()==0) {
        return os<<"Empty";
      }
      os << lb[0];
      for(dimension_type i=1; i!=lb.dimension(); ++i) {
        os << "x" << lb[i];
      }
      return os;
    }
    
    std::ostream& 
    operator<<(std::ostream& os, const LatticeMaskSet& lms) 
    {
      return os << "LatticeMaskSet(\n  block=" << lms.block() << "\n  mask=" << lms.mask() << "\n)\n";
    }
        
    std::ostream& 
    operator<<(std::ostream& os, const LatticeCellListSet& lcls) 
    {
      return write_sequence(os,lcls.begin(),lcls.end(),'[',']',',');
    }
    
    
    
    class chompfstream : public std::ofstream { };
    
    chompfstream& 
    operator<<(chompfstream& cfs, const LatticeMaskSet& lms) 
    {
      std::ostream& os=cfs;
      dimension_type d=lms.dimension();
      for(LatticeMaskSet::const_iterator iter=lms.begin(); iter!=lms.end(); ++iter) {
        LatticeCell lc=*iter;
        for(size_type i=0; i!=d; ++i) {
          uint k=lc.lower_bound(i);
          os << k << " ";
        }
        os << "\n";
      }
      return cfs;
    }
    
    
   
} // namespace Ariadne

