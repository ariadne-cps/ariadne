/***************************************************************************
 *            lattice_map.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it,  Pieter.Collins@cwi.nl
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

#include "debug.h"

#include "base/stlio.h"

#include "combinatoric/lattice_set.h"
#include "combinatoric/lattice_map.h"


Ariadne::dbgstream cdbg(std::cerr,0);

namespace Ariadne {
  namespace Combinatoric {
    
    void
    LatticeMultiMap::adjoin_to_image(const LatticeCell& lc, const LatticeCell& img)
    {
      typedef std::map<LatticeCell,LatticeCellListSet>::iterator map_iterator;
      map_iterator iter=this->_map.find(lc);
      if(iter==this->_map.end()) {
        this->_map.insert(std::make_pair(lc,LatticeCellListSet(this->_result_dimension)));
        iter=this->_map.find(lc);
      }
      iter->second.adjoin(img);
    }
    
    void
    LatticeMultiMap::adjoin_to_image(const LatticeCell& lc, const LatticeCellListSet& img)
    {
      //std::cerr << "LatticeMultiMap::adjoin_to_image(const LatticeCell& lc, const LatticeCellListSet& img)" << std::endl;
      //std::cerr << "result_dimension=" << this->_result_dimension << " lc=" << lc << " img=" << img << std::endl;
      typedef std::map<LatticeCell,LatticeCellListSet>::iterator map_iterator;
      map_iterator iter=this->_map.find(lc);
      if(iter==this->_map.end()) {
        this->_map.insert(std::make_pair(lc,LatticeCellListSet(this->_result_dimension)));
        //std::cerr << *this << std::endl;
        iter=this->_map.find(lc);
        //std::cerr << *iter << std::endl;
        assert(iter!=this->_map.end());
      }
      //std::cerr << iter->first << " " << iter->second.dimension() << " " << iter->second << std::endl;
      iter->second.adjoin(img);
    }
    

    LatticeCellListSet
    LatticeMultiMap::image(const LatticeCell& lc) const
    {
      cdbg << "LatticeMap::image(const LatticeCell& lc) const\n";
      typedef std::map<LatticeCell,LatticeCellListSet>::iterator map_iterator;
      map_iterator iter=this->_map.find(lc);
      if(iter==this->_map.end()) {
        LatticeCellListSet img(this->_result_dimension);
        this->_map.insert(std::make_pair(lc,img));
        iter=this->_map.find(lc);
        assert(iter!=this->_map.end());
      }
      return iter->second;
    }
    
    LatticeCellListSet 
    LatticeMultiMap::image(const LatticeMaskSet& lms) const 
    {
      LatticeCellListSet result(this->_result_dimension);
      for(LatticeMaskSet::const_iterator cell_iter=lms.begin(); cell_iter!=lms.end(); ++cell_iter) {
        result.adjoin(this->apply(*cell_iter));
      }
      return result;
    }

    LatticeCellListSet
    LatticeMultiMap::apply(const LatticeCell& lc) const
    {
      return this->image(lc);
    }

    LatticeCellListSet
    LatticeMultiMap::weak_preimage(const LatticeMaskSet& lms) const
    {
      cdbg << "LatticeMap::weak_preimage(const LatticeCell& lc) const\n";
      typedef std::map<LatticeCell,LatticeCellListSet>::iterator map_iterator;
      LatticeCellListSet result(this->argument_dimension());
      for(map_iterator iter=this->_map.begin(); iter!=this->_map.end(); ++iter) {
        LatticeCellListSet& image=iter->second;
        if(overlap(image,lms)) { 
          result.adjoin(iter->first);
        }
      }
      return result;
    }
    
    LatticeCellListSet
    LatticeMultiMap::strong_preimage(const LatticeMaskSet& lms) const
    {
      cdbg << "LatticeMap::strong_preimage(const LatticeCell& lc) const\n";
      typedef std::map<LatticeCell,LatticeCellListSet>::iterator map_iterator;
      LatticeCellListSet result(this->argument_dimension());
      for(map_iterator iter=this->_map.begin(); iter!=this->_map.end(); ++iter) {
        LatticeCellListSet& image=iter->second;
        if(subset(image,lms)) { 
          result.adjoin(iter->first);
        }
      }
      return result;
    }
    
    LatticeCellListSet
    LatticeMultiMap::operator() (const LatticeCell& lc) const
    {
      cdbg << "LatticeMap::operator() (const LatticeCell& lc) const\n";
      return this->apply(lc);
    }
    
    LatticeCellListSet 
    LatticeMultiMap::operator() (const LatticeBlock& lr) const 
    {
      cdbg << "LatticeMap::operator() (const LatticeBlock& lr) const\n";
      LatticeCellListSet result(this->_result_dimension);
      for(LatticeBlock::const_iterator cell_iter=lr.begin(); cell_iter!=lr.end(); ++cell_iter) {
        result.adjoin(this->apply(*cell_iter));
      }
      return result;
    }
    
    LatticeCellListSet 
    LatticeMultiMap::operator() (const LatticeCellListSet& lcls) const 
    {
      cdbg << "LatticeMap::operator() (const LatticeCellListSet& lcls) const\n";
      LatticeCellListSet result(this->_result_dimension);
      for(LatticeCellListSet::const_iterator cell_iter=lcls.begin(); cell_iter!=lcls.end(); ++cell_iter) {
        result.adjoin(this->apply(*cell_iter));
      }
      cdbg << "Returning" << result << "\n";
      return result;
    }
    
    LatticeCellListSet 
    LatticeMultiMap::operator() (const LatticeBlockListSet& lrls) const 
    {
      cdbg << "LatticeMap::operator() (const LatticeBlockListSet& lrls) const\n";
      LatticeCellListSet result(this->_result_dimension);
      for(LatticeBlockListSet::const_iterator rect_iter=lrls.begin(); rect_iter!=lrls.end(); ++rect_iter) {
        result.adjoin(this->operator()(*rect_iter));
      }
      return result;
    }
    
    LatticeCellListSet 
    LatticeMultiMap::operator() (const LatticeMaskSet& lms) const 
    {
      LatticeCellListSet result(this->_result_dimension);
      for(LatticeMaskSet::const_iterator cell_iter=lms.begin(); cell_iter!=lms.end(); ++cell_iter) {
        result.adjoin(this->apply(*cell_iter));
      }
      return result;
    }
    
    LatticeMultiMap
    LatticeMultiMap::inverse() const 
    {
      LatticeMultiMap result(this->_result_dimension, this->_argument_dimension);
      for(LatticeMultiMap::const_iterator arg_iter=this->begin(); arg_iter!=this->end(); ++arg_iter) {
        const LatticeCellListSet& image=arg_iter->second;
        for(LatticeCellListSet::const_iterator im_iter=image.begin(); im_iter!=image.end(); ++im_iter) {
          result.adjoin_to_image(*im_iter,arg_iter->first);
        }
      }
      return result;
    }
    
    LatticeMultiMap::const_iterator
    LatticeMultiMap::begin() const 
    {
      return _map.begin();
    }
    
    LatticeMultiMap::const_iterator
    LatticeMultiMap::end() const 
    {
      return _map.end();
    }

    
    void
    LatticeSystem::set_control_values(const LatticeCell& lc, const LatticeCellListSet& lcls) 
    {
      assert(lc.dimension()==this->space_dimension());
      assert(lcls.dimension()==this->input_dimension());
      this->_control_map.insert(std::make_pair(lc,lcls));
    }
    
    void
    LatticeSystem::set_noise_values(const LatticeCell& splc, const LatticeCell& inlc,const LatticeCellListSet& lcls) 
    {
      assert(splc.dimension()==this->space_dimension());
      assert(inlc.dimension()==this->input_dimension());
      assert(lcls.dimension()==this->space_dimension());
      LatticeCell lc(this->space_dimension()+this->input_dimension());
      for(dimension_type i=0; i!=this->space_dimension(); ++i) {
        lc[i]=splc[i];
      }
      for(dimension_type i=0; i!=this->input_dimension(); ++i) {
        lc[i+this->space_dimension()]=inlc[i];
      }
      this->set_noise_values(lc,lcls);
    }
    
    void
    LatticeSystem::set_noise_values(const LatticeCell& lc, const LatticeCellListSet& lcls) 
    {
      assert(lc.dimension()==this->space_dimension()+this->input_dimension());
      assert(lcls.dimension()==this->input_dimension());
      this->_noise_map.insert(std::make_pair(lc,lcls));
    }
    
    LatticeCellListSet
    LatticeSystem::preimage(const LatticeMaskSet& lms) const
    {
      LatticeCellListSet result(this->space_dimension());
      
      assert(lms.dimension()==this->space_dimension());
      typedef std::map<LatticeCell,LatticeCellListSet>::iterator map_iterator;
      
      LatticeCell space(IndexArray(this->space_dimension()));
      LatticeCell control(IndexArray(this->input_dimension()));
      
      for(LatticeMultiMap::const_iterator nz_iter=this->_noise_map.begin();
          nz_iter!=this->_noise_map.end(); ++nz_iter)
      {
        if(subset(nz_iter->second,lms)) {
          const LatticeCell& cell=nz_iter->first;
          for(dimension_type i=0; i!=this->space_dimension(); ++i) {
            space[i]=cell[i];
          }
          
          // If we've already found this cell, we can skip the next steps
          if(!subset(space,result)) {
            for(dimension_type i=0; i!=this->input_dimension(); ++i) {
              control[i]=cell[space_dimension()+i];
            }
            map_iterator cntrl_iter=this->_control_map.find(space);
            if(cntrl_iter!=this->_control_map.end()) {
              const LatticeCellListSet& possible_controls(cntrl_iter->second);
              if(subset(control,possible_controls)) {
              result.adjoin(space);
              }
            }
          } 
        }
      }
      return result;
    }
          
    
    
    std::ostream&
    operator<<(std::ostream& os, const LatticeMultiMap& m) 
    {
      return os << "LatticeMap(" << m._map << ")";
    }
    
  }
}
