/***************************************************************************
 *            approximate_set.h
 *
 *  Copyright 2008  Pieter Collins
 * 
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
 
/*! \file approximate_set.h
 *  \brief Approximations to for open, closed, overt and compact subsets of Euclidean space.
 */

#ifndef ARIADNE_APPROXIMATE_SET_H
#define ARIADNE_APPROXIMATE_SET_H

#include <iosfwd>

#include "tribool.h"
#include "set_interface.h"

#include "grid_set.h"
#include "hybrid_set.h"

namespace Ariadne {

//! \brief An inner approximation to an open set.
class InnerApproximation 
    : public OpenSetInterface
{
  public:
    InnerApproximation(const Grid& g) : _concrete_set(g) { }
    InnerApproximation(const GridTreeSet& gts) : _concrete_set(gts) { }
    InnerApproximation* clone() const { return new InnerApproximation(*this); }
    uint dimension() const { return this->_concrete_set.dimension(); }
    tribool overlaps(const Box& bx) const { return !this->_concrete_set.disjoint(bx) || indeterminate; }
    tribool covers(const Box& bx) const { return this->_concrete_set.superset(bx) || indeterminate; }
    std::ostream& write(std::ostream& os) const { return os << this->_concrete_set; }
    friend std::ostream& operator<<(std::ostream& os, const InnerApproximation& ia) { return ia.write(os); };
  private:
    GridTreeSet _concrete_set;
};



//! \brief A lower approximation to an overt set \a S, consisting of a union of boxes, each of which intersects \a S.
class LowerApproximation
    : public OvertSetInterface
{
  public:
    LowerApproximation(const Grid& g) : _concrete_set(g) { }
    LowerApproximation(const GridTreeSet& gts) : _concrete_set(gts) { }
    LowerApproximation* clone() const { return new LowerApproximation(*this); }
    uint dimension() const { return this->_concrete_set.dimension(); }
    tribool overlaps(const Box& bx) const { return !this->_concrete_set.disjoint(bx) || indeterminate; }
    std::ostream& write(std::ostream& os) const { return os << this->_concrete_set; }
    friend std::ostream& operator<<(std::ostream& os, const LowerApproximation& la) { return la.write(os); };
  private:
    GridTreeSet _concrete_set;
};


//! \brief An outer approximation to an overt set \a S, consisting of a union of boxes covering \a S.
class OuterApproximation 
    : public virtual CompactSetInterface
{
  public:
    OuterApproximation(const Grid& g) : _concrete_set(g) { }
    OuterApproximation(const GridTreeSet& gts) : _concrete_set(gts) { }
    OuterApproximation* clone() const { return new OuterApproximation(*this); }
    uint dimension() const { return this->_concrete_set.dimension(); }
    tribool disjoint(const Box& bx) const { return !this->_concrete_set.overlaps(bx) || indeterminate; }
    tribool inside(const Box& bx) const { return this->_concrete_set.subset(bx) || indeterminate; }
    Box bounding_box() const { return this->_concrete_set.bounding_box(); }
    tribool empty() const { return this->_concrete_set.empty(); }
    std::ostream& write(std::ostream& os) const { return os << this->_concrete_set; }
    friend std::ostream& operator<<(std::ostream& os, const OuterApproximation& oa) { return oa.write(os); }
  private:
    GridTreeSet _concrete_set;
};

//! \brief An approximation to a located set \a S, consisting of a union of boxes, each of which intersects \a S and whose union covers \a S.
class MetricApproximation 
    : public virtual LocatedSetInterface
{
  public:
    MetricApproximation(const Grid& g) : _concrete_set(g) { }
    MetricApproximation(const GridTreeSet& gts) : _concrete_set(gts) { }
    MetricApproximation* clone() const { return new MetricApproximation(*this); }
    uint dimension() const { return this->_concrete_set.dimension(); }
    tribool overlaps(const Box& bx) const { return !this->_concrete_set.disjoint(bx) || indeterminate; }
    tribool disjoint(const Box& bx) const { return !this->_concrete_set.overlaps(bx) || indeterminate; }
    tribool inside(const Box& bx) const { return this->_concrete_set.subset(bx) || indeterminate; }
    Box bounding_box() const { return this->_concrete_set.bounding_box(); }
    tribool empty() const { return this->_concrete_set.empty(); }
    std::ostream& write(std::ostream& os) const { return os << this->_concrete_set; }
    friend std::ostream& operator<<(std::ostream& os, const MetricApproximation& oa) { return oa.write(os); }
  private:
    GridTreeSet _concrete_set;
};




class HybridInnerApproximation 
    : public HybridOpenSetInterface
{
  public:
    HybridInnerApproximation() : _concrete_set() { }
    HybridInnerApproximation* clone() const { return new HybridInnerApproximation(*this); }
    HybridSpace space() const { return HybridSpace(this->_concrete_set); }
    InnerApproximation const& operator[](DiscreteState q) const { return this->_concrete_set.find(q)->second; }
    tribool overlaps(const HybridBox& bx) const { return this->_concrete_set.find(bx.first)->second.overlaps(bx.second); }
    tribool covers(const HybridBox& bx) const { return this->_concrete_set.find(bx.first)->second.covers(bx.second); }
    std::ostream& write(std::ostream& os) const { return os << this->_concrete_set; }
  private:
    std::map<DiscreteState,InnerApproximation> _concrete_set;
};



class HybridLowerApproximation
    : public HybridOvertSetInterface
{
  public:
    HybridLowerApproximation() : _concrete_set() { }
    HybridLowerApproximation* clone() const { return new HybridLowerApproximation(*this); }
    HybridSpace space() const { return HybridSpace(this->_concrete_set); }
    LowerApproximation const& operator[](DiscreteState q) const { return this->_concrete_set.find(q)->second; }
    tribool overlaps(const HybridBox& bx) const { return this->_concrete_set.find(bx.first)->second.overlaps(bx.second); }
    std::ostream& write(std::ostream& os) const { return os << this->_concrete_set; }
  private:
    std::map<DiscreteState,LowerApproximation> _concrete_set;
};


//! \brief Interface for overt sets. 
class HybridOuterApproximation 
    : public HybridCompactSetInterface
{
  public:
    HybridOuterApproximation() : _concrete_set() { }
    HybridOuterApproximation* clone() const { return new HybridOuterApproximation(*this); }
    HybridSpace space() const { return HybridSpace(this->_concrete_set); }
    OuterApproximation const& operator[](DiscreteState q) const { return this->_concrete_set.find(q)->second; }
    tribool disjoint(const HybridBox& hbx) const { return this->_concrete_set.find(hbx.first)->second.disjoint(hbx.second); }
    std::ostream& write(std::ostream& os) const { return os << this->_concrete_set; }
    tribool inside(const HybridBoxes& hbxs) const;
    HybridBoxes bounding_box() const;
  private:
    std::map<DiscreteState,OuterApproximation> _concrete_set;
};

class HybridMetricApproximation
    : public virtual HybridLocatedSetInterface
{
  public:
    HybridMetricApproximation() : _concrete_set() { }
    HybridMetricApproximation* clone() const { return new HybridMetricApproximation(*this); }
    HybridSpace space() const { return HybridSpace(this->_concrete_set); }
    MetricApproximation const& operator[](DiscreteState q) const { return this->_concrete_set.find(q)->second; }
    tribool overlaps(const HybridBox& bx) const { return this->_concrete_set.find(bx.first)->second.overlaps(bx.second); }
    tribool disjoint(const HybridBox& hbx) const { return this->_concrete_set.find(hbx.first)->second.disjoint(hbx.second); }
    tribool inside(const HybridBoxes& hbxs) const { return reinterpret_cast<const HybridOuterApproximation&>(*this).inside(hbxs); }
    HybridBoxes bounding_box() const { return reinterpret_cast<const HybridOuterApproximation&>(*this).bounding_box(); };
    std::ostream& write(std::ostream& os) const { return os << this->_concrete_set; }
  private:
    std::map<DiscreteState,MetricApproximation> _concrete_set;
};

    
inline tribool HybridOuterApproximation::inside(const HybridBoxes& hbxs) const { 
    for( std::map<DiscreteState,OuterApproximation>::const_iterator loc_iter = this->_concrete_set.begin(); 
         loc_iter != this->_concrete_set.end(); ++loc_iter ) 
    {
        if( !loc_iter->second.empty() ) { 
            HybridBoxes::const_iterator hbxs_loc_iter = hbxs.find( loc_iter->first ); 
            if( hbxs_loc_iter != hbxs.end() && ! loc_iter->second.inside( hbxs_loc_iter->second ) ) { 
                return false; 
            }
        } 
    }
    return true;
}

inline HybridBoxes HybridOuterApproximation::bounding_box() const {  
    HybridBoxes result;
    for( std::map<DiscreteState,OuterApproximation>::const_iterator loc_iter = this->_concrete_set.begin(); 
         loc_iter != this->_concrete_set.end(); ++loc_iter ) 
    {
        if( !loc_iter->second.empty() ) {
            result.insert( std::make_pair( loc_iter->first, loc_iter->second.bounding_box() ) );
        }
    } 
    return result;
}


}

#endif // ARIADNE_APPROXIMATE_SET
