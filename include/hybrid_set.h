/***************************************************************************
 *            hybrid_set.h
 *
 *  Mon Feb  7 19:57:13 2005
 *  Copyright  2005  Alberto Casagrande
 *  casagrande@dimi.uniud.it
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
 
#ifndef _HYBRID_SET_H
#define _HYBRID_SET_H


#include <discrete_location.h>
#include <denotable_set.h>

template <typename DS>
inline void _print_ds_(std::ostream &os,
                       const DS &p) {
  
  using namespace Ariadne::Geometry::IO_Operators;
  
  os << p;
}

namespace Ariadne {	
namespace HybridDefinitions {
	
template <typename LOC>
class LocationDenotableSet;

template <typename LDS>
class HybridDenotableSet;

template <typename LOC>
std::ostream& operator<<(std::ostream &os, 
                         const LocationDenotableSet<LOC> &A)
{
  os << "[" << std::endl
     << "Location=\""<< (A.location).name() << "\""<<std::endl
     << "Denotable Set = [" << std::endl;
  
  _print_ds_(os, A.set);
  
  os << std::endl << "] ]";
  
  return os;	
}

template <typename LDS>
std::ostream& operator<<(std::ostream &os, 
                         const HybridDenotableSet<LDS> &A)
{
  if (A.size() >0 ) {
    os << "LocationDenotableSet[0]=" << A[0];
  }
  for (size_t i=1; i< A.size(); i++) {
    os << std::endl <<"LocationDenotableSet["<<i<<"]=" 
       << A._location_set[i];
  }
  
  return os;	
}

template < typename LOC >
class LocationDenotableSet {
 public:
  typedef LOC DiscreteLocation;
  typedef typename DiscreteLocation::DenotableSet DenotableSet;
  typedef typename DenotableSet::BasicSet BasicSet;
  typedef typename BasicSet::State State;
  typedef typename State::Real Real;
  
 public:
  
  DiscreteLocation location;
  
  DenotableSet set;
  
  /*! \brief A denotable set constructor. */
  LocationDenotableSet(const DiscreteLocation &loc): 
    location(loc), set((loc.vector_field()).dimension()) {}
  
  /*! \brief A denotable set constructor. */
  LocationDenotableSet(const DiscreteLocation &loc, 
                       const BasicSet &A):
    
    location(loc), set(A) {}
  
  /*! \brief A denotable set constructor. */
  LocationDenotableSet(const DiscreteLocation &location,
                       const DenotableSet &A):
    location(loc), set(A) {}
  
  /*! \brief A denotable set constructor. */
  LocationDenotableSet(
                       const LocationDenotableSet< DiscreteLocation > &A):
    location(A.location), set(A.set) {}
  
  /*! \brief Copy the hybrid denotable set. */
  inline LocationDenotableSet<DiscreteLocation> &operator=(
                                                           const LocationDenotableSet<DiscreteLocation> &A) {
    
    this->location=A.location;
    this->set=A.set;
    
    return *this;
  }
  
  inline const std::string &name() const {
    return (this->location).name();
  }
  
  inline const DiscreteLocationID &id() const {
    return (this->location).id();
  }
  
  template <typename L>
  friend std::ostream& operator<<(std::ostream &os, 
                                  const LocationDenotableSet<L> &A);
  
};

template < typename LDS >
class HybridDenotableSet {
 public:
  typedef LDS LocationDenotableSet;
  typedef typename LocationDenotableSet::DiscreteLocation DiscreteLocation;
  typedef typename DiscreteLocation::DenotableSet DenotableSet;
  typedef typename DenotableSet::BasicSet BasicSet;
  typedef typename BasicSet::State State;
  typedef typename State::Real Real;
 private:
  
  inline bool _insert_unordered(const DiscreteLocation &location, 
                                const DenotableSet &A)	{
    
    for (size_t i=0; i<this->size(); i++) {
      
      if (((*this[i]).name())==location.name()) {
        
        ((*this[i]).set).inplace_union(A);
        
        return true;
      }
    }
    
    return false;
  }
  
  inline bool _insert_unordered(const DiscreteLocation &location, 
                                const BasicSet &A)	{
    
    for (size_t i=0; i<this->size(); i++) {
      
      if (((this->_location_set[i]).name())==location.name()) {
        ((this->_location_set[i]).set).inplace_union(A);
        return true;
      }
    }
    
    return false;
  }
  
  bool _ordered;
  
  std::vector<LocationDenotableSet> _location_set;
  
 public:
  /*! \brief A denotable set constructor. */
  HybridDenotableSet(): _ordered(false) {}
  
  inline void add_location_set(const DiscreteLocation &loc) {
    
    LocationDenotableSet set(loc);
    
    (this->_location_set).push_back(set);
  }
  
  inline void add_location_set(const DiscreteLocation &loc, 
                               const BasicSet &A) {
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    if (this->_ordered) {
      (this->_location_set[loc.id()].set).inplace_union(A);
      
      return;				
    }
    
    if (this->_insert_unordered(loc,A)) { 
      
#ifdef DEBUG
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
      
      return;
    }
    
    LocationDenotableSet set(loc,A);
    
    (this->_location_set).push_back(set);		
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
  }
  
  inline void add_location_set(const DiscreteLocation &location,
                               const DenotableSet &A) {
    
    if (this->_ordered) {
      (this->_location_set[loc.id()].set).inplace_union(A);
      
      return;				
    }
    
    if (this->_insert_unordered(loc,A)) return;
    
    LocationDenotableSet set(loc,A);
    
    (this->_location_set).push_back(set);				
    
  }
  
  inline void add_location_set(const LocationDenotableSet &A){
    
    if (this->_ordered) {
      (this->_location_set[A.id()].set).inplace_union(A.set);
      
      return;				
    }
    
    if (this->_insert_unordered(A.location,A.set)) return;
    
    (this->_location_set).push_back(A);
  }
  
  inline size_t size() const {	
    return (this->_location_set).size();
  }
  
  /*! \brief Copy the hybrid denotable set. */
  inline HybridDenotableSet< LDS > &operator=(
                                              const HybridDenotableSet< LDS > &A) {
    
    for (size_t i=0; i< A.size(); i++) {
      (this->_location_set).push_back(A._location_set[i]);
    }				
    
    return *this;
  }
  
  inline const LDS &operator[](
                               const size_t &index) const{
    
    if (index>=this->size()) {
      throw std::invalid_argument("Index is bigger than vector size.");	
    }
    
    return this->_location_set[index];
  }
  
  
  inline HybridDenotableSet< LDS > &
  delete_element(const size_t &index) {
    size_t i;
    
    std::vector<LocationDenotableSet> l_set=this->_location_set;
    
    if (index>=this->size()) {
      throw std::invalid_argument("Index is bigger than vector size.");	
    }
    
    for (i=0; i< index; i++) {
      (this->_location_set).push_back(l_set[i]);
    }
    
    for (i=index+1; i< l_set.size(); i++) {
      (this->_location_set).push_back(l_set[i]);
    }
    
    return *this;
  }
  
  inline HybridDenotableSet< LDS > &
  order_using(const std::vector<DiscreteLocation> &l_vector) {
    
    size_t i,j,new_size=l_vector.size();
    bool not_found=true;	
    
    std::vector<LocationDenotableSet> l_set;
    
    for (i=0; i<new_size; i++) {
      l_set.push_back(this->_location_set[0]);
    }
    
    for (i=0; i<this->size(); i++) {
      
      for (j=0; j<new_size; j++) {
        
        LocationDenotableSet &loc_ds=this->_location_set[i];
        
        if (l_vector[j].name()==(loc_ds.location).name()) {
          
          (loc_ds.location)._set_id(l_vector[j].id());
          
          l_set[l_vector[j].id()]=loc_ds;
          
          not_found=false;
        }			
      }
      if (not_found) {				
        throw std::invalid_argument("A discrete location is not present into the vector.");
      }
    }
    
    this->_location_set=l_set;
    
    this->_ordered=true;
    
    return *this;
  }
  
  inline bool ordered() const {
    return this->_ordered;
  }
  
  template <typename L>
  friend std::ostream& operator<<(std::ostream &os, 
                                  const HybridDenotableSet<L> &A);
};

}
}

#endif /* _HYBRID_SET_H */
