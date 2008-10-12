#ifndef ARIADNE_HYBRID_SET_H
#define ARIADNE_HYBRID_SET_H

#include <map>
#include "macros.h"
#include "grid.h"

namespace Ariadne {

typedef uint DiscreteState;

class GridSet { 
  template<class S> void adjoin(const S&); 
  template<class S> void adjoin_outer_approximation(const S&); 
};

class HybridGridSet 
  : public std::map<DiscreteState,GridSet>
{
  bool has_location(DiscreteState q) const {
    return this->find(q)!=this->end(); }

  template<class S> void adjoin(DiscreteState q, const S& s) {
    this->operator[](q).adjoin(s); }

  template<class S> void adjoin_outer_approximation(DiscreteState q, const S& s) {
    this->operator[](q).adjoin_outer_approximation(s); }

  const GridSet& operator[](DiscreteState q) const {
    ARIADNE_ASSERT(this->find(q)!=this->end());
    return const_cast<HybridGridSet*>(this)->operator[](q);
  }

};

} // namespace Ariadne

#endif
