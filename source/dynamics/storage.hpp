/***************************************************************************
 *            dynamics/storage.hpp
 *
 *  Copyright  2011-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file dynamics/storage.hpp
 *  \brief Paved sets for continuous systems
 */

#ifndef ARIADNE_STORAGE_HPP
#define ARIADNE_STORAGE_HPP

#include <iosfwd>
#include "../geometry/set_interface.hpp"
#include "../geometry/paving_interface.hpp"
#include "../geometry/grid_paving.hpp"
#include "../output/graphics_interface.hpp"

namespace Ariadne {

//! \ingroup DynamicsModule
//! \brief A storage for part of the reachable or evolved set of a dynamical system.
class Storage
    : public DrawableInterface
{
    GridTreePaving _state_set;
    EffectiveVectorMultivariateFunction _auxiliary_mapping;
  public:
    typedef Grid GridType;
    Storage(Grid const& g) : _state_set(g), _auxiliary_mapping(0,g.dimension()) { }
    Storage(GridTreePaving const& gtp) : _state_set(gtp), _auxiliary_mapping(0,gtp.dimension()) { }
    Storage(Grid const& g, EffectiveVectorMultivariateFunction const& aux)
        : _state_set(g), _auxiliary_mapping(aux) { }
    Storage(GridTreePaving const& gtp, EffectiveVectorMultivariateFunction const& aux)
        : _state_set(gtp), _auxiliary_mapping(aux) { }

    template<class SYS> Void set_system(SYS const&) { }
    Void set_auxiliary(EffectiveVectorMultivariateFunction const& aux) { this->_auxiliary_mapping=aux; }
    EffectiveVectorMultivariateFunction auxiliary() const { return this->_auxiliary_mapping; }

    GridTreePaving& state_set() { return this->_state_set; }
    GridTreePaving const& state_set() const { return this->_state_set; }
    GridTreePaving state_auxiliary_set() const {
        Grid auxiliary_grid(this->_auxiliary_mapping.result_size());
        return outer_skew_product(this->_state_set, auxiliary_grid, this->_auxiliary_mapping); }

    Bool is_empty() const { return this->_state_set.is_empty(); }
    SizeType size() const { return this->_state_set.size(); }
    GridTreePaving::ConstIterator begin() const { return this->_state_set.begin(); }
    GridTreePaving::ConstIterator end() const { return this->_state_set.end(); }

    Grid grid() const { return this->_state_set.grid(); }

    Void clear() { this->_state_set.clear(); }
    Void adjoin(Storage const& other) { return this->_state_set.adjoin(other._state_set); }
    Void restrict(Storage const& other) { return this->_state_set.restrict(other._state_set); }
    Void remove(Storage const& other) { return this->_state_set.remove(other._state_set); }

    Void adjoin_lower_approximation(OvertSetInterface const& set, Nat extent, Nat fineness) {
        this->_state_set.adjoin_lower_approximation(set,extent,fineness); }
    Void adjoin_outer_approximation(CompactSetInterface const& set, Nat fineness) {
        this->_state_set.adjoin_outer_approximation(set,fineness); }
    Void adjoin_outer_approximation(ValidatedCompactSetInterface const& set, Nat fineness) {
        this->_state_set.adjoin_outer_approximation(set,fineness); }
    Void adjoin_outer_approximation(ExactBoxType const& bx, Nat fineness) {
        this->_state_set.adjoin_outer_approximation(bx,fineness); }

    Void restrict_to_extent(Nat level) { this->_state_set.restrict_to_extent(level); }
    Void mince(Nat fineness) { this->_state_set.mince(fineness); }
    Void recombine() { this->_state_set.recombine(); }

    friend Bool subset(Storage const& set1, Storage const& set2) {
        //ARIADNE_PRECONDITION(same(set1._auxiliary_mapping,set2._auxiliary_mapping));
        return subset(set1._state_set,set2._state_set); }
    friend Bool intersect(Storage const& set1, Storage const& set2) {
        return intersect(set1._state_set,set2._state_set); }

    virtual Storage* clone() const override { return new Storage(*this); }
    virtual DimensionType dimension() const override { return this->_state_set.dimension(); }
    virtual Void draw(CanvasInterface& c, const Projection2d& p) const override {
        std::cerr<<"Storage::draw(Canvas, Projection2d): dimension="<<this->dimension()<<", projection="<<p<<"\n";
        GridTreePaving state_aux_set=this->state_auxiliary_set();
        std::cerr<<"state_auxiliary_set.size()="<<state_aux_set.size()<<"\n";
        return this->state_auxiliary_set().draw(c,p); }

    friend OutputStream& operator<<(OutputStream& os, Storage const& set) {
        return os << "Storage(" << set._state_set << "," << set._auxiliary_mapping << ")"; }
};

} //namespace Ariadne

#endif /* ARIADNE_STORAGE_HPP */
