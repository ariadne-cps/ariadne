/***************************************************************************
 *            dynamics/ltl_formula.hpp
 *
 *  Copyright  20  Pieter Collins, Luca Geretti
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

/*! \file dynamics/ltl_formula.hpp
 *  \brief Linear Temporal Logic Formula.
 */

#ifndef ARIADNE_LTL_FORMULA_HPP
#define ARIADNE_LTL_FORMULA_HPP


#include "../utility/variant.hpp"
#include "../symbolic/identifier.hpp"
#include "../symbolic/templates.hpp"

namespace Ariadne {

struct AndOp;
struct OrOp;
struct NotOp;

struct Next { };
struct Always { };
struct Eventually { };
struct Until { };
struct Release { };

struct ForAll { };
struct ThereExists { };

class StateFormula;

class AtomicProposition;
class PathFormula;
class StateFormula;

class LTLOperators {
    friend PathFormula operator!(PathFormula const& phi);
    friend PathFormula operator&&(PathFormula const& phi1, PathFormula const& phi2);
    friend PathFormula operator||(PathFormula const& phi1, PathFormula const& phi2);
    friend PathFormula next(PathFormula const& phi);
    friend PathFormula until(PathFormula const& phi1, PathFormula const& phi2);
};

class AtomicProposition
    : public LTLOperators
{
    Identifier _name;
  public:
    AtomicProposition(Identifier name) : _name(name) { }
    Identifier name() const;
};

template<class OP, class... ARGS> Symbolic<OP,ARGS...> make_symbolic(OP op, ARGS... args) {
    return Symbolic<OP,ARGS...>(op,args...); }


class StateFormula {
    typedef Variant<True,AtomicProposition,
        Symbolic<NotOp,StateFormula>,Symbolic<AndOp,StateFormula,StateFormula>,Symbolic<OrOp,StateFormula,StateFormula>,
        Symbolic<ForAll,PathFormula>,Symbolic<ThereExists,PathFormula>> VariantType;
    SharedPointer<VariantType> _node;
  private:
    explicit StateFormula(SharedPointer<VariantType> const& ptr) : _node(ptr) { };
  public:
    StateFormula(True const& p)
        : _node(std::make_shared<VariantType>(p)) { }
    StateFormula(False p)
        : _node(std::make_shared<VariantType>(make_symbolic(NotOp(),StateFormula(True())))) { }
    StateFormula(AtomicProposition const& p)
        : _node(std::make_shared<VariantType>(p)) { }
    friend StateFormula operator!(StateFormula const& phi) {
        return StateFormula(std::make_shared<VariantType>(make_symbolic(NotOp(),phi))); }
    friend StateFormula operator&&(StateFormula const& phi1, StateFormula const& phi2) {
        return StateFormula(std::make_shared<VariantType>(make_symbolic(AndOp(),phi1,phi2))); }
    friend StateFormula operator||(StateFormula const& phi1, StateFormula const& phi2) {
        return StateFormula(std::make_shared<VariantType>(make_symbolic(OrOp(),phi1,phi2))); }
    friend StateFormula all(PathFormula const& psi);
    friend StateFormula exists(PathFormula const& psi);

};

class PathFormula {
    typedef Variant<StateFormula,
        Symbolic<NotOp,PathFormula>,Symbolic<AndOp,PathFormula,PathFormula>,Symbolic<OrOp,PathFormula,PathFormula>,
        Symbolic<Next,PathFormula>, Symbolic<Always,PathFormula>, Symbolic<Eventually,PathFormula>,
        Symbolic<Until,PathFormula,PathFormula>,Symbolic<Release,PathFormula,PathFormula>> VariantType;
    SharedPointer<VariantType> _node;
  private:
    explicit PathFormula(SharedPointer<VariantType> const& ptr) : _node(ptr) { };
  public:
    PathFormula(True const& p)
        : PathFormula(StateFormula(p)) { }
    PathFormula(False const& p)
        : PathFormula(StateFormula(p)) { }
    PathFormula(AtomicProposition const& p)
        : PathFormula(StateFormula(p)) { }
    PathFormula(StateFormula const& psi)
        : _node(std::make_shared<VariantType>(psi)) { }
    friend PathFormula operator!(PathFormula const& phi) {
        return PathFormula(std::make_shared<VariantType>(make_symbolic(NotOp(),phi))); }
    friend PathFormula operator&&(PathFormula const& phi1, PathFormula const& phi2) {
        return PathFormula(std::make_shared<VariantType>(make_symbolic(AndOp(),phi1,phi2))); }
    friend PathFormula operator||(PathFormula const& phi1, PathFormula const& phi2) {
        return PathFormula(std::make_shared<VariantType>(make_symbolic(OrOp(),phi1,phi2))); }
    friend PathFormula next(PathFormula const& phi) {
        return PathFormula(std::make_shared<VariantType>(make_symbolic(Next(),phi))); }
    friend PathFormula eventually(PathFormula const& phi) {
        // PathFormula tru=True(); return PathFormula(std::make_shared<VariantType>(make_symbolic(Until(),tru,phi))); }
        return PathFormula(std::make_shared<VariantType>(make_symbolic(Eventually(),phi))); }
    friend PathFormula always(PathFormula const& phi) {
        return PathFormula(std::make_shared<VariantType>(make_symbolic(Always(),phi))); }
    friend PathFormula until(PathFormula const& phi1, PathFormula const& phi2) {
        return PathFormula(std::make_shared<VariantType>(make_symbolic(Until(),phi1,phi2))); }
    friend PathFormula release(PathFormula const& phi1, PathFormula const& phi2) {
        return PathFormula(std::make_shared<VariantType>(make_symbolic(Until(),phi1,phi2))); }

    friend StateFormula all(PathFormula const& psi) {
        return StateFormula(std::make_shared<StateFormula::VariantType>(make_symbolic(ForAll(),psi))); }
    friend StateFormula exists(PathFormula const& psi) {
        return StateFormula(std::make_shared<StateFormula::VariantType>(make_symbolic(ForAll(),psi))); }
};


} // namespace Ariadne

#endif // ARIADNE_LTL_FORMULA_HPP
