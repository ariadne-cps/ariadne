/***************************************************************************
 *            hybrid_system.cc
 *
 *  Copyright  2009  Pieter Collins
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

#include <map>

#include "macros.h"
#include "stlio.h"
#include "formula.h"
#include "function_interface.h"
#include "hybrid_time.h"
#include "hybrid_system.h"
#include "grid.h"

namespace Ariadne {

std::vector<std::string> DiscreteEvent::_names=std::vector<std::string>();

HybridSystem::~HybridSystem()
{
}

HybridSystem::HybridSystem()
{
}


std::ostream& operator<<(std::ostream& os, const HybridSystem& sys) {
    os << "HybridSystem(\n  "
       << "algebraic_equations=\n" << sys._algebraic_equations << ",\n"
       << "differential_equations=\n" << sys._differential_equations << ",\n"
       << "discrete_assignments=\n" << sys._discrete_assignments << ",\n"
       << "reset_equations=\n" << sys._update_equations << ",\n"
       << "guard_predicates=\n" << sys._guard_predicates << ",\n"
       << "invariant_predicates=\n" << sys._invariant_predicates << "\n)\n";
    return os;
}

std::ostream& operator<<(std::ostream& os, const HybridSystem::DifferentialEquation& de) {
    os << "\n" << de.loc << " -> dot("<<de.lhs<<")="<<de.rhs;
    return os;
}

std::ostream& operator<<(std::ostream& os, const HybridSystem::AlgebraicEquation& ae) {
    os << "\n" << ae.loc << " -> "<<ae.lhs<<"="<<ae.rhs;
    return os;
}

std::ostream& operator<<(std::ostream& os, const HybridSystem::DiscreteAssignment& da) {
    os << "\n" << da.e << "; " << da.loc << " -> "<<da.lhs.name()<<"="<<da.rhs;
    return os;
}

std::ostream& operator<<(std::ostream& os, const HybridSystem::UpdateEquation& da) {
    os << "\n" << da.e << "; " << da.loc << " -> "<<da.lhs.name()<<"'="<<da.rhs;
    return os;
}

std::ostream& operator<<(std::ostream& os, const HybridSystem::GuardPredicate& g) {
    os << "\n" << g.e << "; " << g.loc << ", "<<g.pred;
    return os;
}

std::ostream& operator<<(std::ostream& os, const HybridSystem::InvariantPredicate& inv) {
    os << "\n" << inv.loc << ", "<<inv.pred;
    return os;
}





}
