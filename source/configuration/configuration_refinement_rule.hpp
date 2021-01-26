/***************************************************************************
 *            configuration/configuration_refinement_rule.hpp
 *
 *  Copyright  2007-20  Luca Geretti
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

/*! \file configuration/configuration_refinement_rule.hpp
 *  \brief A rule for refining a configuration based on a task result.
 */

#ifndef ARIADNE_CONFIGURATION_REFINEMENT_RULE_HPP
#define ARIADNE_CONFIGURATION_REFINEMENT_RULE_HPP

#include <functional>
#include "configurable.hpp"
#include "utility/pointer.hpp"

namespace Ariadne {

template<class R> struct TaskInput;
template<class R> struct TaskOutput;

template<class R>
class ConfigurationRefinementRule {
  public:
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;
    typedef Configuration<R> ConfigurationType;
    typedef std::function<void(InputType const&, OutputType const&, ConfigurationType&)> FunctionType;

    ConfigurationRefinementRule(FunctionType const& function) : _function(function) { }
    void operator()(InputType const& i, OutputType const& o, ConfigurationType& c) { _function(i,o,c); }
    Bool operator<(ConfigurationRefinementRule<R> const& r) const { return &_function < &r._function; }
  private:
    FunctionType const _function;
};


} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_REFINEMENT_RULE_HPP
