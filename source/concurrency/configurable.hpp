/***************************************************************************
 *            concurrency/configurable.hpp
 *
 *  Copyright  2011-20  Luca Geretti
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

/*! \file concurrency/configurable.hpp
 *  \brief Classes for configurable tasks.
 */

#ifndef ARIADNE_CONFIGURABLE_HPP
#define ARIADNE_CONFIGURABLE_HPP

#include <ostream>
#include <type_traits>
#include "utility/writable.hpp"
#include "utility/macros.hpp"
#include "utility/container.hpp"
#include "utility/pointer.hpp"
#include "symbolic/identifier.hpp"
#include "numeric/builtin.hpp"
#include "numeric/real.hpp"
#include "solvers/configuration.hpp"

namespace Ariadne {

class SearchableConfiguration;

class ConfigurableInterface {
public:
    virtual SearchableConfiguration const& searchable_configuration() const = 0;
};

//! \brief Base template class to be specialised while deriving from SearchableConfigurationInterface
template<class C> class Configuration;// : public SearchableConfiguration { };

//! \brief Is-a component that provides a configuration
//! \details Since the configuration returned is const, a Configurable object should be constructed from
//! a pre-set configuration. If the configuration must specify certain properties or if some properties
//! must be coherent with the Configurable (e.g., the system used by a Configurable evolver), then a Builder
//! approach should be used for creation of the configuration, which should become an immutable object.
template<class C>
class Configurable : public ConfigurableInterface {
    friend class Configuration<C>;
public:
    Configurable(Configuration<C> const& config);
    Configuration<C> const& configuration() const;
    SearchableConfiguration const& searchable_configuration() const override;
private:
    SharedPointer<Configuration<C>> _configuration;
};

} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_HPP
