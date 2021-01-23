/***************************************************************************
 *            concurrency/searchable_configuration.hpp
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

/*! \file concurrency/searchable_configuration.hpp
 *  \brief Classes for configuration and related helpers.
 */

#ifndef ARIADNE_SEARCHABLE_CONFIGURATION_HPP
#define ARIADNE_SEARCHABLE_CONFIGURATION_HPP

#include <ostream>
#include <type_traits>
#include "utility/writable.hpp"
#include "utility/macros.hpp"
#include "utility/container.hpp"
#include "utility/pointer.hpp"
#include "symbolic/identifier.hpp"
#include "configuration.hpp"
#include "configuration_property_interface.hpp"

namespace Ariadne {

class TaskSearchSpace;

//! \brief Extension of ConfigurationInterface to deal with search in the properties space
class SearchableConfiguration : public ConfigurationInterface {
  public:
    SearchableConfiguration() = default;
    SearchableConfiguration(SearchableConfiguration const& c);
    SearchableConfiguration& operator=(SearchableConfiguration const& c);
    virtual ~SearchableConfiguration() = default;

    //! \brief Construct a search space from the current configuration
    TaskSearchSpace search_space() const;

    //! \brief If the configuration is made of single values
    Bool is_singleton() const;

    Map<Identifier,SharedPointer<ConfigurationPropertyInterface>>& properties();
    Map<Identifier,SharedPointer<ConfigurationPropertyInterface>> const& properties() const;

    //! \brief Accessors for get and set of a property identified by \a identifier of type \a P
    template<class P> P const& at(Identifier const& identifier) const { return static_cast<P const&>(*_properties.get(identifier)); }
    template<class P> P& at(Identifier const& identifier) { return static_cast<P&>(*_properties.get(identifier)); }

    //! \brief Add a property to the configuration
    void add_property(Identifier const& name, ConfigurationPropertyInterface const& property);
    //! \brief Add a property to the configuration from another configuration
    //! \details This is used to copy a property from a user configuration onto the base class configuration
    void add_shared_property(Identifier const& name, SearchableConfiguration const& c);

    OutputStream& _write(OutputStream& os) const override;
  private:
    Map<Identifier,SharedPointer<ConfigurationPropertyInterface>> _properties;
};

} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_HPP
