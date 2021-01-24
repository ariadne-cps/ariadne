/***************************************************************************
 *            concurrency/configurable.tpl.hpp
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

#ifndef ARIADNE_CONFIGURABLE_TPL_HPP
#define ARIADNE_CONFIGURABLE_TPL_HPP

#include "utility/macros.hpp"
#include "concurrency/task_search_point.hpp"
#include "concurrency/task_search_space.hpp"
#include "configurable.hpp"

namespace Ariadne {

template<class C> Configurable<C>::Configurable(Configuration<C> const& config) : _configuration(new Configuration<C>(config)) { }

template<class C> Configuration<C> const& Configurable<C>::configuration() const {
    return *_configuration;
}

template<class C> SearchableConfiguration const& Configurable<C>::searchable_configuration() const {
    return dynamic_cast<SearchableConfiguration const &>(*_configuration);
}

template<class C> ConfigurableBase<C>::ConfigurableBase(Configuration<C> const& config) : _configuration(new Configuration<C>(config)) { }

template<class C> Configuration<C> const& ConfigurableBase<C>::configuration() const {
    return *_configuration;
}

//! \brief Make a configuration from another configuration \a cfg and a point \a p in the search space
template<class C> Configuration<C> make_singleton(Configuration<C> const& cfg, TaskSearchPoint const& p) {
    ARIADNE_PRECONDITION(not cfg.is_singleton());
    auto result = cfg;
    for (auto param : p.space().parameters()) {
        auto prop_ptr = result.properties().find(param.path().first());
        ARIADNE_ASSERT_MSG(prop_ptr != cfg.properties().end(), "The TaskSearchPoint parameter '" << param.path() << "' is not in the configuration.");
        prop_ptr->second->set_single(param.path().subpath(),p.value(param.path()));
    }
    ARIADNE_ASSERT_MSG(result.is_singleton(),"There are missing parameters in the search point, since the configuration could not be made singleton.");
    return result;
}

} // namespace Ariadne

#endif // ARIADNE_CONFIGURABLE_TPL_HPP
