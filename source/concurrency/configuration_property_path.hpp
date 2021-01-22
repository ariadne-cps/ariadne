/***************************************************************************
 *            concurrency/configuration_property_path.hpp
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

/*! \file concurrency/configuration_property_path.hpp
 *  \brief A path across a hierarchic configuration object.
 */

#ifndef ARIADNE_CONFIGURATION_PROPERTY_PATH_HPP
#define ARIADNE_CONFIGURATION_PROPERTY_PATH_HPP

#include <utility>
#include <deque>
#include "../utility/typedefs.hpp"
#include "../utility/writable.hpp"

namespace Ariadne {

class Identifier;

class ConfigurationPropertyPath : public WritableInterface {
  public:
    ConfigurationPropertyPath() = default;
    ConfigurationPropertyPath(ConfigurationPropertyPath const& path);
    ConfigurationPropertyPath& operator=(ConfigurationPropertyPath const& path);
    Bool operator<(ConfigurationPropertyPath const& path) const;

    Identifier repr() const;

    Bool is_root() const;
    void append(Identifier const& node);
    void prepend(Identifier const& node);
    ConfigurationPropertyPath subpath() const;

    OutputStream& _write(OutputStream& os) const override;
  private:
    std::deque<Identifier> _path;
};

} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_PROPERTY_PATH_HPP
