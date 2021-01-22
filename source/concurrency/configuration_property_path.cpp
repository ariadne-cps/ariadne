/***************************************************************************
 *            concurrency/configuration_property_path.cpp
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

#include "../utility/macros.hpp"
#include "../symbolic/identifier.hpp"
#include "../concurrency/configuration_property_path.hpp"

namespace Ariadne {

ConfigurationPropertyPath::ConfigurationPropertyPath(ConfigurationPropertyPath const& path) {
    _path = path._path;
}

ConfigurationPropertyPath& ConfigurationPropertyPath::operator=(ConfigurationPropertyPath const& path) {
    _path = path._path;
    return *this;
}

Bool ConfigurationPropertyPath::operator<(ConfigurationPropertyPath const& path) const {
    if (_path.size() < path._path.size()) return true;
    else if (_path.size() > path._path.size()) return false;
    else return (repr() < path.repr());
}

Identifier ConfigurationPropertyPath::repr() const {
    std::ostringstream sstream;
    sstream << *this;
    return sstream.str();
}

Bool ConfigurationPropertyPath::is_root() const {
    return _path.empty();
}

void ConfigurationPropertyPath::append(Identifier const& node) {
    ARIADNE_PRECONDITION(not node.empty());
    _path.push_back(node);
}

void ConfigurationPropertyPath::prepend(Identifier const& node) {
    ARIADNE_PRECONDITION(not node.empty());
    _path.push_front(node);
}

ConfigurationPropertyPath ConfigurationPropertyPath::subpath() const {
    auto result = *this;
    result._path.pop_front();
    return result;
}

OutputStream& ConfigurationPropertyPath::_write(OutputStream& os) const {
    auto size = _path.size();
    auto iter = _path.begin();
    os << "./";
    for (SizeType i=0; i<size; ++i) {
        os << *iter << "/";
        ++iter;
    }
    return os;
}

} // namespace Ariadne
