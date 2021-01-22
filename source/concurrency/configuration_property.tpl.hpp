/***************************************************************************
 *            concurrency/configuration_property.tpl.hpp
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

#ifndef ARIADNE_CONFIGURATION_PROPERTY_TPL_HPP
#define ARIADNE_CONFIGURATION_PROPERTY_TPL_HPP

#include <ostream>
#include <type_traits>
#include "utility/writable.hpp"
#include "utility/macros.hpp"
#include "utility/container.hpp"
#include "utility/pointer.hpp"
#include "symbolic/identifier.hpp"
#include "solvers/configuration.hpp"
#include "concurrency/configuration_property.hpp"
#include "concurrency/configuration_property_path.hpp"
#include "concurrency/searchable_configuration.hpp"
#include "concurrency/configurable.hpp"

namespace Ariadne {

template<class T> ConfigurationPropertyBase<T>::ConfigurationPropertyBase(Bool const& is_specified) : _is_specified(is_specified) { }

template<class T> void ConfigurationPropertyBase<T>::set_specified() {
    _is_specified = true;
}

template<class T> Bool ConfigurationPropertyBase<T>::is_specified() const {
    return _is_specified;
}

template<class T> Map<ConfigurationPropertyPath,List<int>> ConfigurationPropertyBase<T>::integer_values() const {
    Map<ConfigurationPropertyPath,List<int>> result;
    result.insert(Pair<ConfigurationPropertyPath,List<int>>(ConfigurationPropertyPath(), local_integer_values()));
    return result;
}

template<class T> OutputStream& ConfigurationPropertyBase<T>::_write(OutputStream& os) const {
    auto vals = values();
    if (vals.empty()) { os << "<unspecified>"; }
    else if (vals.size() == 1) { os << *vals[0]; }
    else {
        os << "{";
        for (SizeType i=0; i<vals.size()-1; ++i) os << *vals[i] << ",";
        os << *vals[vals.size()-1] << "}";
    }
    return os;
}

template<class T> RangeConfigurationProperty<T>::RangeConfigurationProperty(SearchSpaceConverterInterface<T> const& converter) :
        ConfigurationPropertyBase<T>(false), _lower(T()), _upper(T()),
        _converter(SharedPointer<SearchSpaceConverterInterface<T>>(converter.clone())) { }

template<class T> RangeConfigurationProperty<T>::RangeConfigurationProperty(T const& lower, T const& upper, SearchSpaceConverterInterface<T> const& converter) :
        ConfigurationPropertyBase<T>(true), _lower(lower), _upper(upper),
        _converter(SharedPointer<SearchSpaceConverterInterface<T>>(converter.clone())) {
    ARIADNE_PRECONDITION(not possibly(upper < lower));
}

template<class T> RangeConfigurationProperty<T>::RangeConfigurationProperty(T const& value, SearchSpaceConverterInterface<T> const& converter) :
        ConfigurationPropertyBase<T>(true), _lower(value), _upper(value),
        _converter(SharedPointer<SearchSpaceConverterInterface<T>>(converter.clone())) { }

template<class T> T const& RangeConfigurationProperty<T>::get() const {
    ARIADNE_PRECONDITION(this->is_specified());
    ARIADNE_ASSERT_MSG(this->is_single(),"The property should have a single value when actually used. Are you accessing it outside the related task?");
    return _upper;
}

template<class T> Bool RangeConfigurationProperty<T>::is_single() const {
    if (not this->is_specified()) return false;
    else return possibly(_lower == _upper);
}

template<class T> Bool RangeConfigurationProperty<T>::is_metric() const {
    return true;
}

template<class T> Bool RangeConfigurationProperty<T>::is_configurable() const {
    return false;
}

template<class T> SizeType RangeConfigurationProperty<T>::cardinality() const {
    if (is_single()) return 1;
    else if (not this->is_specified()) return 0;
    else return 1+(SizeType)(_converter->to_int(_upper) - _converter->to_int(_lower));
}

template<class T> List<int> RangeConfigurationProperty<T>::local_integer_values() const {
    List<int> result;
    if (this->is_specified()) {
        int min_value = _converter->to_int(_lower);
        int max_value = _converter->to_int(_upper);
        for (int i=min_value; i<=max_value; ++i) result.push_back(i);
    }
    return result;
}

template<class T> void RangeConfigurationProperty<T>::set_single(ConfigurationPropertyPath const& path, int integer_value) {
    ARIADNE_PRECONDITION(path.is_root());
    local_set_single(integer_value);
}

template<class T> void RangeConfigurationProperty<T>::local_set_single(int integer_value) {
    int min_value = _converter->to_int(_lower);
    int max_value = _converter->to_int(_upper);
    ARIADNE_PRECONDITION(not is_single());
    ARIADNE_PRECONDITION(integer_value >= min_value and integer_value <= max_value);
    if (integer_value == min_value) _upper = _lower; // Avoids rounding error
    else if (integer_value == max_value) _lower = _upper; // Avoids rounding error
    else { _lower = _upper = _converter->to_value(integer_value); }
}

template<class T> ConfigurationPropertyInterface* RangeConfigurationProperty<T>::clone() const {
    return new RangeConfigurationProperty(*this);
}

template<class T> void RangeConfigurationProperty<T>::set(T const& lower, T const& upper) {
    ARIADNE_PRECONDITION(not possibly(upper < lower));
    this->set_specified();
    _lower = lower;
    _upper = upper;
}

template<class T> void RangeConfigurationProperty<T>::set(T const& value) {
    this->set_specified();
    _lower = value;
    _upper = value;
}

template<class T> List<SharedPointer<T>> RangeConfigurationProperty<T>::values() const {
    List<SharedPointer<T>> result;
    if (this->is_specified()) {
        result.append(SharedPointer<T>(new T(_lower)));
        if (not is_single()) result.append(SharedPointer<T>(new T(_upper)));
    }
    return result;
}

template<class T> SetConfigurationProperty<T>::SetConfigurationProperty()
    : ConfigurationPropertyBase<T>(false)
{ }

template<class T> SetConfigurationProperty<T>::SetConfigurationProperty(Set<T> const& values)
    : ConfigurationPropertyBase<T>(true), _values(values) {
        ARIADNE_PRECONDITION(values.size()>0);
}

template<class T> SetConfigurationProperty<T>::SetConfigurationProperty(T const& value)
    : ConfigurationPropertyBase<T>(true) {
        _values.insert(value);
}

template<class T> Bool SetConfigurationProperty<T>::is_single() const {
    return (_values.size() == 1);
}

template<class T> Bool SetConfigurationProperty<T>::is_metric() const {
    return false;
}

template<class T> Bool SetConfigurationProperty<T>::is_configurable() const {
    return false;
}

template<class T> SizeType SetConfigurationProperty<T>::cardinality() const {
    return _values.size();
}

template<class T> List<int> SetConfigurationProperty<T>::local_integer_values() const {
    List<int> result;
    for (SizeType i=0; i<_values.size(); ++i) result.push_back(i);
    return result;
}

template<class T> void SetConfigurationProperty<T>::set_single(ConfigurationPropertyPath const& path, int integer_value) {
    ARIADNE_PRECONDITION(path.is_root());
    local_set_single(integer_value);
}

template<class T> void SetConfigurationProperty<T>::local_set_single(int integer_value) {
    ARIADNE_PRECONDITION(not is_single());
    ARIADNE_PRECONDITION(integer_value >= 0 and integer_value < (int)cardinality());
    auto iter = _values.begin();
    for (SizeType i=0;i<(SizeType)integer_value;++i) ++iter;
    T value = *iter;
    _values.clear();
    _values.insert(value);
}

template<class T> ConfigurationPropertyInterface* SetConfigurationProperty<T>::clone() const {
    return new SetConfigurationProperty(*this);
}

template<class T> T const& SetConfigurationProperty<T>::get() const {
    ARIADNE_PRECONDITION(this->is_specified());
    ARIADNE_ASSERT_MSG(this->is_single(),"The property should have a single value when actually used. Are you accessing it outside the related task?");
    return *_values.begin();
}

template<class T> void SetConfigurationProperty<T>::set(T const& value) {
    this->set_specified();
    _values.clear();
    _values.insert(value);
}

template<class T> void SetConfigurationProperty<T>::set(Set<T> const& values) {
    ARIADNE_PRECONDITION(not values.empty());
    this->set_specified();
    _values = values;
}

template<class T> List<SharedPointer<T>> SetConfigurationProperty<T>::values() const {
    List<SharedPointer<T>> result;
    for (auto v : _values) result.push_back(SharedPointer<T>(new T(v)));
    return result;
}

template<class T> ListConfigurationProperty<T>::ListConfigurationProperty() : ConfigurationPropertyBase<T>(false) { }

template<class T> ListConfigurationProperty<T>::ListConfigurationProperty(List<SharedPointer<T>> const& list) : ConfigurationPropertyBase<T>(true), _values(list) {
    ARIADNE_PRECONDITION(list.size()>0);
}

template<class T> ListConfigurationProperty<T>::ListConfigurationProperty(T const& value) : ConfigurationPropertyBase<T>(true) {
    _values.push_back(SharedPointer<T>(value.clone()));
}

template<class T> Bool ListConfigurationProperty<T>::is_single() const {
    return (_values.size() == 1);
}

template<class T> Bool ListConfigurationProperty<T>::is_metric() const {
    return false;
}

template<class T> Bool ListConfigurationProperty<T>::is_configurable() const {
    ARIADNE_ASSERT_MSG(this->is_specified(),"Cannot check if configurable if the property is not specified.");
    ARIADNE_PRECONDITION(is_single());
    auto configurable_interface_ptr = dynamic_cast<ConfigurableInterface*>(_values.back().get());
    return (configurable_interface_ptr != nullptr);
}

template<class T> SizeType ListConfigurationProperty<T>::cardinality() const {
    return _values.size();
}

template<class T> List<int> ListConfigurationProperty<T>::local_integer_values() const {
    List<int> result;
    for (SizeType i=0; i<_values.size(); ++i) result.push_back(i);
    return result;
}

template<class T> void ListConfigurationProperty<T>::local_set_single(int integer_value) {
    ARIADNE_PRECONDITION(not is_single());
    ARIADNE_PRECONDITION(integer_value >= 0 and integer_value < (int)cardinality());
    SharedPointer<T> value = _values[(SizeType)integer_value];
    _values.clear();
    _values.push_back(value);
}

template<class T> void ListConfigurationProperty<T>::set_single(ConfigurationPropertyPath const& path, int integer_value) {
    if (path.is_root()) {
        local_set_single(integer_value);
    } else { // NOTE : we assume that we already checked for being single when getting the integer_values
        Bool been_set = false;
        auto configurable_interface_ptr = dynamic_cast<ConfigurableInterface*>(_values.back().get());
        if (configurable_interface_ptr != nullptr) {
            auto properties = configurable_interface_ptr->searchable_configuration().properties();
            auto p_ptr = properties.find(path.first());
            if (p_ptr != properties.end()) {
                p_ptr->second->set_single(path.subpath(),integer_value);
                been_set = true;
            }
        }
        ARIADNE_ASSERT_MSG(been_set,"A property for " << path << " has not been found.");
    }
}

template<class T> Map<ConfigurationPropertyPath,List<int>> ListConfigurationProperty<T>::integer_values() const {
    Map<ConfigurationPropertyPath,List<int>> result;
    result.insert(Pair<ConfigurationPropertyPath,List<int>>(ConfigurationPropertyPath(),local_integer_values()));
    if (is_single()) { // NOTE: we could extend to multiple values by using indexes
        auto configurable_interface_ptr = dynamic_cast<ConfigurableInterface*>(_values.back().get());
        if (configurable_interface_ptr != nullptr) {
            for (auto p : configurable_interface_ptr->searchable_configuration().properties()) {
                for (auto entry : p.second->integer_values()) {
                    auto prefixed_path = entry.first;
                    prefixed_path.prepend(p.first);
                    result.insert(Pair<ConfigurationPropertyPath,List<int>>(prefixed_path,entry.second));
                }
            }
        }
    }
    return result;
}

template<class T> ConfigurationPropertyInterface* ListConfigurationProperty<T>::clone() const {
    return new ListConfigurationProperty(*this);
}

template<class T> T const& ListConfigurationProperty<T>::get() const {
    ARIADNE_PRECONDITION(this->is_specified());
    ARIADNE_ASSERT_MSG(this->is_single(),"The property should have a single value when actually used. Are you accessing it outside the related task?");
    return *_values.back();
}

template<class T> void ListConfigurationProperty<T>::set(T const& value) {
    this->set_specified();
    _values.clear();
    _values.push_back(SharedPointer<T>(value.clone()));
}

template<class T> void ListConfigurationProperty<T>::set(SharedPointer<T> const& value) {
    ARIADNE_PRECONDITION(value != nullptr);
    this->set_specified();
    _values.clear();
    _values.push_back(value);
}

template<class T> void ListConfigurationProperty<T>::set(List<SharedPointer<T>> const& values) {
    ARIADNE_PRECONDITION(values.size()>0);
    this->set_specified();
    _values = values;
}

template<class T> List<SharedPointer<T>> ListConfigurationProperty<T>::values() const {
    return _values;
}

} // namespace Ariadne


#endif // ARIADNE_CONFIGURATION_PROPERTY_TPL_HPP