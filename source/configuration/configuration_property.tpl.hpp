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
#include "utility/randomiser.hpp"
#include "symbolic/identifier.hpp"
#include "configuration.hpp"
#include "configuration_property.hpp"
#include "configuration_property_path.hpp"
#include "searchable_configuration.hpp"
#include "configurable.hpp"

namespace Ariadne {

using std::min, std::max;

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

template<class T> RangeConfigurationProperty<T>::RangeConfigurationProperty(SearchSpaceConverterInterface<T> const& converter, ConfigurationPropertyRefinerInterface<T> const& refiner) :
        ConfigurationPropertyBase<T>(false), _lower(T()), _upper(T()), _refined(T()), _is_refined(false),
        _converter(SharedPointer<SearchSpaceConverterInterface<T>>(converter.clone())),
        _refiner(SharedPointer<ConfigurationPropertyRefinerInterface<T>>(refiner.clone())) { }

template<class T> RangeConfigurationProperty<T>::RangeConfigurationProperty(T const& lower, T const& upper,
        SearchSpaceConverterInterface<T> const& converter, ConfigurationPropertyRefinerInterface<T> const& refiner) :
        ConfigurationPropertyBase<T>(true), _lower(lower), _upper(upper), _refined(T()), _is_refined(false),
        _converter(SharedPointer<SearchSpaceConverterInterface<T>>(converter.clone())),
        _refiner(SharedPointer<ConfigurationPropertyRefinerInterface<T>>(refiner.clone())) {
    ARIADNE_PRECONDITION(not possibly(upper < lower));
}

template<class T> RangeConfigurationProperty<T>::RangeConfigurationProperty(T const& value,
        SearchSpaceConverterInterface<T> const& converter, ConfigurationPropertyRefinerInterface<T> const& refiner) :
                RangeConfigurationProperty(value,value,converter,refiner) { }

template<class T> T const& RangeConfigurationProperty<T>::get() const {
    ARIADNE_PRECONDITION(this->is_specified());
    ARIADNE_ASSERT_MSG(this->is_single() or is_refined(),"The property should have a single value or be under refinement when actually used. Are you accessing it outside the related task?");
    return (is_refined() ? _refined : _upper);
}

template<class T> Bool RangeConfigurationProperty<T>::is_single() const {
    if (not this->is_specified()) return false;
    else return possibly(_lower == _upper);
}

template<class T> Bool RangeConfigurationProperty<T>::is_refined() const {
    return _is_refined;
}

template<class T> Bool RangeConfigurationProperty<T>::is_metric(ConfigurationPropertyPath const& path) const {
    ARIADNE_PRECONDITION(path.is_root());
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
        if (not this->is_refined()) {
            int min_value = _converter->to_int(_lower);
            int max_value = _converter->to_int(_upper);
            ARIADNE_ASSERT_MSG(
                    not(max_value == std::numeric_limits<int>::max() and min_value < std::numeric_limits<int>::max()),
                    "An upper bounded range is required.");
            ARIADNE_ASSERT_MSG(
                    not(min_value == std::numeric_limits<int>::min() and max_value > std::numeric_limits<int>::min()),
                    "A lower bounded range is required.");
            if (min_value == max_value) result.push_back(min_value); // Necessary to address the +inf case
            else for (int i = min_value; i <= max_value; ++i) result.push_back(i);
        } else result.push_back(0);
    }
    return result;
}

template<class T> void RangeConfigurationProperty<T>::set_single(ConfigurationPropertyPath const& path, int integer_value) {
    ARIADNE_PRECONDITION(path.is_root());
    local_set_single(integer_value);
}

template<class T> void RangeConfigurationProperty<T>::refine_init(ConfigurationPropertyPath const& path) {
    ARIADNE_PRECONDITION(path.is_root());
    ARIADNE_ASSERT_MSG(not is_single(),"The property value must be in an interval in order to be refined.");
    ARIADNE_ASSERT_MSG(not is_refined(),"The property has already been initialised for refinement.");
    _refined = Randomiser<T>::get(_lower,_upper);
    _is_refined = true;
}

template<class T> void RangeConfigurationProperty<T>::refine_value(ConfigurationPropertyPath const& path, double ratio) {
    ARIADNE_PRECONDITION(path.is_root());
    ARIADNE_PRECONDITION(is_refined());
    ARIADNE_ASSERT_MSG(not is_single(),"The property value must be in an interval in order to be refined.");
    _refined = max(_lower,min(_upper,_refiner->apply(ratio,_refined)));
}

template<class T> void RangeConfigurationProperty<T>::local_set_single(int integer_value) {
    ARIADNE_PRECONDITION(not is_refined());
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

template<class T> ConfigurationPropertyInterface* RangeConfigurationProperty<T>::at(ConfigurationPropertyPath const& path) {
    ARIADNE_ASSERT_MSG(path.is_root(),"The path " << path << " is not a root but a range property can't have configurable objects below.");
    return this;
}

template<class T> void RangeConfigurationProperty<T>::set(T const& lower, T const& upper) {
    ARIADNE_PRECONDITION(not possibly(upper < lower));
    ARIADNE_ASSERT_MSG(not is_refined(),"The value is under refinement, setting it is not allowed.");
    this->set_specified();
    _lower = lower;
    _upper = upper;
}

template<class T> void RangeConfigurationProperty<T>::set(T const& value) {
    ARIADNE_ASSERT_MSG(not is_refined(),"The value is under refinement, setting it is not allowed.");
    this->set_specified();
    _lower = value;
    _upper = value;
}

template<class T> List<SharedPointer<T>> RangeConfigurationProperty<T>::values() const {
    List<SharedPointer<T>> result;
    if (this->is_specified()) {
        if (this->is_refined()) {
            result.append(std::make_shared<T>(_refined));
        } else {
            result.append(std::make_shared<T>(_lower));
            if (not is_single()) result.append(std::make_shared<T>(_upper));
        }
    }
    return result;
}

template<class T> EnumConfigurationProperty<T>::EnumConfigurationProperty()
        : ConfigurationPropertyBase<T>(false) {
    ARIADNE_PRECONDITION(std::is_enum<T>::value);
}

template<class T> EnumConfigurationProperty<T>::EnumConfigurationProperty(Set<T> const& values)
        : ConfigurationPropertyBase<T>(true), _values(values) {
    ARIADNE_PRECONDITION(std::is_enum<T>::value);
    ARIADNE_PRECONDITION(values.size()>0);
}

template<class T> EnumConfigurationProperty<T>::EnumConfigurationProperty(T const& value)
        : ConfigurationPropertyBase<T>(true) {
    ARIADNE_PRECONDITION(std::is_enum<T>::value);
    _values.insert(value);
}

template<class T> Bool EnumConfigurationProperty<T>::is_single() const {
    return (_values.size() == 1);
}

template<class T> Bool EnumConfigurationProperty<T>::is_metric(ConfigurationPropertyPath const& path) const {
    ARIADNE_PRECONDITION(path.is_root());
    return false;
}

template<class T> Bool EnumConfigurationProperty<T>::is_configurable() const {
    return false;
}

template<class T> SizeType EnumConfigurationProperty<T>::cardinality() const {
    return _values.size();
}

template<class T> List<int> EnumConfigurationProperty<T>::local_integer_values() const {
    List<int> result;
    for (SizeType i=0; i<_values.size(); ++i) result.push_back(i);
    return result;
}


template<class T> void EnumConfigurationProperty<T>::set_single(ConfigurationPropertyPath const& path, int integer_value) {
    ARIADNE_PRECONDITION(path.is_root());
    local_set_single(integer_value);
}

template<class T> void EnumConfigurationProperty<T>::local_set_single(int integer_value) {
    ARIADNE_PRECONDITION(not is_single());
    ARIADNE_PRECONDITION(integer_value >= 0 and integer_value < (int)cardinality());
    auto iter = _values.begin();
    for (SizeType i=0;i<(SizeType)integer_value;++i) ++iter;
    T value = *iter;
    _values.clear();
    _values.insert(value);
}

template<class T> void EnumConfigurationProperty<T>::refine_init(ConfigurationPropertyPath const& path) {
    ARIADNE_PRECONDITION(path.is_root());
    ARIADNE_ERROR("The value of an enum property cannot be refined.");
}

template<class T> void EnumConfigurationProperty<T>::refine_value(ConfigurationPropertyPath const& path, double ratio) {
    ARIADNE_PRECONDITION(path.is_root());
    ARIADNE_ERROR("The value of an enum property cannot be refined.");
}

template<class T> ConfigurationPropertyInterface* EnumConfigurationProperty<T>::clone() const {
    return new EnumConfigurationProperty(*this);
}

template<class T> ConfigurationPropertyInterface* EnumConfigurationProperty<T>::at(ConfigurationPropertyPath const& path) {
    ARIADNE_ASSERT_MSG(path.is_root(),"The path " << path << " is not a root but an enum property can't have configurable objects below.");
    return this;
}

template<class T> T const& EnumConfigurationProperty<T>::get() const {
    ARIADNE_PRECONDITION(this->is_specified());
    ARIADNE_ASSERT_MSG(this->is_single(),"The property should have a single value when actually used. Are you accessing it outside the related task?");
    return *_values.begin();
}

template<class T> void EnumConfigurationProperty<T>::set(T const& value) {
    this->set_specified();
    _values.clear();
    _values.insert(value);
}

template<class T> void EnumConfigurationProperty<T>::set(Set<T> const& values) {
    ARIADNE_PRECONDITION(not values.empty());
    this->set_specified();
    _values = values;
}

template<class T> List<SharedPointer<T>> EnumConfigurationProperty<T>::values() const {
    List<SharedPointer<T>> result;
    for (auto v : _values) result.push_back(SharedPointer<T>(new T(v)));
    return result;
}

template<class T> ListConfigurationProperty<T>::ListConfigurationProperty()
    : ConfigurationPropertyBase<T>(false)
{ }

template<class T> ListConfigurationProperty<T>::ListConfigurationProperty(List<T> const& values)
    : ConfigurationPropertyBase<T>(true), _values(values) {
        ARIADNE_PRECONDITION(values.size()>0);
}

template<class T> ListConfigurationProperty<T>::ListConfigurationProperty(T const& value)
    : ConfigurationPropertyBase<T>(true) {
        _values.push_back(value);
}

template<class T> Bool ListConfigurationProperty<T>::is_single() const {
    return (_values.size() == 1);
}

template<class T> Bool ListConfigurationProperty<T>::is_metric(ConfigurationPropertyPath const& path) const {
    ARIADNE_PRECONDITION(path.is_root());
    return false;
}

template<class T> Bool ListConfigurationProperty<T>::is_configurable() const {
    return false;
}

template<class T> SizeType ListConfigurationProperty<T>::cardinality() const {
    return _values.size();
}

template<class T> List<int> ListConfigurationProperty<T>::local_integer_values() const {
    List<int> result;
    for (SizeType i=0; i<_values.size(); ++i) result.push_back(i);
    return result;
}

template<class T> void ListConfigurationProperty<T>::set_single(ConfigurationPropertyPath const& path, int integer_value) {
    ARIADNE_PRECONDITION(path.is_root());
    local_set_single(integer_value);
}

template<class T> void ListConfigurationProperty<T>::local_set_single(int integer_value) {
    ARIADNE_PRECONDITION(not is_single());
    ARIADNE_PRECONDITION(integer_value >= 0 and integer_value < (int)cardinality());
    T value = _values[(SizeType)integer_value];
    _values.clear();
    _values.push_back(value);
}

template<class T> void ListConfigurationProperty<T>::refine_init(ConfigurationPropertyPath const& path) {
    ARIADNE_PRECONDITION(path.is_root());
    ARIADNE_ERROR("The value of a list property cannot be refined.");
}

template<class T> void ListConfigurationProperty<T>::refine_value(ConfigurationPropertyPath const& path, double ratio) {
    ARIADNE_PRECONDITION(path.is_root());
    ARIADNE_ERROR("The value of a list property cannot be refined.");
}

template<class T> ConfigurationPropertyInterface* ListConfigurationProperty<T>::clone() const {
    return new ListConfigurationProperty(*this);
}

template<class T> ConfigurationPropertyInterface* ListConfigurationProperty<T>::at(ConfigurationPropertyPath const& path) {
    ARIADNE_ASSERT_MSG(path.is_root(),"The path " << path << " is not a root but a list property can't have configurable objects below.");
    return this;
}

template<class T> T const& ListConfigurationProperty<T>::get() const {
    ARIADNE_PRECONDITION(this->is_specified());
    ARIADNE_ASSERT_MSG(this->is_single(),"The property should have a single value when actually used. Are you accessing it outside the related task?");
    return _values.back();
}

template<class T> void ListConfigurationProperty<T>::set(T const& value) {
    this->set_specified();
    _values.clear();
    _values.push_back(value);
}

template<class T> void ListConfigurationProperty<T>::set(List<T> const& values) {
    ARIADNE_PRECONDITION(not values.empty());
    this->set_specified();
    _values = values;
}

template<class T> List<SharedPointer<T>> ListConfigurationProperty<T>::values() const {
    List<SharedPointer<T>> result;
    for (auto v : _values) result.push_back(SharedPointer<T>(new T(v)));
    return result;
}

template<class T> InterfaceConfigurationProperty<T>::InterfaceConfigurationProperty() : ConfigurationPropertyBase<T>(false) { }

template<class T> InterfaceConfigurationProperty<T>::InterfaceConfigurationProperty(List<SharedPointer<T>> const& list) : ConfigurationPropertyBase<T>(true), _values(list) {
    ARIADNE_PRECONDITION(list.size()>0);
}

template<class T> InterfaceConfigurationProperty<T>::InterfaceConfigurationProperty(T const& value) : ConfigurationPropertyBase<T>(true) {
    _values.push_back(SharedPointer<T>(value.clone()));
}

template<class T> Bool InterfaceConfigurationProperty<T>::is_single() const {
    return (_values.size() == 1);
}

template<class T> Bool InterfaceConfigurationProperty<T>::is_metric(ConfigurationPropertyPath const& path) const {
    if (path.is_root()) return false;
    if (is_configurable()) {
        ARIADNE_PRECONDITION(is_single());
        auto properties = dynamic_cast<ConfigurableInterface*>(_values.back().get())->searchable_configuration().properties();
        auto p_ptr = properties.find(path.first());
        if (p_ptr != properties.end()) {
            return p_ptr->second->is_metric(path.subpath());
        } else {
            ARIADNE_FAIL_MSG("A property for " << path << " has not been found.");
        }
    } else {
        ARIADNE_FAIL_MSG("The object is not configurable, a property for " << path << " could not been found.");
    }
}

template<class T> Bool InterfaceConfigurationProperty<T>::is_configurable() const {
    ARIADNE_ASSERT_MSG(this->is_specified(),"Cannot check if configurable if the property is not specified.");
    ARIADNE_PRECONDITION(is_single());
    auto configurable_interface_ptr = dynamic_cast<ConfigurableInterface*>(_values.back().get());
    return (configurable_interface_ptr != nullptr);
}

template<class T> SizeType InterfaceConfigurationProperty<T>::cardinality() const {
    return _values.size();
}

template<class T> List<int> InterfaceConfigurationProperty<T>::local_integer_values() const {
    List<int> result;
    for (SizeType i=0; i<_values.size(); ++i) result.push_back(i);
    return result;
}

template<class T> void InterfaceConfigurationProperty<T>::local_set_single(int integer_value) {
    ARIADNE_PRECONDITION(not is_single());
    ARIADNE_PRECONDITION(integer_value >= 0 and integer_value < (int)cardinality());
    SharedPointer<T> value = _values[(SizeType)integer_value];
    _values.clear();
    _values.push_back(value);
}

template<class T> void InterfaceConfigurationProperty<T>::set_single(ConfigurationPropertyPath const& path, int integer_value) {
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

template<class T> void InterfaceConfigurationProperty<T>::refine_init(ConfigurationPropertyPath const& path) {
    ARIADNE_ASSERT_MSG(not path.is_root(),"The value of an interface configuration property can not be refined.");
    ARIADNE_ASSERT_MSG(is_configurable(),"The path is incorrect since the object in the interface configuration property is not configurable.");
    ARIADNE_ASSERT_MSG(is_single(),"Can't reach the path since there are multiple values and multiple configurable objects are not allowed");
    auto properties = dynamic_cast<ConfigurableInterface*>(_values.back().get())->searchable_configuration().properties();
    auto p_ptr = properties.find(path.first());
    if (p_ptr != properties.end()) p_ptr->second->refine_init(path.subpath());
    else { ARIADNE_ERROR("The path " << path << " has not been found."); }
}

template<class T> void InterfaceConfigurationProperty<T>::refine_value(ConfigurationPropertyPath const& path, double ratio) {
    ARIADNE_ASSERT_MSG(not path.is_root(),"The value of an interface configuration property can not be refined.");
    ARIADNE_ASSERT_MSG(is_configurable(),"The path is incorrect since the object in the interface configuration property is not configurable.");
    ARIADNE_ASSERT_MSG(is_single(),"Can't reach the path since there are multiple values and multiple configurable objects are not allowed");
    auto properties = dynamic_cast<ConfigurableInterface*>(_values.back().get())->searchable_configuration().properties();
    auto p_ptr = properties.find(path.first());
    if (p_ptr != properties.end()) p_ptr->second->refine_value(path.subpath(),ratio);
    else { ARIADNE_ERROR("The path " << path << " has not been found."); }
}

template<class T> Map<ConfigurationPropertyPath,List<int>> InterfaceConfigurationProperty<T>::integer_values() const {
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

template<class T> ConfigurationPropertyInterface* InterfaceConfigurationProperty<T>::clone() const {
    List<SharedPointer<T>> values;
    for (auto ptr : _values) values.push_back(SharedPointer<T>(ptr->clone()));
    return new InterfaceConfigurationProperty(values);
}

template<class T> ConfigurationPropertyInterface* InterfaceConfigurationProperty<T>::at(ConfigurationPropertyPath const& path) {
    if (path.is_root()) return this;
    else {
        ARIADNE_ASSERT_MSG(is_configurable(),"The object held is not configurable, path error.");
        ARIADNE_ASSERT_MSG(is_single(),"Cannot retrieve properties if the list has multiple objects.");
        auto configurable_ptr = dynamic_cast<ConfigurableInterface*>(_values.back().get());
        auto properties = configurable_ptr->searchable_configuration().properties();
        auto prop_ptr = properties.find(path.first());
        ARIADNE_ASSERT_MSG(prop_ptr != properties.end(),"The property '" << path.first() << "' was not found in the configuration.");
        return prop_ptr->second->at(path.subpath());
    }
}

template<class T> T const& InterfaceConfigurationProperty<T>::get() const {
    ARIADNE_PRECONDITION(this->is_specified());
    ARIADNE_ASSERT_MSG(this->is_single(),"The property should have a single value when actually used. Are you accessing it outside the related task?");
    return *_values.back();
}

template<class T> void InterfaceConfigurationProperty<T>::set(T const& value) {
    this->set_specified();
    _values.clear();
    _values.push_back(SharedPointer<T>(value.clone()));
}

template<class T> void InterfaceConfigurationProperty<T>::set(SharedPointer<T> const& value) {
    ARIADNE_PRECONDITION(value != nullptr);
    this->set_specified();
    _values.clear();
    _values.push_back(value);
}

template<class T> void InterfaceConfigurationProperty<T>::set(List<SharedPointer<T>> const& values) {
    ARIADNE_PRECONDITION(values.size()>0);
    this->set_specified();
    _values = values;
}

template<class T> List<SharedPointer<T>> InterfaceConfigurationProperty<T>::values() const {
    return _values;
}

} // namespace Ariadne


#endif // ARIADNE_CONFIGURATION_PROPERTY_TPL_HPP