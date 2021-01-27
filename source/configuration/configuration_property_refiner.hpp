/***************************************************************************
 *            configuration/configuration_property_refiner.hpp
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

/*! \file configuration/configuration_property_refiner.hpp
 *  \brief Classes for refining a value of a configuration property
 */

#ifndef ARIADNE_CONFIGURATION_PROPERTY_REFINER_HPP
#define ARIADNE_CONFIGURATION_PROPERTY_REFINER_HPP

namespace Ariadne {

//! \brief Interface for conversion from/into the integer search space
template<class T>
struct ConfigurationPropertyRefinerInterface {
    virtual T get(T const& current, double coefficient) const = 0;

    virtual ConfigurationPropertyRefinerInterface* clone() const = 0;
    virtual ~ConfigurationPropertyRefinerInterface() = default;
};

template<class T> struct ProportionalRefiner;

template<> struct ProportionalRefiner<ExactDouble> : public ConfigurationPropertyRefinerInterface<ExactDouble> {
    ProportionalRefiner(double gain) : _gain(gain) { }
    ExactDouble get(ExactDouble const& current, double coefficient) const override {
        return ExactDouble(current.get_d()*(1.0 - _gain*coefficient));
    }
    ConfigurationPropertyRefinerInterface* clone() const override { return new ProportionalRefiner(*this); }
  private:
    double _gain;
};

template<> struct ProportionalRefiner<DegreeType> : public ConfigurationPropertyRefinerInterface<DegreeType> {
    ProportionalRefiner(double gain) : _gain(gain) { }
    DegreeType get(DegreeType const& current, double coefficient) const override {
        return DegreeType(current*(1.0 - _gain*coefficient));
    }
    ConfigurationPropertyRefinerInterface* clone() const override { return new ProportionalRefiner(*this); }
private:
    double _gain;
};

} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_PROPERTY_REFINER_HPP
