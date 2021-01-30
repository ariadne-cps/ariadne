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

//! \brief Interface for refining a value in the search space (cast to double)
class ConfigurationPropertyRefinerInterface {
  public:
    //! \brief Apply the amount to the current value using the refinement rule chosen
    virtual double apply(double amount, double current) const = 0;

    virtual ConfigurationPropertyRefinerInterface* clone() const = 0;
    virtual ~ConfigurationPropertyRefinerInterface() = default;
};

//! \brief Proportional control refiner: next = current*(1-G*amount) with G a gain constant
class ProportionalRefiner : public ConfigurationPropertyRefinerInterface {
  public:
    ProportionalRefiner(double Kp) : _Kp(Kp) { }
    double apply(double amount, double current) const override { return current + _Kp*amount; }
    ConfigurationPropertyRefinerInterface* clone() const override { return new ProportionalRefiner(*this); }
  private:
    double _Kp;
};

/*
template<> struct PIDRefiner<ExactDouble> : public ConfigurationPropertyRefinerInterface<ExactDouble> {
    PIDRefiner(double Kp, double Ki, double Kd) : _gain(gain) { }
    ExactDouble apply(double amount, ExactDouble const& current) const override {
        return ExactDouble(current.get_d()*(1.0 - _gain*amount));
    }
    ConfigurationPropertyRefinerInterface* clone() const override { return new ProportionalRefiner(*this); }
private:
    double _Kp;
    double _Ki;
    double _Kd;
};
*/

} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_PROPERTY_REFINER_HPP
