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
    //! \brief Apply the \a error with a given \a progress measure to the \a current value using the refinement chosen
    //! \details It may change the internal state of the object.
    virtual double apply(double error, double progress, double current) = 0;

    virtual ConfigurationPropertyRefinerInterface* clone() const = 0;
    virtual ~ConfigurationPropertyRefinerInterface() = default;
};

//! \brief Proportional control refiner
class ProportionalRefiner : public ConfigurationPropertyRefinerInterface {
  public:
    ProportionalRefiner(double Kp) : _Kp(Kp) { }
    double apply(double error, double step, double current) override { return current + _Kp*error; }
    ConfigurationPropertyRefinerInterface* clone() const override { return new ProportionalRefiner(*this); }
  private:
    double _Kp;
};

//! \brief Proportional-Integrative-Derivative control refiner
class PIDRefiner : public ConfigurationPropertyRefinerInterface {
  public:
    PIDRefiner(double Kp, double Ki, double Kd) : _Kp(Kp), _Ki(Ki), _Kd(Kd), _error_k_2(0), _error_k_1(0), _step_k_1(1e-20) { }
    double apply(double error_k, double step_k, double current) override {
        if (error_k == std::numeric_limits<double>::min()) return std::numeric_limits<double>::max();
        auto result = current + (_Kp+_Ki*step_k+_Kd/step_k)*error_k - (_Kp+_Kd*(step_k+_step_k_1)/(step_k*_step_k_1))*_error_k_1 + _Kd/_step_k_1*_error_k_2;
        _error_k_2 = _error_k_1;
        _error_k_1 = error_k;
        _step_k_1 = step_k;
        return result;
    }
    ConfigurationPropertyRefinerInterface* clone() const override { return new PIDRefiner(*this); }
  private:
    const double _Kp;
    const double _Ki;
    const double _Kd;
    double _error_k_2;
    double _error_k_1;
    double _step_k_1;
};

} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_PROPERTY_REFINER_HPP
