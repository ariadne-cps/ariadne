/***************************************************************************
 *            noisy-benchmark.cpp
 *
 *  Copyright  2008-18 Luca Geretti
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

#include "higgins-selkov.hpp"
#include "chemical-reactor.hpp"
#include "lotka-volterra.hpp"
#include "jet-engine.hpp"
#include "pi-controller.hpp"
#include "jerk21.hpp"
#include "lorenz-attractor.hpp"
#include "rossler-attractor.hpp"
#include "jerk16.hpp"
#include "dc-dc.hpp"
#include "harmonic-oscillator.hpp"

#include "noisy-utilities.hpp"

using namespace Ariadne;


int main()
{
    List<SystemType> systems = {HS(),CR(),LV(),JE(),PI(),J21(),LA(),RA(),J16(),DC(),HO()};

    for (SystemType s : systems) {
        std::cout << std::get<0>(s) << std::endl;
        run_noisy_system(s);
    }
}
