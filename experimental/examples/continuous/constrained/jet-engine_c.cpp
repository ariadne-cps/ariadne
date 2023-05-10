/***************************************************************************
 *            jet-engine_c.cpp
 *
 *  Copyright  2023  Luca Geretti
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

#include "ariadne_main.hpp"
#include "jet-engine_c.hpp"

void ariadne_main() {
    auto spec = JET_c();
    auto configuration = get_configuration();
    auto constraints_prescriptions = generate_ellipsoidal_constraints(100,spec,configuration);
    CONCLOG_PRINTLN(frequencies(constraints_prescriptions))

    constrained_execution(spec,configuration,constraints(constraints_prescriptions));
}
