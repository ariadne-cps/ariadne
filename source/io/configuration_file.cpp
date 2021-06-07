/***************************************************************************
 *            io/configuration_file.cpp
 *
 *  Copyright  2020-21  Luca Geretti
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

#include "config.hpp"
#include "configuration_file.hpp"
#include "utility/string.hpp"

namespace Ariadne {

static String CONFIGURATION_FILE_NAME = ".ariadne.yml";

Void ConfigurationFile::load() const {

}

#ifdef HAVE_YAML_H

#endif // HAVE_YAML_H

} // namespace Ariadne

