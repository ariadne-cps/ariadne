/***************************************************************************
 *            io/configuration_file.hpp
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

#ifndef ARIADNE_CONFIGURATION_FILE_HPP
#define ARIADNE_CONFIGURATION_FILE_HPP

#include "utility/typedefs.hpp"

namespace Ariadne {

//! \brief YAML file for configuration preferences
//! \details Preferences loaded from the file will be superseded by CLI preferences.
class ConfigurationFile {
  private:
    ConfigurationFile();
  public:
    ConfigurationFile(ConfigurationFile const&) = delete;
    Void operator=(ConfigurationFile const&) = delete;

    static ConfigurationFile& instance() {
        static ConfigurationFile instance;
        return instance;
    }

    //! \brief Load the preferences.
    Void load() const;
};

} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_FILE_HPP