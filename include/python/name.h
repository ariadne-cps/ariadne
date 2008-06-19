/***************************************************************************
 *            python/name.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file python/name.h
 *  Commonly used inline methods for the Python interface.
 */
 
#ifndef ARIADNE_PYTHON_NAME_H
#define ARIADNE_PYTHON_NAME_H

#include <cstring>
#include <functional>

#include "base/types.h"
#include "numeric/declarations.h"
#include "python/float.h"

namespace Ariadne {
namespace Python {

  template<class R> inline std::string python_name(const std::string& bn);

  template<> inline std::string python_name<bool>(const std::string& bn) { return "Boolean"+bn; }
  template<> inline std::string python_name<index_type>(const std::string& bn) { return "Index"+bn; }
  template<> inline std::string python_name<size_type>(const std::string& bn) { return "Size"+bn; }
  template<> inline std::string python_name<Integer>(const std::string& bn) { return "Z"+bn; }
  template<> inline std::string python_name<Rational>(const std::string& bn) { return "Q"+bn; }

  #if defined PYTHON_FLOAT64
  template<> inline std::string python_name<Float64>(const std::string& bn) { return ""+bn; }
  template<> inline std::string python_name<FloatMP>(const std::string& bn) { return "MPF"+bn; }
  #elif defined PYTHON_FLOATMP
  template<> inline std::string python_name<Float64>(const std::string& bn) { return "F64"+bn; }
  template<> inline std::string python_name<FloatMP>(const std::string& bn) { return ""+bn; }
  #else
  template<> inline std::string python_name<Float64>(const std::string& bn) { return "F64"+bn; }
  template<> inline std::string python_name<FloatMP>(const std::string& bn) { return "MPF"+bn; }
  #endif


} // namespace Python
} // namespace Ariadne

#endif /* ARIADNE_PYTHON_NAME_H */
