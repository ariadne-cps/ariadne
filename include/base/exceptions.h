/***************************************************************************
 *            exceptions.h
 *
 *  Copyright  2005-6  Pieter Collins, Alberto Casagrande
 *  Email  Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
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
 
/*! \file exceptions.h
 *  \brief Exceptions, error handling and assertions.
 */

#ifndef _ARIADNE_EXCEPTION_H
#define _ARIADNE_EXCEPTION_H

#include <stdexcept>
#include <iosfwd>

namespace Ariadne {


class NotImplemented : public std::logic_error {
 public:
  NotImplemented(const std::string& str) : std::logic_error(str) { }
};


class invalid_input : public std::runtime_error {
 public:
  invalid_input(const std::string& str) : std::runtime_error(str) { }
};



namespace Evaluation {

/*! \brief %Base class for exceptions in the Evaluation module. */
class EvaluationException
  : public std::exception 
{
 public:
  EvaluationException(const std::string& s) : _what(s) { }
  ~EvaluationException() throw () { }
  const char* what() const throw () { return this->_what.c_str(); }
 private:
  std::string _what;
};

} // namespace Evaluation

} // namespace Ariadne

#endif /* _ARIADNE_EXCEPTION_H */
