/***************************************************************************
 *            evaluation/exceptions.h
 *
 *  Copyright  2005-7  Pieter Collins, Alberto Casagrande
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
 
/*! \file evaluation/exceptions.h
 *  \brief Exceptions, error handling and assertions.for the Evaluation module.
 */

#ifndef ARIADNE_EVALUATION_EXCEPTIONS_H
#define ARIADNE_EVALUATION_EXCEPTIONS_H

#include <stdexcept>
#include <iosfwd>

namespace Ariadne {
  namespace Evaluation {
    
    //@{ \name Exceptions
    /*! \brief %Base class for exceptions in the Evaluation module. */
    struct EvaluationException : public std::runtime_error {
      EvaluationException(const std::string& s) : std::runtime_error(s) { }
    };
    /*! \brief The set appears to cross a constraint non-transversely. */
    class NonTransverseCrossingException : public std::exception { };
    /*! \brief The set appears to cross two or more constraints. */
    class CornerCollisionException : public std::exception { };
    /*! \brief A constraint is crossed during a time interval which does not contain the integration time window. */
    class PartiallyEnabledConstraintException : public std::exception { };

    //@}

  } 
}

#endif /* ARIADNE_EVALUATION_EXCEPTIONS_H */
