/***************************************************************************
 *            linear_progreamming/exceptions.h
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
 
/*! \file linear_algebra/exceptions.h
 *  \brief Exceptions, error handling and assertions for the Linear Algebra module.
 */

#ifndef ARIADNE_LINEAR_PROGRAMMING_EXCEPTIONS_H
#define ARIADNE_LINEAR_PROGRAMMING_EXCEPTIONS_H

#include <stdexcept>
#include <iosfwd>

#include "../throw.h"
#include "../base/types.h"
#include "../base/exceptions.h"
#include "../linear_algebra/exceptions.h"

namespace Ariadne {
  namespace LinearProgramming {
    
    //@{ \name Exceptions
    
    /*! \brief The problem is unbound. */
    struct UnboundProblem : public std::runtime_error {
      UnboundProblem(const std::string& str) : std::runtime_error(str) { }
    };

    /*! \brief The problem has no feasible origin. */
    struct InfeasibleOrigin : public std::runtime_error {
      InfeasibleOrigin(const std::string& str) : std::runtime_error(str) { }
    };

    /*! \brief The problem has no feasible solution. */
    struct InfeasibleSolution : public std::runtime_error {
      InfeasibleSolution(const std::string& str) : std::runtime_error(str) { }
    };

    //@}

  }
}


#endif /* ARIADNE_LINEAR_PROGRAMMING_EXCEPTIONS_H */
