/***************************************************************************
 *            function/variable.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
 
/*! \file function/variable.h
 *  \brief Named variables of functions
 */
 
#ifndef ARIADNE_FUNCTION_VARIABLE_H
#define ARIADNE_FUNCTION_VARIABLE_H

#include <iostream>

namespace Ariadne {

    struct Variable { std::string name; bool array_flag; uint size; };
    struct FunctionVariable : public Variable { enum Type { OUTPUT=0,INPUT=1,INTERMEDIATE=2,CONSTANT=3 }; Type type; int start; };

    std::ostream& operator<<(std::ostream& os, const Variable& var);
    std::ostream& operator<<(std::ostream& os, const FunctionVariable& var);

} // namespace Ariadne

#endif /* ARIADNE_FUNCTION_VARIABLE_H */
