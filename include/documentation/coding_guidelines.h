/***************************************************************************
 *            coding_guidelines.h
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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

/*! 

\file coding_guidelines.h
\brief Standards and guidelines for writing C++ and Python code



\page codingguidelines Coding Guidelines

\section cppcodingstandards C++ Coding Standards
 - Use C++ casts rather than C-style casts. i.e. use static_cast<T>(x) and not (T)x.

 - Remember that identifiers beginning with two underscores __XXX, __xxx or an underscore and a capital letter _Xxx are reserved for the compiler.

 - For most classes
    - Header files xxx.h should contain the interface and Doxygen documentation.
    - Header files xxx.inline.h should contain inline functions.
    - Header files xxx.template.h should contain non-inline template functions for template classes and functions which are not designed to be explicitly instantiated.
    - Header files xxx.code.h should contain non-inline template functions for template classes and functions which are designed to be explicitly instantiated.
    - Source files xxx.cc should contain non-template, non-inline code and explicit template instantiations.

 - For wrappers and classes with only inline functions, the function definitions may be in the class body or in the file xxx.h.

 - To ease template instantiation, non-member functions related to a class should not require explicit instantiation. This can be forced by providing an _instantiate() method in the class which calls all related functions.

<b>C++ Indentation</b>

 - The basic indentation size is two spaces.
 - Tabs must not be used for indentation
 - Namespaces
 - Scope classifiers (public, protected, private) should be indented by half the basic indentation.

*/
