/***************************************************************************
 *            coding_guidelines.h
 *
 *  Copyright  2007  Pieter Collins
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

/*! 

\file coding_guidelines.h
\brief Standards and guidelines for writing C++ and Python code



\page codingguidelines Coding Guidelines

\section cppcodingstandards C++ Coding Standards
 - Use C++ casts rather than C-style casts. i.e. use static_cast<T>(x) and not (T)x.

 - Remember that identifiers beginning with two underscores __XXX, __xxx or an underscore and a capital letter _Xxx are reserved for the compiler.

 - Each class should be in a single header file, which should include strongly related classes and methods.
    - Header files xxx.h should contain the interface, inline functions, templated code which is not explicitly instantiated and Doxygen documentation.
    - Source files xxx.cc should contain ordinary non-inline source code, templated code which is meant to be explicitly instantiated, internal code only used in the source file and explicit template instantiations.

 - For wrappers and classes with only inline functions, the function definitions may be in the class body or in the file xxx.h.

 - To ease template instantiation, non-member functions related to a class should not require explicit instantiation. 

 - The use of a private static _instantiate() method to enforce template instantiation does not work with all compilers and all optimization settings, and must not be used.

<b>C++ Indentation</b>

 - The basic indentation size is four spaces.
 - Tabs must not be used for indentation
 - Namespaces should not be indented
 - Scope classifiers (public, protected, private) should be indented by half the basic indentation.

<b>Common mistakes</b>
 - A class static constant should be declared inside the class and initialised outside the class body. This initialisation should go in a source code file (.cc or .code.h or .template.h). An exception is an integer constant, which may be initialised in the class definition.
  
\section pythoncodingstandards Python Coding Standards
 - Every executable Python code file should begin with <c>#!/usr/bin/python</c>. This is the proper location for the Python interpreter following the Linux <a href="http://www.pathname.com/fhs/">Filesystem Hierarchy Standard</a> (FHS).
    

*/
