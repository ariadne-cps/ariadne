/***************************************************************************
 *            lptst.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it Pieter.Collins@cwi.nl
 ****************************************************************************/
/*
 * Based on the linear programming algorithms in PPL-0.8
 *   Copyright (C) 2001-2006 Roberto Bagnara <bagnara@cs.unipr.it>
 */

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
 
/*! \file lptst.h
 *  \brief Linear programming test.
 */

#ifndef ARIADNE_LPTST_H
#define ARIADNE_LPTST_H

namespace Ariadne { 
  namespace LinearProgramming {

    
    /*! \ingroup LinearProgramming
     *  \brief Test for \f$Ax\leq b\f$
     *
     * \return True if a proposed solution is surely valid; false if it is surely not, indeterminate otherwise.
     * \param A is an \f$m\times n\f$ matrix.
     * \param b is an \f$m\f$-vector.
     * \param x is an \f$n\f$-vector.
     */
    template<class R, class AP>
    tribool 
    lptst(const LinearAlgebra::Matrix<R>& A, 
          const LinearAlgebra::Vector<R>& b, 
          const LinearAlgebra::Vector<AP>& x);


    /*! \ingroup LinearProgramming
     *  \brief Test for \f$Ax\leq b\f$
     *
     * \return True if a proposed solution is surely optimal; false if it is surely not, indeterminate otherwise.
     * \param A is an \f$m\times n\f$ matrix.
     * \param b is an \f$m\f$-vector.
     * \param c is an \f$n\f$-vector.
     * \param x is an \f$m\f$-vector.
     * \param y is an \f$n\f$-vector.
     */
    template<class R, class AP>
    tribool 
    lptstopt(const LinearAlgebra::Matrix<R>& A, 
             const LinearAlgebra::Vector<R>& b, 
             const LinearAlgebra::Vector<R>& c, 
             const LinearAlgebra::Vector<AP>& x, 
             const LinearAlgebra::Vector<AP>& y);


  } // namespace LinearProgramming
} // namespace Ariadne

#include "lptst.template.h"

#endif /* ARIADNE_LPTST_H */
