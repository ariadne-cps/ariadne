/***************************************************************************
 *            reducer_interface.h
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
 
/*! \file reducer_interface.h
 *  \brief Interface for reducing the description of an enclosure set
 */

#ifndef ARIADNE_REDUCER_INTERFACE_H
#define ARIADNE_REDUCER_INTERFACE_H

namespace Ariadne {
  

    /*! \brief Interface for over-approximating enclosure sets by sets with a reduced (simpler) description. 
     *  \ingroup EvaluatorInterfaces \ingroup Approximators
     */
    template<class ES> 
    class ReducerInterface 
    { 
     public:
      /*! \brief Virtual destructor. */
      virtual ~ReducerInterface() { }
      /*! \brief Create a dynamically-allocated copy. */
      virtual ReducerInterface<ES>* clone() const = 0;
      /*! \brief Computes an over-approximation of the encloser \a es of a reduced complexity. */
      virtual ES over_approximate(const ES& es) const = 0;
    };

  
} // namespace Ariadne

#endif /* ARIADNE_REDUCER_INTERFACE_H */
