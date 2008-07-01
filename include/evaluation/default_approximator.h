/***************************************************************************
 *            default_approximator.h
 *
 *  Copyright  2008  Pieter Collins
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
 
/*! \file default_evaluators.h
 *  \brief Default evaluator classes.
 */

#ifndef ARIADNE_DEFAULT_APPROXIMATOR_H
#define ARIADNE_DEFAULT_APPROXIMATOR_H

#include "evaluation/standard_approximator.h"
#include "evaluation/evolution_parameters.h"

namespace Ariadne {
  


template<class R, class ES>
ApproximatorInterface<GridApproximationScheme<R>,ES>*
default_approximator(const Box<R>& bs, const ES& es, const EvolutionParameters<R>& p)
{
  return new StandardApproximator<ES>();
}


template<class R, class ES>
ApproximatorInterface<GridApproximationScheme<R>,ES>*
default_approximator(const Box<R>& bs, const ES& es)
{
  EvolutionParameters<R> p;
  return default_approximator(bs,es,p);
}


template<class BS, class ES>
ApproximatorInterface<GridApproximationScheme<typename BS::real_type>,ES>*
default_approximator()
{
  BS* bs=0; ES* es=0;
  return default_approximator(*bs,*es);
}

template<class ES>
ApproximatorInterface<GridApproximationScheme<typename ES::real_type>,ES>*
default_approximator()
{
  Box<typename ES::real_type>* bs=0; ES* es=0;
  return default_approximator(*bs,*es);
}



} // namespace Ariadne



#endif /* ARIADNE_DEFAULT_APPROXIMATOR_H */
