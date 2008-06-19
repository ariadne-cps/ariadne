/***************************************************************************
 *            standard_subdivider.h
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
 
/*! \file standard_subdivider.h
 *  \brief Methods for subdividing basic sets
 */

#ifndef ARIADNE_STANDARD_SUBDIVIDER_H
#define ARIADNE_STANDARD_SUBDIVIDER_H

#include "base/types.h"
#include "base/declarations.h"
#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"
#include "geometry/declarations.h"

#include "geometry/list_set.h"
#include "subdivider_interface.h"

namespace Ariadne {
  

    /*! \brief Methods for subdividing a basic set.
     *  \ingroup Approximation
     */
    template<class ES> class StandardSubdivider
      : public SubdividerInterface<ES> 
    {
      typedef typename ES::real_type R;
     public:
      virtual StandardSubdivider<ES>* clone() const {
        return new StandardSubdivider<ES>(*this); }
      virtual R radius(const ES& es) const {
        return es.radius(); }
      virtual ListSet<ES> split(const ES& es) const {
        return Ariadne::split(es); }
      virtual ListSet<ES> subdivide(const ES& es, const R& r) const {
        ListSet<ES> result; ListSet<ES> working(es); ES set;
        while(working.size()!=0) { set=working.pop(); 
          if(set.radius()<r) { result.adjoin(set); } else { working.adjoin(this->split(set)); } }
        return result; }
    };




  
} // namespace Ariadne

#endif /* ARIADNE_STANDARD_SUBDIVIDER_H */
