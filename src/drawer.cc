/***************************************************************************
 *            drawer.cc
 *
 *  Copyright 2011-12  Pieter Collins
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

#include "drawer.h"

#include "macros.h"
#include "logging.h"
#include "polynomial.h"
#include "function.h"
#include "taylor_function.h"
#include "procedure.h"
#include "function_set.h"
#include "affine_set.h"
#include "paving_interface.h"
#include "grid_set.h"
#include "nonlinear_programming.h"
#include "constraint_solver.h"
#include "paver.h"
#include "affine_set.h"

#include "graphics_interface.h"

namespace Ariadne {

Void BoxDrawer::draw(CanvasInterface& cnvs, const Projection2d& proj, const ValidatedConstrainedImageSet& set) { ARIADNE_NOT_IMPLEMENTED; }
Void SubdivisionDrawer::draw(CanvasInterface& cnvs, const Projection2d& proj, const ValidatedConstrainedImageSet& set) { ARIADNE_NOT_IMPLEMENTED; }
Void GridDrawer::draw(CanvasInterface& cnvs, const Projection2d& proj, const ValidatedConstrainedImageSet& set) { ARIADNE_NOT_IMPLEMENTED; }

Void affine_draw(CanvasInterface& cnvs, const Projection2d& proj, const ValidatedConstrainedImageSet& set, Int depth)
{
    if( depth==0) {
        set.affine_over_approximation().draw(cnvs,proj);
    } else {
        Pair<ValidatedConstrainedImageSet,ValidatedConstrainedImageSet> split=set.split();
        affine_draw(cnvs,proj,split.first,depth-1u);
        affine_draw(cnvs,proj,split.second,depth-1u);
    }
}

Void AffineDrawer::draw(CanvasInterface& cnvs, const Projection2d& proj, const ValidatedConstrainedImageSet& set)
{
    affine_draw(cnvs,proj,set,this->_depth);
}

} // namespace Ariadne


