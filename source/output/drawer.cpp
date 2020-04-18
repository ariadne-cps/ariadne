/***************************************************************************
 *            output/drawer.cpp
 *
 *  Copyright  2011-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "../function/functional.hpp"
#include "../config.hpp"

#include "../output/drawer.hpp"

#include "../utility/macros.hpp"
#include "../output/logging.hpp"
#include "../geometry/function_set.hpp"
#include "../geometry/affine_set.hpp"
#include "../geometry/grid_paving.hpp"

#include "../output/graphics_interface.hpp"

namespace Ariadne {

Void box_draw(CanvasInterface& cnvs, const Projection2d& proj, const ValidatedConstrainedImageSet& set);
Void affine_draw(CanvasInterface& cnvs, const Projection2d& proj, const ValidatedConstrainedImageSet& set, Nat splittings_remaining);

Pair<Nat,FloatDP> nonlinearity_index_and_error(const ValidatedVectorMultivariateFunction& function, const ExactBoxType& domain);


Void SubdivisionDrawer::draw(CanvasInterface& cnvs, const Projection2d& proj, const ValidatedConstrainedImageSet& set) const { ARIADNE_NOT_IMPLEMENTED; }

OutputStream& BoxDrawer::_write(OutputStream& os) const {
    return os << "BoxDrawer()"; }
OutputStream& AffineDrawer::_write(OutputStream& os) const {
    return os << "AffineDrawer(splittings=" << this->_splittings << ")"; }
OutputStream& EnclosureAffineDrawer::_write(OutputStream& os) const {
    return os << "EnclosureAffineDrawer(accuracy=" << this->_accuracy << ")"; }
OutputStream& SubdivisionDrawer::_write(OutputStream& os) const {
    return os << "SubdivisionDrawer()"; }
OutputStream& GridDrawer::_write(OutputStream& os) const {
    return os << "GridDrawer(fineness=" << this->_fineness << ")"; }

Void BoxDrawer::draw(CanvasInterface& cnvs, const Projection2d& proj, const ValidatedConstrainedImageSet& set) const { box_draw(cnvs,proj,set); }

Void AffineDrawer::draw(CanvasInterface& cnvs, const Projection2d& proj, const ValidatedConstrainedImageSet& set) const
{
    affine_draw(cnvs,proj,set,this->_splittings);
}

Void EnclosureAffineDrawer::draw(CanvasInterface& canvas, const Projection2d& projection, const ValidatedConstrainedImageSet& set) const {
    Nat accuracy = this->_accuracy;

    ARIADNE_ASSERT_MSG(Ariadne::subset(set.reduced_domain(),set.domain()),set);

    // Bound the maximum number of splittings allowed to draw a particular set.
    // Note that this gives rise to possibly 2^MAX_DEPTH split sets!!
    static const Int MAXIMUM_DEPTH = 16;

    // The basic approximation error when plotting with accuracy=0
    static const double BASIC_ERROR = 0.0625;

    const double max_error=BASIC_ERROR/(1<<accuracy);

    // If the reduced domain is empty, then the set is empty; abort
    if(set.reduced_domain().is_empty()) {
        return;
    }

    ValidatedVectorMultivariateFunction fg(2u+set.number_of_constraints(),set.domain());
    fg[0]=set.function()[projection.i];
    fg[1]=set.function()[projection.i];
    for(Nat i=0; i!=set.constraints().size(); ++i) { fg[i+2u]=set.constraints()[i].function(); }
    Projection2d identity(2, 0,1);
//    ValidatedVectorMultivariateFunctionModelDP fg=join(set.state_function(),set.time_function(),set.constraint_function());

    List<ExactBoxType> subdomains;
    List<ExactBoxType> unsplitdomains;
    List<ExactBoxType> splitdomains;
    unsplitdomains.append(set.reduced_domain());
    ExactBoxType splitdomain1,splitdomain2;
    for(Int i=0; i!=MAXIMUM_DEPTH; ++i) {
        //std::cerr<<"i="<<i<<"\nsubdomains="<<subdomains<<"\nunsplitdomains="<<unsplitdomains<<"\n\n";
        for(Nat n=0; n!=unsplitdomains.size(); ++n) {
            Nat k; FloatDP err;
            make_lpair(k,err)=nonlinearity_index_and_error(fg,unsplitdomains[n]);
            //std::cerr<<"  domain="<<unsplitdomains[n]<<" k="<<k<<" err="<<err<<" max_err="<<max_error<<"\n";
            if(k==set.number_of_parameters() || err < max_error) {
                subdomains.append(unsplitdomains[n]);
            } else {
                make_lpair(splitdomain1,splitdomain2)=unsplitdomains[n].split(k);
                splitdomains.append(splitdomain1);
                splitdomains.append(splitdomain2);
            }
        }
        unsplitdomains.swap(splitdomains);
        splitdomains.clear();
        if(unsplitdomains.empty()) { break; }
    }
    subdomains.concatenate(unsplitdomains);
    if(!unsplitdomains.empty()) {
        ARIADNE_WARN("Cannot obtain desired accuracy in drawing "<<set<<" without excessive splitting.");
    }

    for(Nat n=0; n!=subdomains.size(); ++n) {
        try {
            set.restriction(subdomains[n]).affine_over_approximation().draw(canvas,projection);
        } catch(const std::runtime_error& e) {
            ARIADNE_WARN("ErrorTag "<<e.what()<<" in EnclosureAffineDrawer::draw(...) for "<<set<<"\n");
            set.restriction(subdomains[n]).box_draw(canvas,projection);
        }
    }
}

Void GridDrawer::draw(CanvasInterface& canvas, const Projection2d& projection, const ValidatedConstrainedImageSet& set) const {
    set.outer_approximation(Grid(set.dimension()),this->_fineness).draw(canvas,projection);
}

Void box_draw(CanvasInterface& cnvs, const Projection2d& proj, const ValidatedConstrainedImageSet& set)
{
    cast_exact_box(apply(set.function(),set.domain())).draw(cnvs,proj);
}

Void affine_draw(CanvasInterface& cnvs, const Projection2d& proj, const ValidatedConstrainedImageSet& set, Nat splittings_remaining)
{
    if(splittings_remaining==0) {
        set.affine_over_approximation().draw(cnvs,proj);
    } else {
        Pair<ValidatedConstrainedImageSet,ValidatedConstrainedImageSet> split=set.split();
        affine_draw(cnvs,proj,split.first,splittings_remaining-1u);
        affine_draw(cnvs,proj,split.second,splittings_remaining-1u);
    }
}

} // namespace Ariadne


