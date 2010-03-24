/***************************************************************************
 *            taylor_set.cc
 *
 *  Copyright 2008  Pieter Collins
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

#include "macros.h"
#include "exceptions.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"
#include "differential.h"
#include "function.h"

#include "geometry.h"
#include "taylor_model.h"
#include "taylor_function.h"

#include "box.h"
#include "zonotope.h"
#include "polytope.h"

#include "function_set.h"
#include "taylor_set.h"
#include "affine_set.h"

#include "list_set.h"
#include "grid_set.h"

#include "constraint_solver.h"
#include "nonlinear_programming.h"

#include "graphics_interface.h"

#include "discrete_event.h"

namespace Ariadne {

template<class T> std::string str(const T& t) { std::stringstream ss; ss<<t; return ss.str(); }

typedef Vector<Float> FloatVector;
typedef Vector<Interval> IntervalVector;

TaylorImageSet::TaylorImageSet(uint d, uint ng)
    : _models(d,ng)
{
}


TaylorImageSet::TaylorImageSet(const VectorFunction& f, const Vector<Interval>& d)
    : _models(f.result_size())
{
    Vector<TaylorModel> x=TaylorModel::scalings(d);
    this->_models=f.evaluate(x);
}

TaylorImageSet::TaylorImageSet(const Vector<Interval>& bx)
    : _models(TaylorModel::scalings(bx))
{
}


TaylorImageSet::TaylorImageSet(const Vector<TaylorModel>& tvs)
    : _models(tvs)
{
    if(this->_models.size()>0) {
        Vector<Interval> unit_domain(this->_models[0].argument_size(),Interval(-1,+1));
        for(size_t i=0; i!=this->_models.size(); ++i) {
            ARIADNE_ASSERT(this->_models[i].domain()==unit_domain);
        }
    }
}

TaylorImageSet::TaylorImageSet(const Vector< Expansion<Float> >& f, const Vector<Float>& e)
    : _models(f.size())
{
    ARIADNE_ASSERT(f.size()==e.size());
    if(f.size()>0) {
        for(uint i=0; i!=f.size(); ++i) { (*this)[i]=TaylorModel(f[i],e[i]); }
    }
}

TaylorImageSet::TaylorImageSet(uint rs, uint as, uint deg, double x0, ...)
    : _models(rs)
{
    double x=x0;
    va_list args;
    va_start(args,x0);
    for(uint i=0; i!=rs; ++i) {
        (*this)[i]=TaylorModel(as);
        for(MultiIndex j(as); j.degree()<=deg; ++j) {
            if(x!=0.0) { (*this)[i].expansion().append(j,x); }
            x=va_arg(args,double);
        }
        (*this)[i].error()=x;
        x=va_arg(args,double);
        (*this)[i].expansion().cleanup();
        (*this)[i].clean();
    }
    va_end(args);
}

bool
operator==(const TaylorImageSet& ts1,const TaylorImageSet& ts2)
{
    return ts1.models()==ts2.models();
}


void
TaylorImageSet::set_accuracy(shared_ptr<TaylorModel::Accuracy> acc_ptr)
{
    for(uint i=0; i!=this->dimension(); ++i) {
        this->_models[i].set_accuracy(acc_ptr);
    }
}

shared_ptr<TaylorModel::Accuracy>
TaylorImageSet::accuracy_ptr() const
{
    return this->_models[0].accuracy_ptr();
}


Vector<Float>
TaylorImageSet::centre() const
{
    Vector<Float> result(this->dimension());
    for(uint i=0; i!=this->dimension(); ++i) {
        result[i]=this->_models[i].value();
    }
    return result;
}



Float
TaylorImageSet::radius() const
{
    Float result=0.0;
    for(uint i=0; i!=this->dimension(); ++i) {
        result=max(result,this->models()[i].range().radius());
    }
    return result;
}


tribool
TaylorImageSet::disjoint(const Box& bx) const
{
    ARIADNE_ASSERT_MSG(this->dimension()==bx.dimension(),"Dimension of "<<*this<<" is different from dimension of "<<bx);

    // Expand the bounding box by the error term
    Box ebx(bx);
    for(uint i=0; i!=bx.dimension(); ++i) {
        ebx[i]+=Interval(-this->_models[i].error(),+this->_models[i].error());
    }

    // Copy the box eliminating error terms
    TaylorImageSet efts(*this);
    for(uint i=0; i!=efts.dimension(); ++i) {
        efts._models[i].clobber();
    }

    Box bbx=efts.bounding_box();
    if(bbx.disjoint(ebx)) {
        return true;
    } else if(bbx.subset(ebx)) {
        return false;
    } else if(ebx.contains(efts.centre())) {
        return false;
    } else if(bbx.radius() * 4 < ebx.radius()) {
        return indeterminate;
    } else {
        std::pair<TaylorImageSet,TaylorImageSet> split_efts=this->split();
        return split_efts.first.disjoint(ebx) && split_efts.second.disjoint(ebx);
    }
}



tribool
TaylorImageSet::overlaps(const Box& bx) const
{
    return !this->disjoint(bx);
}


tribool
TaylorImageSet::inside(const Box& bx) const
{
    Vector<Interval> bb=this->bounding_box();
    return Ariadne::inside(bb,bx) || indeterminate;
}


Box
TaylorImageSet::bounding_box() const
{
    Box r(this->dimension());
    for(uint i=0; i!=this->dimension(); ++i) {
        r[i]=(*this)[i].range();
    }
    return r.bounding_box();
}


TaylorImageSet
TaylorImageSet::linearise() const
{
    TaylorImageSet res(*this);
    for(uint i=0; i!=this->dimension(); ++i) {
        res[i].truncate(1);
    }
    return res;
}


void
_discretise_step(ListSet<Box>& result, const Vector<TaylorModel>& models, const Vector<Interval>& errors, const Box& domain, const Float& eps)
{
    Box range=evaluate(models,domain);
    if(range.radius()<eps) {
        result.adjoin(range+errors);
    } else {
        std::pair<Box,Box> subdomains=split(domain);
        _discretise_step(result,models,errors,subdomains.first,eps);
        _discretise_step(result,models,errors,subdomains.second,eps);
    }
}

ListSet<Box>
TaylorImageSet::discretise(const Float& eps) const
{
    ListSet<Box> result;
    Box domain=this->domain();
    Vector<Interval> errors(this->dimension());
    Vector<TaylorModel> models=this->models();
    for(uint i=0; i!=errors.size(); ++i) {
        errors[i]=models[i].error()*Interval(-1,+1);
        models[i].set_error(0.0);
    }
    _discretise_step(result,models,errors,domain,eps);
    return result;
}

GridTreeSet
TaylorImageSet::discretise(const Grid& g, uint d) const
{
    GridTreeSet gts(g);
    return this->discretise(gts,d);
}

GridTreeSet&
TaylorImageSet::discretise(GridTreeSet& gts, uint d) const
{
    adjoin_outer_approximation(gts,*this,d);
    gts.recombine();
    return gts;
}


TaylorModel apply(const ScalarFunction& f, const TaylorImageSet& ts)
{
    return f.evaluate(ts.models());
}

TaylorImageSet apply(const VectorFunction& f, const TaylorImageSet& ts)
{
    return f.evaluate(ts.models());
}

TaylorModel apply(const ScalarTaylorFunction& tf, const TaylorImageSet& ts)
{
    ARIADNE_ASSERT_MSG(subset(ts.range(),tf.domain()),"tf="<<tf<<" ts="<<ts);
    return unchecked_compose(tf.model(),unscale(ts.models(),tf.domain()));
}

TaylorImageSet apply(const VectorTaylorFunction& tf, const TaylorImageSet& ts)
{
    ARIADNE_ASSERT_MSG(possibly(subset(ts.range(),tf.domain())),
        std::setprecision(18)<<"\n  tf="<<tf<<"\n  ts="<<ts<<"\n  ts.range() ="<<ts.range()<<"\n  tf.domain()="<<tf.domain());
    return unchecked_compose(tf.models(),unscale(ts.models(),tf.domain()));
}

TaylorModel unchecked_apply(const ScalarTaylorFunction& tf, const TaylorImageSet& ts)
{
    return unchecked_compose(tf.model(),unscale(ts.models(),tf.domain()));
}

TaylorImageSet unchecked_apply(const VectorTaylorFunction& tf, const TaylorImageSet& ts)
{
    return unchecked_compose(tf.models(),unscale(ts.models(),tf.domain()));
}


Zonotope
zonotope(const TaylorImageSet& ts)
{
    uint d=ts.dimension();
    uint ng=ts.generators_size();
    Vector<Float> c(d);
    Matrix<Float> G(d,ng);
    Vector<Float> e(d);
    for(uint i=0; i!=d; ++i) {
        c[i]=ts[i].value();
        for(uint j=0; j!=ng; ++j) {
            G[i][j]=ts[i][MultiIndex::unit(ng,j)];
        }
        e[i]=ts[i].error();
    }

    set_rounding_upward();
    for(uint i=0; i!=d; ++i) {
        for(TaylorModel::const_iterator iter=ts[i].begin(); iter!=ts[i].end(); ++iter) {
            if(iter->key().degree()>=2) {
                e[i]+=abs(iter->data());
            }
        }
    }
    set_rounding_to_nearest();
    return Zonotope(c,G,e);
}




pair<TaylorImageSet,TaylorImageSet>
TaylorImageSet::split(uint j) const
{
    pair< Vector<TaylorModel>, Vector<TaylorModel> > s=Ariadne::split(this->_models,j);
    return make_pair( TaylorImageSet(s.first),
                      TaylorImageSet(s.second) );
}

pair<TaylorImageSet,TaylorImageSet>
TaylorImageSet::split() const
{
    uint d=this->dimension();
    uint ng=this->generators_size();
    Float rmax=0.0;
    uint armax=0;
    for(uint j=0; j!=this->generators_size(); ++j) {
        Float r=0.0;
        for(uint i=0; i!=this->dimension(); ++i) {
            r+=abs(this->models()[i][MultiIndex::unit(ng,j)]);
        }
        if(r>rmax) {
            rmax=r;
            armax=j;
        }
    }

    Float emax=0.0;
    uint aemax=0;
    for(uint i=0; i!=d; ++i) {
        Float e=this->models()[i].error();
        if(e>emax) {
            emax=e;
            aemax=i;
        }
    }

    if(emax>rmax) {
        pair<TaylorImageSet,TaylorImageSet> result=make_pair(*this,*this);
        TaylorModel& model1=result.first._models[aemax];
        TaylorModel& model2=result.second._models[aemax];
        model1.set_error(emax/2);
        model1-=emax/2;
        model2.set_error(emax/2);
        model2+=emax/2;
        return result;
    } else {
        return this->split(armax);
    }
}



void
_adjoin_outer_approximation1(GridTreeSet& grid_set,
                             const Vector<TaylorModel>& models, const Vector<Interval>& errors,
                             const Box& domain, Float eps, uint depth)
{
    uint d=models.size();

    Box range=evaluate(models,domain);
    if(range.radius()<eps) {
        grid_set.adjoin_over_approximation(range+errors,depth);
    } else {
        Box domain1(d),domain2(d);
        make_lpair(domain1,domain2)=split(domain);
        _adjoin_outer_approximation1(grid_set,models,errors,domain1,eps,depth);
        _adjoin_outer_approximation1(grid_set,models,errors,domain2,eps,depth);
    }
}


void
adjoin_outer_approximation1(GridTreeSet& grid_set, const TaylorImageSet& set, uint depth)
{
    Box domain(set.generators_size(),Interval(-1,+1));
    Float eps=1.0/(1<<(depth/set.dimension()));
    Vector<TaylorModel> models=set.models();
    Vector<Interval> errors(models.size());
    for(uint i=0; i!=models.size(); ++i) {
        errors[i]=models[i].error()*Interval(-1,+1);
        models[i].set_error(0.0);
    }
    _adjoin_outer_approximation1(grid_set,models,errors,domain,eps,depth);
    grid_set.recombine();
}

void
_adjoin_outer_approximation2(GridTreeSet& grid_set,
                             const TaylorImageSet& set, const Vector<Interval>& errors,
                             Float eps, uint depth)
{
    Box range=set.range();
    if(range.radius()<eps) {
        grid_set.adjoin_over_approximation(range+errors,depth);
    } else {
        std::pair<TaylorImageSet,TaylorImageSet> subsets=set.split();
        _adjoin_outer_approximation2(grid_set,subsets.first,errors,eps,depth);
        _adjoin_outer_approximation2(grid_set,subsets.second,errors,eps,depth);
    }
}


void
adjoin_outer_approximation2(GridTreeSet& grid_set, const TaylorImageSet& set, uint depth)
{
    Box domain(set.generators_size(),Interval(-1,+1));
    Float eps=1.0/(1<<(depth/set.dimension()));
    TaylorImageSet error_free_set(set);
    Vector<Interval> errors(set.dimension());
    for(uint i=0; i!=set.dimension(); ++i) {
        errors[i]=set.models()[i].error()*Interval(-1,+1);
        const_cast<TaylorModel&>(error_free_set.models()[i]).set_error(0.0);
    }
    _adjoin_outer_approximation2(grid_set,error_free_set,errors,eps,depth);
    grid_set.recombine();
}

void
_adjoin_outer_approximation3(GridTreeSet& grid_set, const TaylorImageSet& set, Box& domain,
			     GridOpenCell cell, uint depth)
{
    // Compute an over-approximation to the set
    Box range=evaluate(set.models(),domain);

    // Find an open cell which is a good over-approximation to the range
    while(true) {
        GridOpenCell subcell = cell.split(false);
        if(!subset(range,subcell.box())) {
            subcell=cell.split(true);
            if(!subset(range,subcell.box())) {
                subcell=cell.split(indeterminate);
                if(!subset(range,subcell.box())) {
                    break;
                }
            }
        }
        cell=subcell;
    }

    if(cell.depth()>=int(depth)) {
        grid_set.adjoin(cell.closure());
        return;
    }

    std::pair<Box,Box> subdomains=split(domain);
    _adjoin_outer_approximation3(grid_set,set,subdomains.first,cell,depth);
    _adjoin_outer_approximation3(grid_set,set,subdomains.second,cell,depth);
    //Box subdomain1=split(domain,left);
    //Box subdomain2=split(domain,right);
    //_adjoin_outer_approximation3(grid_set,set,subdomain1,cell,depth);
    //_adjoin_outer_approximation3(grid_set,set,subdomain2,cell,depth);
}

void
adjoin_outer_approximation3(GridTreeSet& grid_set, const TaylorImageSet& set, uint depth)
{
    TaylorImageSet error_free_set=set.subsume();
    Box domain=error_free_set.domain();
    GridOpenCell cell=GridOpenCell::outer_approximation(error_free_set.bounding_box(),grid_set.grid());
    _adjoin_outer_approximation3(grid_set,error_free_set,domain,cell,depth);
    grid_set.recombine();
}


void
_adjoin_outer_approximation4(GridTreeSet& grid_set, const TaylorImageSet& set,
                             GridOpenCell cell, uint depth)
{
    // Compute an over-approximation to the set
    Box range=set.range();

    // Find an open cell which is a good over-approximation to the range
    while(true) {
        GridOpenCell subcell = cell.split(false);
        if(!subset(range,subcell.box())) {
            subcell=cell.split(true);
            if(!subset(range,subcell.box())) {
                subcell=cell.split(indeterminate);
                if(!subset(range,subcell.box())) {
                    break;
                }
            }
        }
        cell=subcell;
    }

    // If the cell is sufficiently small, adjoin it to the grid set
    if(cell.depth()>=int(depth)) {
        grid_set.adjoin(cell.closure());
        return;
    }

    // Subdivide the TaylorImageSet and try again
    std::pair<TaylorImageSet,TaylorImageSet> subsets=set.split();
    _adjoin_outer_approximation4(grid_set,subsets.first,cell,depth);
    _adjoin_outer_approximation4(grid_set,subsets.second,cell,depth);
}

void
adjoin_outer_approximation4(GridTreeSet& grid_set, const TaylorImageSet& set, uint depth)
{
    GridOpenCell cell=GridOpenCell::outer_approximation(set.bounding_box(),grid_set.grid());
    TaylorImageSet error_free_set=set.subsume();
    _adjoin_outer_approximation4(grid_set,error_free_set,cell,depth);
    grid_set.recombine();
}

// Profiling suggests that adjoin_outer_approximation1 is most efficient.
// adjoin_outer_approximation3 (using subdivision of the set) is better than
// adjoin_outer_approximation2 (using subdivision of the domain).
void
adjoin_outer_approximation(GridTreeSet& grid_set, const TaylorImageSet& set, uint depth)
{
    adjoin_outer_approximation1(grid_set,set,depth);
}


GridTreeSet
discretise1(const TaylorImageSet& ts, const Grid& g, uint d)
{
    GridTreeSet gts(g);
    adjoin_outer_approximation1(gts,ts,d);
    gts.recombine();
    return gts;
}

GridTreeSet
discretise2(const TaylorImageSet& ts, const Grid& g, uint d)
{
    GridTreeSet gts(g);
    adjoin_outer_approximation2(gts,ts,d);
    gts.recombine();
    return gts;
}

GridTreeSet
discretise3(const TaylorImageSet& ts, const Grid& g, uint d)
{
    GridTreeSet gts(g);
    adjoin_outer_approximation3(gts,ts,d);
    gts.recombine();
    return gts;
}

GridTreeSet
discretise4(const TaylorImageSet& ts, const Grid& g, uint d)
{
    GridTreeSet gts(g);
    adjoin_outer_approximation4(gts,ts,d);
    gts.recombine();
    return gts;
}


Matrix<Float>
TaylorImageSet::jacobian() const
{
    Matrix<Float> J(this->dimension(),this->generators_size());
    for(uint i=0; i!=this->dimension(); ++i) {
        for(uint j=0; j!=this->generators_size(); ++j) {
            J[i][j]=this->models()[i][MultiIndex::unit(this->generators_size(),j)];
        }
    }
    return J;
}

TaylorImageSet
TaylorImageSet::subsume() const
{
    return this->subsume(0.0);
}


TaylorImageSet
TaylorImageSet::subsume(double eps) const
{
    uint d=this->dimension();
    uint ng=this->generators_size();

    // Compute the number of terms whose error is greater the eps
    uint ne=0;
    for(uint i=0; i!=d; ++i) {
        if(this->models()[i].error()>eps) { ++ne; }
    }

    TaylorImageSet result(embed(this->models(),ne));
    uint k=ng; // The in independent variable corresponding to the error term
    for(uint i=0; i!=d; ++i) {
        if(result._models[i].error()>eps) {
            MultiIndex a=MultiIndex::unit(ng+ne,k);
            result._models[i][a]=result._models[i].error();
            result._models[i].set_error(0.0);
            ++k;
        }
    }
    return result;
}


GridTreeSet
outer_approximation(const TaylorImageSet& set, const Grid& grid, uint depth)
{
    ARIADNE_ASSERT(set.dimension()==grid.dimension());
    GridTreeSet grid_set(grid);
    adjoin_outer_approximation(grid_set,set,depth);
    return grid_set;
}



TaylorImageSet
TaylorImageSet::recondition() const
{
    Matrix<Float> T=triangular_multiplier(this->jacobian());
    Vector<TaylorModel> scal(T.row_size(),T.column_size());
    for(uint i=0; i!=T.row_size(); ++i) {
        for(uint j=0; j!=T.column_size(); ++j) {
            scal[i][MultiIndex::unit(this->generators_size(),j)]=T[i][j];
        }
    }
    return TaylorImageSet(compose(this->models(),scal));
}

std::ostream&
TaylorImageSet::write(std::ostream& os) const
{
    return os << "TaylorImageSet(" << this->_models << ")";
    os << "TaylorImageSet(\n";
    os << "  dimension=" << this->dimension() << ",\n" << std::flush;
    os << "  models=" << this->_models << ",\n" << std::flush;
    os << ")\n";
    return os;
}

Vector<Interval> evaluate(const TaylorImageSet& ts, const Vector<Interval>& iv) {
    Vector<Interval> r(ts.dimension());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=ts[i].evaluate(iv);
    }
    return r;
}

Vector<Interval> evaluate(const TaylorImageSet& ts, const Vector<Float>& v) {
    Vector<Interval> r(ts.dimension());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=ts[i].evaluate(v);
    }
    return r;
}

Vector<Interval> error(const TaylorImageSet& ts) {
    Vector<Interval> r(ts.dimension());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=ts[i].error();
    }
    return r;
}


void box_draw(CanvasInterface& fig, const TaylorImageSet& ts) {
    ts.bounding_box().draw(fig);
}

void boxes_draw(CanvasInterface& fig, const TaylorImageSet& ts) {
    static const double resolution=1.0/8;
    double line_width=fig.get_line_width();
    fig.set_line_width(0.0);
    //fig.set_line_style(false);
    //fig.set_line_colour(fig.get_fill_colour());
    ListSet<Box> boxes=ts.discretise(resolution);
    for(uint i=0; i!=boxes.size(); ++i) { boxes[i].draw(fig); }
    //fig.set_line_colour(black);
    fig.set_line_width(line_width);
    //fig.set_line_style(true);
}

void zonotopes_draw_subdivide(ListSet<Zonotope>& zonotopes, const TaylorImageSet& set, const double resolution)
{
    if(set.bounding_box().radius()>=resolution) {
        std::pair<TaylorImageSet,TaylorImageSet> subdivisions=set.split();
        zonotopes_draw_subdivide(zonotopes,subdivisions.first,resolution);
        zonotopes_draw_subdivide(zonotopes,subdivisions.second,resolution);
    } else {
        zonotopes.adjoin(zonotope(set));
    }
}

void zonotopes_draw(CanvasInterface& fig, const TaylorImageSet& ts) {
    static const double resolution=1.0/8;
    //double old_line_width=fig.get_line_width();
    //Colour old_line_colour=fig.get_line_colour();
    //fig.set_line_width(0.0);
    //fig.set_line_style(false);
    //fig.set_line_colour(fig.get_fill_colour());
    ListSet<Zonotope> zonotopes;
    zonotopes_draw_subdivide(zonotopes,ts,resolution);
    for(uint i=0; i!=zonotopes.size(); ++i) { zonotopes[i].draw(fig); }
    //fig.set_line_colour(old_line_colour);
    //fig.set_line_width(old_line_width);
    //fig.set_line_style(true);
}

void affine_draw(CanvasInterface& fig, const TaylorImageSet& ts) {
    zonotope(ts).draw(fig);
}

void curve_draw(CanvasInterface& fig, const TaylorImageSet& ts) {
    assert(ts.dimension()==2 && ts.generators_size()==2);
    std::vector<Point> pts;
    Vector<Float> s(2);
    s[0]=-1; s[1]=-1;
    for(uint i=0; i!=16; ++i) {
        pts.push_back(midpoint(evaluate(ts,s)));
        s[0]+=1./8;
    }
    for(uint i=0; i!=16; ++i) {
        pts.push_back(midpoint(evaluate(ts,s)));
        s[1]+=1./8;
    }
    for(uint i=0; i!=16; ++i) {
        pts.push_back(midpoint(evaluate(ts,s)));
        s[0]-=1./8;
    }
    for(uint i=0; i!=16; ++i) {
        pts.push_back(midpoint(evaluate(ts,s)));
        s[1]-=1./8;
    }
    fig.move_to(pts[0][0],pts[0][1]);
    for(uint i=0; i!=pts.size(); ++i) {
        fig.line_to(pts[i][0],pts[i][1]);
    }
}


void grid_draw(CanvasInterface& fig, const TaylorImageSet& ts)
{
    uint depth=12;
    Float rad=1./8;
    GridTreeSet gts(Grid(2));
    for(Float i=-1; i!=+1; i+=rad) {
        for(Float j=-1; j!=+1; j+=rad) {
            Vector<Interval> pt(2); pt[0]=Interval(i,i+rad); pt[1]=Interval(j,j+rad);
            gts.adjoin_outer_approximation(Box(evaluate(ts.models(),pt)),depth);
        }
    }
    gts.recombine();
    gts.draw(fig);
}

void standard_draw(CanvasInterface& fig, const TaylorImageSet& ts) {
    affine_draw(fig,ts);
    //box_draw(fig,ts);
}

void TaylorImageSet::draw(CanvasInterface& fig) const {
    Ariadne::standard_draw(fig,*this);
}






} // namespace Ariadne


#include "config.h"

#ifdef HAVE_CAIRO_H

#include <cairo/cairo.h>

namespace Ariadne {

void plot(const char* filename, const Box& bbx, const TaylorImageSet& set)
{
    Box bb=set.bounding_box();

    //cairo_format_t format=CAIRO_FORMAT_ARGB32;
    cairo_format_t format=CAIRO_FORMAT_RGB24;

    int width=256;
    int height=256;
    int stride=cairo_format_stride_for_width(format,width);

    int* data = new int[width*height];

    int white=0x00FFFFFF;
    int black=0x00000000;
    //int red=0x00FF0000;
    //int green=0x0000FF00;
    int blue=0x000000FF;
    //int background=0xFFFFFFFF;
    //int foreground=0x00FF00FF;

    int i,j;
    for(i=0; i!=width; ++i) {
        for(j=0; j!=height; ++j) {
            data[(i*height+j)]=white;
        }
    }

    Vector<Float> v(2);
    Vector<Float> w(2);
    double rad=1.0/256;
    for(double x=-1; x< +1; x+=rad) {
        for(double y=-1; y< +1; y+=rad) {
            v[0]=x; v[1]=y;
            w=midpoint(evaluate(set,v));
            i=((w[0]-bb[0].lower())/2/bb[0].radius())*(width-1);
            assert(0<=i); assert(i<=width); //assert(i<width);
            j=(height-1)-((w[1]-bb[1].lower())/2/bb[1].radius())*(height-1);
            assert(0<=j); assert(j<=height); //assert(j<height);
            data[j*height+i]=blue;
        }
    }

    for(double x=-1; x<=+1; x+=rad) {
        for(double y=-1; y<=+1; y+=2) {
            v[0]=x; v[1]=y;
            w=midpoint(evaluate(set,v));
            i=((w[0]-bb[0].lower())/2/bb[0].radius())*(width-1);
            j=(height-1)-((w[1]-bb[1].lower())/2/bb[1].radius())*(height-1);
            data[j*height+i]=black;
            v[0]=y; v[1]=x;
            w=midpoint(evaluate(set,v));
            i=((w[0]-bb[0].lower())/2/bb[0].radius())*(width-1);
            j=(height-1)-((w[1]-bb[1].lower())/2/bb[1].radius())*(height-1);
            data[j*height+i]=black;
        }
    }


    cairo_surface_t* surface=cairo_image_surface_create_for_data((unsigned char*)data, format, width, height, stride);
    //cairo_t* cr = cairo_create (surface);
    cairo_surface_write_to_png (surface, filename);
    cairo_surface_destroy(surface);

}


TaylorModel& operator-=(TaylorModel& tm, const MultiIndex& a) {
    for(TaylorModel::iterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        iter->key()-=a;
    }
    return tm;
}

// TODO: Make more efficient
inline void assign_all_but_last(MultiIndex& r, const MultiIndex& a) {
    for(uint i=0; i!=r.size(); ++i) { r[i]=a[i]; }
}


void TaylorConstrainedImageSet::_check() const {
    ARIADNE_ASSERT_MSG(this->_function.argument_size()==this->domain().size(),*this);
    for(List<ScalarTaylorFunction>::const_iterator iter=this->_constraints.begin(); iter!=this->_constraints.end(); ++iter) {
        ARIADNE_ASSERT_MSG(iter->argument_size()==this->domain().size(),*this);
    }
    for(List<ScalarTaylorFunction>::const_iterator iter=this->_equations.begin(); iter!=this->_equations.end(); ++iter) {
        ARIADNE_ASSERT_MSG(iter->argument_size()==this->domain().size(),*this);
    }
}

// FIXME: What if solving for constraint leaves domain?
void TaylorConstrainedImageSet::_solve_zero_constraints() {
    this->_check();
    for(List<ScalarTaylorFunction>::iterator iter=this->_equations.begin(); iter!=this->_equations.end(); ) {
        const Vector<Interval>& domain=this->domain();
        const TaylorModel& model=iter->model();
        const uint k=model.argument_size()-1u;
        TaylorModel zeroth_order(k);
        TaylorModel first_order(k);
        bool is_zeroth_order=true;
        bool is_first_order=true;
        MultiIndex r(k);
        // Try linear approach in last coefficient
        for(TaylorModel::const_iterator tmiter=model.begin(); tmiter!=model.end(); ++tmiter) {
            if(tmiter->key()[k]==0) {
                assign_all_but_last(r,tmiter->key());
                zeroth_order.expansion().append(r,tmiter->data());
            } else if(tmiter->key()[k]==1) {
                is_zeroth_order=false;
                assign_all_but_last(r,tmiter->key());
                first_order.expansion().append(r,tmiter->data());
            } else {
                is_first_order=false; break;
            }
        }
        if(is_first_order && !is_zeroth_order) {
            const Vector<Interval> new_domain=project(domain,range(0,k));
            TaylorModel substitution_model=-zeroth_order/first_order;
            this->_function=VectorTaylorFunction(new_domain,Ariadne::substitute(this->_function.models(),k,substitution_model));
            for(List<ScalarTaylorFunction>::iterator constraint_iter=this->_constraints.begin();
                    constraint_iter!=this->_constraints.end(); ++constraint_iter) {
                ScalarTaylorFunction& constraint=*constraint_iter;
                constraint=ScalarTaylorFunction(new_domain,Ariadne::substitute(constraint.model(),k,substitution_model));
            }
            for(List<ScalarTaylorFunction>::iterator constraint_iter=this->_equations.begin();
                    constraint_iter!=this->_equations.end(); ++constraint_iter) {
                ScalarTaylorFunction& constraint=*constraint_iter;
                constraint=ScalarTaylorFunction(new_domain,Ariadne::substitute(constraint.model(),k,substitution_model));
            }
            // Since we are using an std::vector, assign iterator to next element
            iter=this->_equations.erase(iter);
            this->_check();
        } else {
            ARIADNE_WARN("No method for solving constraint "<<*iter<<" currently implemented.");
            ++iter;
        }
    }
}


TaylorConstrainedImageSet::TaylorConstrainedImageSet()
    : _domain(), _function()
{
}

TaylorConstrainedImageSet* TaylorConstrainedImageSet::clone() const
{
    return new TaylorConstrainedImageSet(*this);
}

TaylorConstrainedImageSet::TaylorConstrainedImageSet(Box box)
{
    // Ensure domain elements have nonempty radius
    const double min=std::numeric_limits<double>::min();
    IntervalVector domain=box;
    for(uint i=0; i!=domain.size(); ++i) {
        if(domain[i].radius()==0) {
            domain[i]+=Interval(-min,+min);
        }
    }
    this->_function=VectorTaylorFunction::identity(box);
}


TaylorConstrainedImageSet::TaylorConstrainedImageSet(Box box, VectorFunction function)
{
    ARIADNE_ASSERT_MSG(box.size()==function.argument_size(),"domain="<<box<<", function="<<function);
    const double min=std::numeric_limits<double>::min();
    IntervalVector domain=box;
    for(uint i=0; i!=domain.size(); ++i) {
        if(domain[i].radius()==0) {
            domain[i]+=Interval(-min,+min);
        }
    }

    this->_function=VectorTaylorFunction(domain,function);
}

TaylorConstrainedImageSet::TaylorConstrainedImageSet(Box box, VectorFunction function, List<NonlinearConstraint> constraints)
{
    ARIADNE_ASSERT_MSG(box.size()==function.argument_size(),"domain="<<box<<", function="<<function);
    const double min=std::numeric_limits<double>::min();
    IntervalVector domain=box;
    for(uint i=0; i!=domain.size(); ++i) {
        if(domain[i].radius()==0) {
            domain[i]+=Interval(-min,+min);
        }
    }

    this->_function=VectorTaylorFunction(domain,function);

    for(uint i=0; i!=constraints.size(); ++i) {
        ARIADNE_ASSERT_MSG(box.size()==constraints[i].function().argument_size(),"domain="<<box<<", constraint="<<constraints[i]);
        if(constraints[i].bounds().singleton()) {
            this->new_equality_constraint(constraints[i].function()-constraints[i].bounds().midpoint());
        } else {
            if(constraints[i].bounds().lower()>-inf<Float>()) {
                this->new_negative_constraint(constraints[i].bounds().lower()-constraints[i].function());
            }
            if(constraints[i].bounds().upper()<+inf<Float>()) {
                this->new_negative_constraint(constraints[i].function()-constraints[i].bounds().upper());
            }
        }
    }

}

TaylorConstrainedImageSet::TaylorConstrainedImageSet(Box box, VectorFunction function, NonlinearConstraint constraint)
{
    *this=TaylorConstrainedImageSet(box,function,make_list(constraint));
}

TaylorConstrainedImageSet::TaylorConstrainedImageSet(Box box, const List<ScalarFunction>& components)
{
    const double min=std::numeric_limits<double>::min();
    IntervalVector domain=box;
    for(uint i=0; i!=domain.size(); ++i) {
        if(domain[i].radius()==0) {
            domain[i]+=Interval(-min,+min);
        }
    }

    this->_function=VectorTaylorFunction(components.size(),domain);
    for(uint i=0; i!=components.size(); ++i) {
        this->_function.set(i,ScalarTaylorFunction(domain,components[i]));
    }

    this->_check();
}


TaylorConstrainedImageSet::TaylorConstrainedImageSet(const VectorTaylorFunction& function)
{
    this->_function=function;
}



tribool TaylorConstrainedImageSet::satisfies(ScalarFunction constraint) const
{
    Interval constraint_range=constraint(this->codomain());
    if(constraint_range.lower()>0.0) { return true; }
    else if(constraint_range.upper()<0.0) { return false; }
    else { return indeterminate; }
}


void TaylorConstrainedImageSet::substitute(uint j, ScalarTaylorFunction v)
{
    ARIADNE_ASSERT_MSG(v.argument_size()+1u==this->number_of_parameters(),
                       "number_of_parameters="<<this->number_of_parameters()<<", variable="<<v);
    this->_function = Ariadne::substitute(this->_function,j,v);
    for(List<ScalarTaylorFunction>::iterator iter=this->_constraints.begin(); iter!=this->_constraints.end(); ++iter) {
        *iter = Ariadne::substitute(*iter,j,v);
    }
    for(List<ScalarTaylorFunction>::iterator iter=this->_equations.begin(); iter!=this->_equations.end(); ++iter) {
        *iter = Ariadne::substitute(*iter,j,v);
    }

    this->_check();
}

void TaylorConstrainedImageSet::apply_map(VectorFunction map)
{
    ARIADNE_ASSERT_MSG(map.argument_size()==this->dimension(),"dimension="<<this->dimension()<<", map="<<map);
    VectorTaylorFunction& function=this->_function;
    function=compose(map,function);
    this->_check();
}

void TaylorConstrainedImageSet::apply_map(VectorTaylorFunction map)
{
    ARIADNE_ASSERT_MSG(map.argument_size()==this->dimension(),"dimension="<<this->dimension()<<", map="<<map);
    VectorTaylorFunction& function=this->_function;
    function=compose(map,function);
    this->_check();
}

void TaylorConstrainedImageSet::apply_flow(VectorFunction flow, Interval time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    this->_function=compose(flow,combine(this->_function,VectorTaylorFunction::identity(Vector<Interval>(1u,time))));
    for(List<ScalarTaylorFunction>::iterator iter=this->_constraints.begin(); iter!=this->_constraints.end(); ++iter) {
        *iter=embed(*iter,time);
    }
    for(List<ScalarTaylorFunction>::iterator iter=this->_equations.begin(); iter!=this->_equations.end(); ++iter) {
        *iter=embed(*iter,time);
    }
    this->_check();
}

void TaylorConstrainedImageSet::apply_flow(VectorTaylorFunction flow, Float time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    this->_function=compose(flow,join(this->_function,ScalarTaylorFunction::constant(this->_function.domain(),time)));
    for(List<ScalarTaylorFunction>::iterator iter=this->_constraints.begin(); iter!=this->_constraints.end(); ++iter) {
        *iter=embed(*iter,time);
    }
    for(List<ScalarTaylorFunction>::iterator iter=this->_equations.begin(); iter!=this->_equations.end(); ++iter) {
        *iter=embed(*iter,time);
    }
    this->_check();
}

void TaylorConstrainedImageSet::apply_flow(VectorTaylorFunction flow, Interval time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    this->_function=compose(flow,combine(this->_function,ScalarTaylorFunction::identity(time)));
    for(List<ScalarTaylorFunction>::iterator iter=this->_constraints.begin(); iter!=this->_constraints.end(); ++iter) {
        *iter=embed(*iter,time);
    }
    for(List<ScalarTaylorFunction>::iterator iter=this->_equations.begin(); iter!=this->_equations.end(); ++iter) {
        *iter=embed(*iter,time);
    }
    this->_check();
}


void TaylorConstrainedImageSet::new_negative_constraint(ScalarFunction constraint) {
    ARIADNE_ASSERT_MSG(constraint.argument_size()==this->domain().size(),"domain="<<this->domain()<<", constraint="<<constraint);
    this->_constraints.append(ScalarTaylorFunction(this->domain(),constraint));
}

void TaylorConstrainedImageSet::new_negative_constraint(ScalarTaylorFunction constraint) {
    ARIADNE_ASSERT_MSG(constraint.domain()==this->domain(),"domain="<<this->domain()<<", constraint="<<constraint);
    this->_constraints.append(constraint);
}

void TaylorConstrainedImageSet::new_equality_constraint(ScalarFunction constraint) {
    ARIADNE_ASSERT_MSG(constraint.argument_size()==this->domain().size(),"domain="<<this->domain()<<", constraint="<<constraint);
    this->_equations.append(ScalarTaylorFunction(this->domain(),constraint));
}





IntervalVector TaylorConstrainedImageSet::domain() const {
    return this->_function.domain();
}

IntervalVector TaylorConstrainedImageSet::codomain() const {
    return Box(this->_function.range()).bounding_box();
}

VectorFunction TaylorConstrainedImageSet::function() const {
    return VectorFunction(Vector< Polynomial<Real> >(this->_function.polynomial()));
}

VectorTaylorFunction const& TaylorConstrainedImageSet::taylor_function() const {
    return this->_function;
}

uint TaylorConstrainedImageSet::dimension() const {
    return this->_function.result_size();
}

uint TaylorConstrainedImageSet::number_of_parameters() const {
    return this->_function.argument_size();
}

Box TaylorConstrainedImageSet::bounding_box() const {
    return Box(this->_function.codomain()).bounding_box();
}

Float TaylorConstrainedImageSet::radius() const {
    return this->bounding_box().radius();
}

Point TaylorConstrainedImageSet::centre() const {
    return this->bounding_box().centre();
}


tribool TaylorConstrainedImageSet::bounded() const {
    return Box(this->domain()).bounded() || indeterminate;
}

tribool TaylorConstrainedImageSet::empty() const {
    return this->disjoint(this->bounding_box());
}

tribool TaylorConstrainedImageSet::disjoint(Box bx) const {
    static bool warn=true;
    if(warn) {
        ARIADNE_WARN("TaylorConstrainedImageSet::disjoint(Box) is not correctly implemented");
        warn=false;
    }
    return this->affine_approximation().disjoint(bx);
}


Pair<TaylorConstrainedImageSet,TaylorConstrainedImageSet>
TaylorConstrainedImageSet::split(uint d) const
{
    Vector<Interval> subdomain1,subdomain2;
    make_lpair(subdomain1,subdomain2)=Ariadne::split(this->_function.domain(),d);

    VectorTaylorFunction function1,function2;
    make_lpair(function1,function2)=Ariadne::split(this->_function,d);

    Pair<TaylorConstrainedImageSet,TaylorConstrainedImageSet>
    result=make_pair(TaylorConstrainedImageSet(function1),TaylorConstrainedImageSet(function2));
    TaylorConstrainedImageSet& result1=result.first;
    TaylorConstrainedImageSet& result2=result.second;

    ScalarTaylorFunction constraint1,constraint2;
    for(List<ScalarTaylorFunction>::const_iterator iter=this->_constraints.begin();
        iter!=this->_constraints.end(); ++iter)
    {
        const ScalarTaylorFunction& constraint=*iter;
        make_lpair(constraint1,constraint2)=Ariadne::split(constraint,d);
        result1._constraints.append(constraint1);
        result2._constraints.append(constraint2);
    }

    ScalarTaylorFunction equation1,equation2;
    for(List<ScalarTaylorFunction>::const_iterator iter=this->_equations.begin();
        iter!=this->_equations.end(); ++iter)
    {
        const ScalarTaylorFunction& equation=*iter;
        make_lpair(equation1,equation1)=Ariadne::split(equation,d);
        result1._equations.append(equation1);
        result2._equations.append(equation1);
    }

    return make_pair(result1,result2);
}






ScalarFunction make_function(const ScalarTaylorFunction& stf) {
    return ScalarFunction(stf.polynomial())+Real(Interval(-stf.error(),+stf.error()));
}

void adjoin_outer_approximation_to(GridTreeSet& r, const Box& d, const VectorFunction& fg, const Box& c, const GridCell& b, Point& x, Point& y, int e)
{
    // When making a new starting primal point, need to move components away from zero
    // This constant shows how far away from zero the points are
    static const double XSIGMA=0.125;
    static const double TERR=-1.0/((1<<e)*1024.0);

    uint verbosity=0u;

    const uint m=fg.argument_size();
    const uint n=fg.result_size();
    ARIADNE_LOG(2,"\nadjoin_outer_approximation(...)\n");
    ARIADNE_LOG(2,"  dom="<<d<<" cnst="<<c<<" cell="<<b.box()<<" dpth="<<b.tree_depth()<<" e="<<e<<"\n");

    ConstraintSolver solver;
    NonlinearInteriorPointOptimiser optimiser;

    Float t;
    Point z(x.size());

    if(subset(b,r)) {
        return;
    }

    Box bx=join(static_cast<const IntervalVector&>(b.box()),static_cast<const IntervalVector&>(c));

    optimiser.compute_tz(d,fg,bx,y,t,z);
    for(uint i=0; i!=12; ++i) {
        ARIADNE_LOG(4," t="<<t);
        optimiser.linearised_feasibility_step(d,fg,bx,x,y,z,t);
        if(t>0) { break; }
    }
    ARIADNE_LOG(4,"\n  t="<<t<<"\n  y="<<y<<"\n    x="<<x<<"\n    z="<<z<<"\n");

    if(t<TERR) {
        // Probably disjoint, so try to prove this
        Box nd=d;

        // Use the computed dual variables to try to make a scalar function which is negative over the entire domain.
        // This should be easier than using all constraints separately
        ScalarFunction xg=ScalarFunction::constant(m,0);
        Interval cnst=0.0;
        for(uint j=0; j!=n; ++j) {
            xg = xg - (Real(x[j])-x[n+j])*fg[j];
            cnst += (bx[j].upper()*x[j]-bx[j].lower()*x[n+j]);
        }
        for(uint i=0; i!=m; ++i) {
            xg = xg - (x[2*n+i]-x[2*n+m+i])*ScalarFunction::coordinate(m,i);
            cnst += (d[i].upper()*x[2*n+i]-d[i].lower()*x[2*n+m+i]);
        }
        xg = Real(cnst) + xg;

        ARIADNE_LOG(4,"    xg="<<xg<<"\n");
        ScalarTaylorFunction txg(d,xg);
        ARIADNE_LOG(4,"    txg="<<txg.polynomial()<<"\n");

        xg=ScalarFunction(txg.polynomial());
        NonlinearConstraint constraint=(xg>=0.0);

        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        solver.hull_reduce(constraint,nd);
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        if(nd.empty()) {
            ARIADNE_LOG(4,"  Proved disjointness using hull reduce\n");
            return;
        }

        for(uint i=0; i!=m; ++i) {
            solver.box_reduce(constraint,nd,i);
            ARIADNE_LOG(8,"  dom="<<nd<<"\n");
            if(nd.empty()) { ARIADNE_LOG(4,"  Proved disjointness using box reduce\n"); return; }
        }
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");

        //Pair<Box,Box> sd=solver.split(List<NonlinearConstraint>(1u,constraint),d);
        ARIADNE_LOG(4,"  Splitting domain\n");
        Pair<Box,Box> sd=d.split();
        Point nx = (1.0-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());
        Point ny = midpoint(sd.first);
        adjoin_outer_approximation_to(r, sd.first, fg, c, b, nx, ny, e);
        nx = (1.0-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());
        ny = midpoint(sd.second);
        adjoin_outer_approximation_to(r, sd.second, fg, c, b, x, ny, e);
    }

    if(b.tree_depth()>=e*int(b.dimension())) {
        ARIADNE_LOG(4,"  Adjoining cell "<<b.box()<<"\n");
        r.adjoin(b);
    } else {
        ARIADNE_LOG(4,"  Splitting cell; t="<<t<<"\n");
        Pair<GridCell,GridCell> sb = b.split();
        Point sx = (1-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());
        Point sy = y;
        adjoin_outer_approximation_to(r,d,fg,c,sb.first,sx,sy,e);
        sx = (1-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());
        sy = y;
        adjoin_outer_approximation_to(r,d,fg,c,sb.second,sx,sy,e);
    }


}


void adjoin_outer_approximation_to(GridTreeSet&, const Box& domain, const VectorFunction& function, const VectorFunction& negative_constraints, int depth);

void TaylorConstrainedImageSet::adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const
{
    this->affine_adjoin_outer_approximation_to(paving,depth);
    //this->this->_adjoin_outer_approximation_to(paving,this->bounding_box(),depth);
}

void TaylorConstrainedImageSet::subdivision_adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const
{
    Vector<Float> errors(paving.dimension());
    for(uint i=0; i!=errors.size(); ++i) {
        errors[i]=paving.grid().lengths()[i]/(1<<depth);
    }
    this->_subdivision_adjoin_outer_approximation_to(paving,this->domain(),depth,errors);
}

void TaylorConstrainedImageSet::affine_adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const
{
    this->affine_approximation().adjoin_outer_approximation_to(paving,depth);
}

void TaylorConstrainedImageSet::constraint_adjoin_outer_approximation_to(GridTreeSet& p, int e) const
{
    ARIADNE_ASSERT(p.dimension()==this->dimension());
    const Box& d=this->domain();
    const VectorFunction& f=this->function();

    VectorFunction g(this->_constraints.size(),d.size());
    uint i=0;
    for(List<ScalarTaylorFunction>::const_iterator citer=this->_constraints.begin(); citer!=this->_constraints.end(); ++citer) {
        g.set(i,make_function(*citer));
        ++i;
    }

    GridCell b=GridCell::smallest_enclosing_primary_cell(g(d),p.grid());
    Box c=intersection(g(d)+IntervalVector(g.result_size(),Interval(-1,1)),Box(g.result_size(),Interval(-inf<Float>(),0.0)));

    Point y=midpoint(d);
    const uint l=(d.size()+f.result_size()+g.result_size())*2;
    Point x(l); for(uint k=0; k!=l; ++k) { x[k]=1.0/l; }

    VectorFunction fg=join(f,g);

    Ariadne::adjoin_outer_approximation_to(p,d,fg,c,b,x,y,e);

}

GridTreeSet TaylorConstrainedImageSet::outer_approximation(const Grid& grid, int depth) const
{
    GridTreeSet paving(grid);
    this->adjoin_outer_approximation_to(paving,depth);
    return paving;
}

GridTreeSet TaylorConstrainedImageSet::subdivision_outer_approximation(const Grid& grid, int depth) const
{
    GridTreeSet paving(grid);
    this->subdivision_adjoin_outer_approximation_to(paving,depth);
    return paving;
}

GridTreeSet TaylorConstrainedImageSet::affine_outer_approximation(const Grid& grid, int depth) const
{
    GridTreeSet paving(grid);
    this->affine_adjoin_outer_approximation_to(paving,depth);
    return paving;
}



TaylorConstrainedImageSet TaylorConstrainedImageSet::restriction(const Vector<Interval>& subdomain) const
{
    ARIADNE_ASSERT_MSG(subdomain.size()==this->number_of_parameters(),"set="<<*this<<", subdomain="<<subdomain);
    ARIADNE_ASSERT_MSG(subset(subdomain,this->domain()),"set.domain()="<<this->domain()<<", subdomain="<<subdomain);
    TaylorConstrainedImageSet result(*this);
    result._function=Ariadne::restrict(result._function,subdomain);
    ScalarTaylorFunction new_constraint;
    for(List<ScalarTaylorFunction>::iterator iter=result._constraints.begin();
        iter!=result._constraints.end(); ++iter)
    {
        ScalarTaylorFunction& constraint=*iter;
        constraint=Ariadne::restrict(constraint,subdomain);
    }
    for(List<ScalarTaylorFunction>::iterator iter=result._equations.begin();
        iter!=result._equations.end(); ++iter)
    {
        ScalarTaylorFunction& equation=*iter;
        equation=Ariadne::restrict(equation,subdomain);
    }
    return result;
}



void TaylorConstrainedImageSet::draw(CanvasInterface& canvas) const {
    this->affine_draw(canvas,0u);
}

void TaylorConstrainedImageSet::box_draw(CanvasInterface& canvas) const {
    this->bounding_box().draw(canvas);
}

void TaylorConstrainedImageSet::affine_draw(CanvasInterface& canvas, uint accuracy) const {
    static const int depth=accuracy;
    List<Box> subdomains;
    List<Box> splitdomains;
    subdomains.append(this->domain());
    Box splitdomain1,splitdomain2;
    for(int i=0; i!=depth; ++i) {
        for(uint k=0; k!=this->number_of_parameters(); ++k) {
            for(uint n=0; n!=subdomains.size(); ++n) {
                make_lpair(splitdomain1,splitdomain2)=subdomains[n].split(k);
                splitdomains.append(splitdomain1);
                splitdomains.append(splitdomain2);
            }
            subdomains.swap(splitdomains);
            splitdomains.clear();
        }
    }

    for(uint n=0; n!=subdomains.size(); ++n) {
        this->restriction(subdomains[n]).affine_over_approximation().draw(canvas);
    }
};


Map<List<DiscreteEvent>,ScalarFunction> pretty(const Map<List<DiscreteEvent>,ScalarTaylorFunction>& constraints) {
    Map<List<DiscreteEvent>,ScalarFunction> result;
    for(Map<List<DiscreteEvent>,ScalarTaylorFunction>::const_iterator iter=constraints.begin();
    iter!=constraints.end(); ++iter) {
        result.insert(iter->first,iter->second.function());
    }
    return result;
}


template<class K, class V> Map<K,V> filter(const Map<K,V>& m, const Set<K>& s) {
    Map<K,V> r;
    for(typename Set<K>::const_iterator iter=s.begin(); iter!=s.end(); ++iter) {
        r.insert(*m.find(*iter));
    }
    return r;
}

std::ostream& TaylorConstrainedImageSet::write(std::ostream& os) const {
    os << "TaylorConstrainedImageSet";
    os << "(\n  domain=" << this->domain();
    os << ",\n  range=" << this->bounding_box();
    os << ",\n  function=" << this->taylor_function();
    os << ",\n  constraints=" << this->_constraints;
    os << ",\n  equations=" << this->_equations;
    os << "\n)\n";
    return os;
}



void TaylorConstrainedImageSet::
_subdivision_adjoin_outer_approximation_to(GridTreeSet& gts, const IntervalVector& subdomain,
                                           uint depth, const FloatVector& errors) const
{
    // How small an over-approximating box needs to be relative to the cell size
    static const double RELATIVE_SMALLNESS=0.5;

    typedef List<ScalarTaylorFunction>::const_iterator const_iterator;

    for(const_iterator iter=this->_constraints.begin(); iter!=this->_constraints.end(); ++iter) {
        const ScalarTaylorFunction& constraint=*iter;;
        Interval constraint_range=constraint.evaluate(subdomain);
        if(constraint_range.lower() > 0.0) {
            return;
        }
    }
    for(const_iterator iter=this->_equations.begin(); iter!=this->_equations.end(); ++iter) {
        const ScalarTaylorFunction& constraint=*iter;
        Interval constraint_range=constraint.evaluate(subdomain);
        if(constraint_range.lower() > 0.0 || constraint_range.upper() < 0.0 ) {
            return;
        }
    }

    Box range=evaluate(this->_function,subdomain);
    bool small=true;
    for(uint i=0; i!=range.size(); ++i) {
        if(range[i].radius()>errors[i]*RELATIVE_SMALLNESS) {
            small=false;
            break;
        }
    }

    if(small) {
        gts.adjoin_outer_approximation(range,depth);
    } else {
        Vector<Interval> subdomain1,subdomain2;
        make_lpair(subdomain1,subdomain2)=Ariadne::split(subdomain);
        this->_subdivision_adjoin_outer_approximation_to(gts,subdomain1,depth,errors);
        this->_subdivision_adjoin_outer_approximation_to(gts,subdomain2,depth,errors);
    }
}


AffineSet
TaylorConstrainedImageSet::affine_approximation() const
{
    this->_check();
    typedef List<ScalarTaylorFunction>::const_iterator const_iterator;

    const uint nx=this->dimension();
    const uint np=this->number_of_parameters();

    TaylorConstrainedImageSet set(*this);

    if(set._equations.size()>0) {
        set._solve_zero_constraints();
    }
    this->_check();

    Vector<Float> h(nx);
    Matrix<Float> G(nx,np);
    for(uint i=0; i!=nx; ++i) {
        ScalarTaylorFunction component=set._function[i];
        h[i]=component.model().value();
        for(uint j=0; j!=np; ++j) {
            G[i][j]=component.model().gradient(j);
        }
    }
    AffineSet result(G,h);

    Vector<Float> a(np);
    Float b;

    for(const_iterator iter=set._constraints.begin();
            iter!=set._constraints.end(); ++iter) {
        const ScalarTaylorFunction& constraint=*iter;
        b=-constraint.model().value();
        for(uint j=0; j!=np; ++j) { a[j]=constraint.model().gradient(j); }
        result.new_inequality_constraint(a,b);
    }

    for(const_iterator iter=set._equations.begin();
            iter!=set._equations.end(); ++iter) {
        const ScalarTaylorFunction& constraint=*iter;
        b=-constraint.model().value();
        for(uint j=0; j!=np; ++j) { a[j]=constraint.model().gradient(j); }
        result.new_equality_constraint(a,b);
    }
    return result;
}

AffineSet
TaylorConstrainedImageSet::affine_over_approximation() const
{
    this->_check();
    typedef List<ScalarTaylorFunction>::const_iterator const_iterator;

    const uint nx=this->dimension();
    const uint np=this->number_of_parameters();

    TaylorConstrainedImageSet set(*this);

    if(set._equations.size()>0) {
        set._solve_zero_constraints();
    }
    this->_check();

    for(uint i=0; i!=nx; ++i) {
        const_cast<TaylorModel&>(set._function.models()[i]).truncate(1u);
    }

    Vector<Float> h(nx);
    Matrix<Float> G(nx,np+nx);
    for(uint i=0; i!=nx; ++i) {
        ScalarTaylorFunction component=set._function[i];
        h[i]=component.model().value();
        for(uint j=0; j!=np; ++j) {
            G[i][j]=component.model().gradient(j);
        }
        G[i][np+i]=component.model().error();
    }
    AffineSet result(G,h);

    Vector<Float> a(np+nx);
    Float b;

    for(const_iterator iter=set._constraints.begin();
            iter!=set._constraints.end(); ++iter) {
        const ScalarTaylorFunction& constraint=*iter;
        b=-constraint.model().value();
        for(uint j=0; j!=np; ++j) { a[j]=constraint.model().gradient(j); }
        result.new_inequality_constraint(a,b);
    }

    for(const_iterator iter=set._equations.begin();
            iter!=set._equations.end(); ++iter) {
        const ScalarTaylorFunction& constraint=*iter;
        b=-constraint.model().value();
        for(uint j=0; j!=np; ++j) { a[j]=constraint.model().gradient(j); }
        result.new_equality_constraint(a,b);
    }
    return result;
}

TaylorConstrainedImageSet product(const TaylorConstrainedImageSet& set, const Interval& ivl) {
    typedef List<ScalarTaylorFunction>::const_iterator const_iterator;

    VectorTaylorFunction new_function=combine(set.taylor_function(),ScalarTaylorFunction::identity(ivl));

    TaylorConstrainedImageSet result(new_function);
    for(const_iterator iter=set._constraints.begin(); iter!=set._constraints.end(); ++iter) {
        result._constraints.append(embed(*iter,ivl));
    }
    for(const_iterator iter=set._equations.begin(); iter!=set._equations.end(); ++iter) {
        result._equations.append(embed(*iter,ivl));
    }

    return result;
}

} // namespace Ariadne


#else // Not HAVE_CAIRO_H

namespace Ariadne {

void plot(const char* filename, const Box& bbx, const TaylorImageSet& set)
{
    throw std::runtime_error("No facilities for drawing graphics are available.");
}

}

#endif // HAVE_CAIRO_H
