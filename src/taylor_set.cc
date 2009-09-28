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
#include "taylor_set.h"

#include "list_set.h"
#include "grid_set.h"

#include "zonotope.h"
#include "polytope.h"
#include "graphics_interface.h"

namespace Ariadne {

TaylorSet::TaylorSet(uint d, uint ng)
    : _models(d,ng)
{
}


TaylorSet::TaylorSet(const VectorFunction& f, const Vector<Interval>& d)
    : _models(f.result_size())
{
    Vector<TaylorModel> x=TaylorModel::scalings(d);
    this->_models=f.evaluate(x);
}

TaylorSet::TaylorSet(const Vector<Interval>& bx)
    : _models(TaylorModel::scalings(bx))
{
}


TaylorSet::TaylorSet(const Vector<TaylorModel>& tvs)
    : _models(tvs)
{
    if(this->_models.size()>0) {
        Vector<Interval> unit_domain(this->_models[0].argument_size(),Interval(-1,+1));
        for(size_t i=0; i!=this->_models.size(); ++i) {
            ARIADNE_ASSERT(this->_models[i].domain()==unit_domain);
        }
    }
}

TaylorSet::TaylorSet(const Vector< Expansion<Float> >& f, const Vector<Float>& e)
    : _models(f.size())
{
    ARIADNE_ASSERT(f.size()==e.size());
    if(f.size()>0) {
        for(uint i=0; i!=f.size(); ++i) { (*this)[i]=TaylorModel(f[i],e[i]); }
    }
}

TaylorSet::TaylorSet(uint rs, uint as, uint deg, double x0, ...)
    : _models(rs)
{
    double x=x0;
    va_list args;
    va_start(args,x0);
    for(uint i=0; i!=rs; ++i) {
        (*this)[i]=TaylorModel(as);
        for(MultiIndex j(as); j.degree()<=deg; ++j) {
            //std::cerr<<"  "<<int(j.degree())<<":"<<j<<" "<<x<<"\n";
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
operator==(const TaylorSet& ts1,const TaylorSet& ts2)
{
    return ts1.models()==ts2.models();
}


void
TaylorSet::set_accuracy(shared_ptr<TaylorModel::Accuracy> acc_ptr)
{
    for(uint i=0; i!=this->dimension(); ++i) {
        this->_models[i].set_accuracy(acc_ptr);
    }
}

shared_ptr<TaylorModel::Accuracy>
TaylorSet::accuracy_ptr() const
{
    return this->_models[0].accuracy_ptr();
}


Vector<Float>
TaylorSet::centre() const
{
    Vector<Float> result(this->dimension());
    for(uint i=0; i!=this->dimension(); ++i) {
        result[i]=this->_models[i].value();
    }
    return result;
}



Float
TaylorSet::radius() const
{
    Float result=0.0;
    for(uint i=0; i!=this->dimension(); ++i) {
        result=max(result,this->models()[i].range().radius());
    }
    return result;
}


tribool
TaylorSet::disjoint(const Box& bx) const
{
    ARIADNE_ASSERT_MSG(this->dimension()==bx.dimension(),"Dimension of "<<*this<<" is different from dimension of "<<bx);

    // Expand the bounding box by the error term
    Box ebx(bx);
    for(uint i=0; i!=bx.dimension(); ++i) {
        ebx[i]+=Interval(-this->_models[i].error(),+this->_models[i].error());
    }

    // Copy the box eliminating error terms
    TaylorSet efts(*this);
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
        std::pair<TaylorSet,TaylorSet> split_efts=this->split();
        return split_efts.first.disjoint(ebx) && split_efts.second.disjoint(ebx);
    }
}



tribool
TaylorSet::overlaps(const Box& bx) const
{
    return !this->disjoint(bx);
}


tribool
TaylorSet::inside(const Box& bx) const
{
    Vector<Interval> bb=this->bounding_box();
    return Ariadne::inside(bb,bx) || indeterminate;
}


Box
TaylorSet::bounding_box() const
{
    Box r(this->dimension());
    for(uint i=0; i!=this->dimension(); ++i) {
        r[i]=(*this)[i].range();
    }
    return r;
}


TaylorSet
TaylorSet::linearise() const
{
    TaylorSet res(*this);
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
TaylorSet::discretise(const Float& eps) const
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
TaylorSet::discretise(const Grid& g, uint d) const
{
    GridTreeSet gts(g);
    return this->discretise(gts,d);
}

GridTreeSet&
TaylorSet::discretise(GridTreeSet& gts, uint d) const
{
    adjoin_outer_approximation(gts,*this,d);
    gts.recombine();
    return gts;
}


TaylorModel apply(const ScalarFunction& f, const TaylorSet& ts)
{
    return f.evaluate(ts.models());
}

TaylorSet apply(const VectorFunction& f, const TaylorSet& ts)
{
    return f.evaluate(ts.models());
}

TaylorModel apply(const ScalarTaylorFunction& tf, const TaylorSet& ts)
{
    ARIADNE_ASSERT_MSG(subset(ts.range(),tf.domain()),"tf="<<tf<<" ts="<<ts);
    return unchecked_compose(tf.model(),unscale(ts.models(),tf.domain()));
}

TaylorSet apply(const VectorTaylorFunction& tf, const TaylorSet& ts)
{
    ARIADNE_ASSERT_MSG(possibly(subset(ts.range(),tf.domain())),
        std::setprecision(18)<<"\n  tf="<<tf<<"\n  ts="<<ts<<"\n  ts.range() ="<<ts.range()<<"\n  tf.domain()="<<tf.domain());
    return unchecked_compose(tf.models(),unscale(ts.models(),tf.domain()));
}

TaylorModel unchecked_apply(const ScalarTaylorFunction& tf, const TaylorSet& ts)
{
    return unchecked_compose(tf.model(),unscale(ts.models(),tf.domain()));
}

TaylorSet unchecked_apply(const VectorTaylorFunction& tf, const TaylorSet& ts)
{
    return unchecked_compose(tf.models(),unscale(ts.models(),tf.domain()));
}


Zonotope
zonotope(const TaylorSet& ts)
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




pair<TaylorSet,TaylorSet>
TaylorSet::split(uint j) const
{
    pair< Vector<TaylorModel>, Vector<TaylorModel> > s=Ariadne::split(this->_models,j);
    return make_pair( TaylorSet(s.first),
                      TaylorSet(s.second) );
}

pair<TaylorSet,TaylorSet>
TaylorSet::split() const
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
        pair<TaylorSet,TaylorSet> result=make_pair(*this,*this);
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
    //std::cerr<<"range="<<range<<"\n";
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
adjoin_outer_approximation1(GridTreeSet& grid_set, const TaylorSet& set, uint depth)
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
                             const TaylorSet& set, const Vector<Interval>& errors,
                             Float eps, uint depth)
{
    Box range=set.range();
    //std::cerr<<"range="<<range<<"\n";
    if(range.radius()<eps) {
        grid_set.adjoin_over_approximation(range+errors,depth);
    } else {
        std::pair<TaylorSet,TaylorSet> subsets=set.split();
        _adjoin_outer_approximation2(grid_set,subsets.first,errors,eps,depth);
        _adjoin_outer_approximation2(grid_set,subsets.second,errors,eps,depth);
    }
}


void
adjoin_outer_approximation2(GridTreeSet& grid_set, const TaylorSet& set, uint depth)
{
    Box domain(set.generators_size(),Interval(-1,+1));
    Float eps=1.0/(1<<(depth/set.dimension()));
    TaylorSet error_free_set(set);
    Vector<Interval> errors(set.dimension());
    for(uint i=0; i!=set.dimension(); ++i) {
        errors[i]=set.models()[i].error()*Interval(-1,+1);
        const_cast<TaylorModel&>(error_free_set.models()[i]).set_error(0.0);
    }
    _adjoin_outer_approximation2(grid_set,error_free_set,errors,eps,depth);
    grid_set.recombine();
}

void
_adjoin_outer_approximation3(GridTreeSet& grid_set, const TaylorSet& set, Box& domain,
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
adjoin_outer_approximation3(GridTreeSet& grid_set, const TaylorSet& set, uint depth)
{
    TaylorSet error_free_set=set.subsume();
    Box domain=error_free_set.domain();
    GridOpenCell cell=GridOpenCell::outer_approximation(error_free_set.bounding_box(),grid_set.grid());
    _adjoin_outer_approximation3(grid_set,error_free_set,domain,cell,depth);
    grid_set.recombine();
}


void
_adjoin_outer_approximation4(GridTreeSet& grid_set, const TaylorSet& set,
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

    // Subdivide the TaylorSet and try again
    std::pair<TaylorSet,TaylorSet> subsets=set.split();
    _adjoin_outer_approximation4(grid_set,subsets.first,cell,depth);
    _adjoin_outer_approximation4(grid_set,subsets.second,cell,depth);
}

void
adjoin_outer_approximation4(GridTreeSet& grid_set, const TaylorSet& set, uint depth)
{
    GridOpenCell cell=GridOpenCell::outer_approximation(set.bounding_box(),grid_set.grid());
    TaylorSet error_free_set=set.subsume();
    _adjoin_outer_approximation4(grid_set,error_free_set,cell,depth);
    grid_set.recombine();
}

// Profiling suggests that adjoin_outer_approximation1 is most efficient.
// adjoin_outer_approximation3 (using subdivision of the set) is better than
// adjoin_outer_approximation2 (using subdivision of the domain).
void
adjoin_outer_approximation(GridTreeSet& grid_set, const TaylorSet& set, uint depth)
{
    adjoin_outer_approximation1(grid_set,set,depth);
}


GridTreeSet
discretise1(const TaylorSet& ts, const Grid& g, uint d)
{
    GridTreeSet gts(g);
    adjoin_outer_approximation1(gts,ts,d);
    gts.recombine();
    return gts;
}

GridTreeSet
discretise2(const TaylorSet& ts, const Grid& g, uint d)
{
    GridTreeSet gts(g);
    adjoin_outer_approximation2(gts,ts,d);
    gts.recombine();
    return gts;
}

GridTreeSet
discretise3(const TaylorSet& ts, const Grid& g, uint d)
{
    GridTreeSet gts(g);
    adjoin_outer_approximation3(gts,ts,d);
    gts.recombine();
    return gts;
}

GridTreeSet
discretise4(const TaylorSet& ts, const Grid& g, uint d)
{
    GridTreeSet gts(g);
    adjoin_outer_approximation4(gts,ts,d);
    gts.recombine();
    return gts;
}


Matrix<Float>
TaylorSet::jacobian() const
{
    Matrix<Float> J(this->dimension(),this->generators_size());
    for(uint i=0; i!=this->dimension(); ++i) {
        for(uint j=0; j!=this->generators_size(); ++j) {
            J[i][j]=this->models()[i][MultiIndex::unit(this->generators_size(),j)];
        }
    }
    return J;
}

TaylorSet
TaylorSet::subsume() const
{
    return this->subsume(0.0);
}


TaylorSet
TaylorSet::subsume(double eps) const
{
    uint d=this->dimension();
    uint ng=this->generators_size();

    // Compute the number of terms whose error is greater the eps
    uint ne=0;
    for(uint i=0; i!=d; ++i) {
        if(this->models()[i].error()>eps) { ++ne; }
    }

    TaylorSet result(embed(this->models(),ne));
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
outer_approximation(const TaylorSet& set, const Grid& grid, uint depth)
{
    ARIADNE_ASSERT(set.dimension()==grid.dimension());
    GridTreeSet grid_set(grid);
    adjoin_outer_approximation(grid_set,set,depth);
    return grid_set;
}



TaylorSet
TaylorSet::recondition() const
{
    Matrix<Float> T=triangular_multiplier(this->jacobian());
    Vector<TaylorModel> scal(T.row_size(),T.column_size());
    for(uint i=0; i!=T.row_size(); ++i) {
        for(uint j=0; j!=T.column_size(); ++j) {
            scal[i][MultiIndex::unit(this->generators_size(),j)]=T[i][j];
        }
    }
    return TaylorSet(compose(this->models(),scal));
}

std::ostream&
TaylorSet::write(std::ostream& os) const
{
    return os << "TaylorSet(" << this->_models << ")";
    os << "TaylorSet(\n";
    os << "  dimension=" << this->dimension() << ",\n" << std::flush;
    os << "  models=" << this->_models << ",\n" << std::flush;
    os << ")\n";
    return os;
}

Vector<Interval> evaluate(const TaylorSet& ts, const Vector<Interval>& iv) {
    Vector<Interval> r(ts.dimension());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=ts[i].evaluate(iv);
    }
    return r;
}

Vector<Interval> evaluate(const TaylorSet& ts, const Vector<Float>& v) {
    Vector<Interval> r(ts.dimension());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=ts[i].evaluate(v);
    }
    return r;
}

Vector<Interval> error(const TaylorSet& ts) {
    Vector<Interval> r(ts.dimension());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=ts[i].error();
    }
    return r;
}


void box_draw(CanvasInterface& fig, const TaylorSet& ts) {
    ts.bounding_box().draw(fig);
}

void boxes_draw(CanvasInterface& fig, const TaylorSet& ts) {
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

void zonotopes_draw_subdivide(ListSet<Zonotope>& zonotopes, const TaylorSet& set, const double resolution)
{
    if(set.bounding_box().radius()>=resolution) {
        std::pair<TaylorSet,TaylorSet> subdivisions=set.split();
        zonotopes_draw_subdivide(zonotopes,subdivisions.first,resolution);
        zonotopes_draw_subdivide(zonotopes,subdivisions.second,resolution);
    } else {
        zonotopes.adjoin(zonotope(set));
    }
}

void zonotopes_draw(CanvasInterface& fig, const TaylorSet& ts) {
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

void affine_draw(CanvasInterface& fig, const TaylorSet& ts) {
    zonotope(ts).draw(fig);
}

void curve_draw(CanvasInterface& fig, const TaylorSet& ts) {
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


void grid_draw(CanvasInterface& fig, const TaylorSet& ts)
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

void standard_draw(CanvasInterface& fig, const TaylorSet& ts) {
    affine_draw(fig,ts);
    //box_draw(fig,ts);
}

void TaylorSet::draw(CanvasInterface& fig) const {
    Ariadne::standard_draw(fig,*this);
}

} // namespace Ariadne


#include "config.h"

#ifdef HAVE_CAIRO_H

#include <cairo/cairo.h>

namespace Ariadne {

void plot(const char* filename, const Box& bbx, const TaylorSet& set)
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

} // namespace Ariadne


#else // Not HAVE_CAIRO_H

namespace Ariadne {

void plot(const char* filename, const Box& bbx, const TaylorSet& set)
{
    throw std::runtime_error("No facilities for drawing graphics are available.");
}

}

#endif // HAVE_CAIRO_H
