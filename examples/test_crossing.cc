/***************************************************************************
 *           test_crossing.h
 *
 *  Copyright  2011  Luca Geretti
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


#include "ariadne.h"
#include "taylor_calculus.h"

using namespace Ariadne;

Matrix<Interval>
jacobian2_range(const Vector<TaylorModel>& f)
{
    uint rs=f.size();
    uint fas=f[0].argument_size();
    uint has=fas-rs;
    Matrix<Interval> J(rs,rs);
    for(uint i=0; i!=rs; ++i) {
        for(TaylorModel::const_iterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            for(uint k=0; k!=rs; ++k) {
                const uint c=iter->key()[has+k];
                if(c>0) {
                    const double& x=iter->data();
                    if(iter->key().degree()==1) { J[i][k]+=x; }
                    else { J[i][k]+=Interval(-1,1)*x*c; }
                    //std::cerr<<"  J="<<J<<" i="<<i<<" a="<<iter->key()<<" k="<<k<<" c="<<c<<" x="<<x<<std::endl;
                }
            }
        }
    }
    return J;
}

int main(int argc, char* argv[])
{
    // System variables
    RealVariable x("x");    
    RealVariable y("y");   
    List<RealVariable> varlist;
    varlist.append(x);
    varlist.append(y);

    RealExpression xdot = 1.0;
    RealExpression ydot = 0.0;

    // Dynamics at the different modes
    List<RealExpression> exprlist;
    exprlist.append(xdot);
    exprlist.append(ydot);
    VectorFunction dyn(exprlist, varlist);

    // Guards are true when f(x) >= 0
    RealExpression guard_expr = x+y*y;
    ScalarFunction guard_sf(guard_expr, varlist);
    VectorFunction guard(1u,guard_sf);

    Box graphics_bx(2,-2.0,2.0,-2.0,2.0);
    int depth = 4;

    Grid grid(2);
    GridTreeSet gts(grid);
    gts.adjoin_outer_approximation(graphics_bx,depth);
    gts.mince(depth);

    TaylorCalculus tc(4,6,1e-10);

		// Sets up the figure
	Figure crossed_fig;
	Figure positively_fig;
	array<uint> xy(2,0,1);
	crossed_fig.set_projection_map(ProjectionFunction(xy,2));
	crossed_fig.set_bounding_box(graphics_bx);
	positively_fig.set_projection_map(ProjectionFunction(xy,2));
	positively_fig.set_bounding_box(graphics_bx);

    for (GridTreeSet::const_iterator cell_it=gts.begin();cell_it!=gts.end();++cell_it) {
        Box bx = cell_it->box();
        tribool crossed = tc.active(guard,bx);

        ScalarFunction derivative=lie_derivative(guard_sf,dyn);
        Interval derivative_range = derivative.evaluate(bx);

        tribool positively_crossing;
        if (derivative_range.lower() > 0)
        	positively_crossing = true;
        else if (derivative_range.upper() < 0)
        	positively_crossing = false;
        else positively_crossing = indeterminate;


		if (definitely(crossed)) {
			crossed_fig.set_fill_colour(Colour(0.0,0.83,0.0));
		}
		else if (indeterminate(crossed)) {
			crossed_fig.set_fill_colour(Colour(1.0,1.0,0.0));
		}
		else {
			crossed_fig.set_fill_colour(Colour(1.0,0.34,0.34));
		}
		if (definitely(positively_crossing)) {
			positively_fig.set_fill_colour(Colour(0.0,0.83,0.0));
		}
		else if (indeterminate(positively_crossing)) {
			positively_fig.set_fill_colour(Colour(1.0,1.0,0.0));
		}
		else {
			positively_fig.set_fill_colour(Colour(1.0,0.34,0.34));
		}

   		crossed_fig.draw(bx);
   		positively_fig.draw(bx);
    }

    crossed_fig.write("test_crossing-crossed");
    positively_fig.write("test_crossing-positively");
}
