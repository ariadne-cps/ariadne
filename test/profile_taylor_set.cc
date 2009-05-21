/***************************************************************************
 *            profile_taylor_set.cc
 *
 *  Copyright 2009  Pieter Collins
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

#include <cassert>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include <string>
using std::string;

#include <boost/timer.hpp>
#include <boost/progress.hpp>

#include "numeric.h"
#include "taylor_set.h"
#include "taylor_expression.h"
#include "taylor_function.h"
#include "list_set.h"
#include "taylor_set.h"
#include "grid_set.h"
#include "box.h"
#include "graphics.h"

namespace Ariadne {


struct ProfileDiscretise {
    typedef GridTreeSet (* FunctionPtr)(const TaylorSet&, const Grid&,uint);
    ProfileDiscretise(FunctionPtr f_, TaylorSet& s_, Grid& g_, uint d_)  : f(f_), s(s_), g(g_), d(d_) { }
    FunctionPtr f; TaylorSet s; Grid g; uint d; typedef GridTreeSet Result;
    Result run() const { GridTreeSet r = f(s,g,d); return r; }
};

struct ProfileBoxDiscretise {
    ProfileBoxDiscretise(TaylorSet& s_, double e_)  : s(s_), e(e_) { }
    TaylorSet s; double e; typedef ListSet<Box> Result;
    Result run() const { ListSet<Box> r = s.discretise(e); return r; }
};

struct ProfileBoxGridDiscretise {
    ProfileBoxGridDiscretise(TaylorSet& s_, double e_, Grid& g_, uint d_)  : s(s_), e(e_), g(g_), d(d_) { }
    TaylorSet s; double e; Grid g; uint d; typedef GridTreeSet Result;
    Result run() const {
        ListSet<Box> ls = s.discretise(e);
        GridTreeSet r=outer_approximation(ls,g,d);
        return r; }
};


template<class Test>
void profile(const char* name, const Test& test, unsigned int tries)
{

    boost::timer tm; double t=0;

    typename Test::Result res=test.run();
    std::string filename=std::string("profile_taylor_set-")+name;
    Box bounding_box=Box(2, -2.,4., -2.,4.);
    //res.mince(13+res.cell().height());
    uint nc=res.size();
    plot(filename.c_str(),PlanarProjectionMap(2,0,1),bounding_box,
         Colour(1,0,0),res);


    tm.restart();
    for(uint i=0; i!=tries; ++i) {
        test.run();
    }
    t=tm.elapsed();

    

    std::cout << std::fixed
              << std::setw(23) << std::left << name << std::right
              << std::setw(11) << std::setprecision(1) << 1000*(t/tries)
              << std::setw(6) << nc
              << std::setw(9) << std::setprecision(4) << std::fixed << t
              << "/"  << tries
              << std::endl;
}

GridTreeSet discretise1(const TaylorSet&, const Grid&, const uint);
GridTreeSet discretise2(const TaylorSet&, const Grid&, const uint);
GridTreeSet discretise3(const TaylorSet&, const Grid&, const uint);
GridTreeSet discretise4(const TaylorSet&, const Grid&, const uint);

int main(int argc, const char* argv[])
{

    int n=1;

    TaylorSet set(2,2,2, 0.0,1.0,0.25,0.0,0.0,0.0, 0.03125, 0.0,0.5,1.0,1.0,0.0,0.0, 0.0);
    //TaylorSet set(2,2,2, 0.0,1.0,0.25,0.125,-0.375,0.35, 0.03125, 0.0,0.5,1.0,1.0,0.125,-0.125, 0.0);
    std::cerr<<set<<"\n";

    Grid grid(2);

    std::cout << std::setw(23) << std::left << "name" << std::right
              << std::setw(11) << "time(ms)"
              << std::setw(6) << "size"
              << std::setw(9) << "time(s)" << "/" << "tries\n";

/*
    profile("discretise-box-02",ProfileBoxDiscretise(set,1.0/ 2),n*100);
    profile("discretise-box-04",ProfileBoxDiscretise(set,1.0/ 4),n*10);
    profile("discretise-box-06",ProfileBoxDiscretise(set,1.0/ 8),n*10);
    profile("discretise-box-08",ProfileBoxDiscretise(set,1.0/16),n);
    profile("discretise-box-10",ProfileBoxDiscretise(set,1.0/32),n);

    profile("discretise-box-grid-02",ProfileBoxGridDiscretise(set,1.0/ 2,grid, 2),n*10);
    profile("discretise-box-grid-04",ProfileBoxGridDiscretise(set,1.0/ 4,grid, 4),n*5);
    profile("discretise-box-grid-06",ProfileBoxGridDiscretise(set,1.0/ 8,grid, 6),n);
    profile("discretise-box-grid-08",ProfileBoxGridDiscretise(set,1.0/16,grid, 8),n);
    profile("discretise-box-grid-10",ProfileBoxGridDiscretise(set,1.0/32,grid,10),n);

    profile("discretise1-grid-02",ProfileDiscretise(&discretise1,set,grid, 2),n*10);
    profile("discretise1-grid-04",ProfileDiscretise(&discretise1,set,grid, 4),n*5);
    profile("discretise1-grid-06",ProfileDiscretise(&discretise1,set,grid, 6),n);
    profile("discretise1-grid-08",ProfileDiscretise(&discretise1,set,grid, 8),n);
    profile("discretise1-grid-10",ProfileDiscretise(&discretise1,set,grid,10),n);
    profile("discretise2-grid-02",ProfileDiscretise(&discretise2,set,grid, 2),n*10);
    profile("discretise2-grid-04",ProfileDiscretise(&discretise2,set,grid, 4),n*5);
    profile("discretise2-grid-06",ProfileDiscretise(&discretise2,set,grid, 6),n);
    profile("discretise2-grid-08",ProfileDiscretise(&discretise2,set,grid, 8),n);
    profile("discretise2-grid-10",ProfileDiscretise(&discretise2,set,grid,10),n);
    profile("discretise3-grid-02",ProfileDiscretise(&discretise3,set,grid, 2),n*10);
    profile("discretise3-grid-04",ProfileDiscretise(&discretise3,set,grid, 4),n*5);
    profile("discretise3-grid-06",ProfileDiscretise(&discretise3,set,grid, 6),n);
    profile("discretise3-grid-08",ProfileDiscretise(&discretise3,set,grid, 8),n);
    profile("discretise3-grid-10",ProfileDiscretise(&discretise3,set,grid,10),n);
*/

    profile("discretise4-grid-02",ProfileDiscretise(&discretise4,set,grid, 2),n*10);
    profile("discretise4-grid-04",ProfileDiscretise(&discretise4,set,grid, 4),n*5);
    profile("discretise4-grid-06",ProfileDiscretise(&discretise4,set,grid, 6),n);
    profile("discretise4-grid-08",ProfileDiscretise(&discretise4,set,grid, 8),n);
    profile("discretise4-grid-10",ProfileDiscretise(&discretise4,set,grid,10),n);


    return 0;
}

} // namespace Ariadne


int main(int argc, const char* argv[]) {
    Ariadne::main(argc,argv);
}
