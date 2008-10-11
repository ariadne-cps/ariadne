/***************************************************************************
 *            test_serialization.cc
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "ariadne.h"
#include "test_float.h"

#include "geometry/grid_mask_set.h"
#include "geometry/hybrid_set.h"
#include "output/serialization.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

template<class R>
class TestSerialization
{
 public:
  int test() {
    test_array();
    test_grid_mask_set();
    test_hybrid_grid_mask_set();
    return 0;
  }

  int test_array() {
    double data[]={5,23,42,111,242};
    const array<double> oary1(data,data+5);
    const array<double> oary2(data,data+3);
    ofstream ofs("test_serialization-array.txt");
    boost::archive::text_oarchive txtoa(ofs);
    txtoa << oary1 << oary2;
    ofs.close();
    
    array<double> iary1(100),iary2(1);
    ifstream ifs("test_serialization-array.txt");
    boost::archive::text_iarchive txtia(ifs);
    txtia >> iary1 >> iary2;
    ofs.close();
    
    ARIADNE_TEST_EQUAL(oary1,iary1);
    ARIADNE_TEST_EQUAL(oary2,iary2);

    return 0;
  }

  int test_grid_mask_set() {
    // Construct a GridMaskSet
    Grid<R> grid(2,0.5);
    Box<R> bbox("[-3,3]x[-2,2]");
    GridMaskSet<R> gms(grid,bbox);
    gms.adjoin_over_approximation(Box<R>("[-0.1,2.9]x[1.4,1.9]"));

    // Save to an output archive. Note that the boost 
    // library insists that the output is a constant object
    const GridMaskSet<R>& ogms=gms;
    ofstream ofs("test_serialization-gms.txt");
    boost::archive::text_oarchive txtoa(ofs);
    txtoa << ogms;
    ofs.close();

    // Load from an input archive. 
    GridMaskSet<R> igms(Grid<R>(1,1.0),Box<R>(1));
    ifstream ifs("test_serialization-gms.txt");
    boost::archive::text_iarchive txtia(ifs);
    txtia >> igms;
    ifs.close();

    ARIADNE_TEST_EQUAL(igms,ogms);
    ARIADNE_ASSERT(igms==ogms);
    std::cout << ogms << std::endl;

    return 1;
  }  
    
  int test_hybrid_grid_mask_set() {
    DiscreteState loc1(1);
    DiscreteState loc2(2);
    GridMaskSet<R> gms1(Grid<R>(2,0.25),Box<R>("[-1,1]x[-1,1]"));
    GridMaskSet<R> gms2(Grid<R>(3,0.5),Box<R>("[-2,0.5]x[-1.5,2.5]"));
    HybridGridMaskSet<R> hgms;
    hgms.new_location(loc1,gms1);
    hgms.new_location(loc2,gms2);
    const HybridGridMaskSet<R>& ohgms=hgms;

    ofstream ofs("test_serialization-hgms.txt");
    boost::archive::text_oarchive txtoa(ofs);
    txtoa << ohgms;
    ofs.close();

    HybridGridMaskSet<R> ihgms;
    ifstream ifs("test_serialization-hgms.txt");
    boost::archive::text_iarchive txtia(ifs);
    txtia >> ihgms;
    ofs.close();

    ARIADNE_TEST_EQUAL(ihgms[loc1],hgms[loc1]);
    ARIADNE_TEST_EQUAL(ihgms[loc2],hgms[loc2]);
    return 0;
  }
};

  
int main() {
  return TestSerialization<Flt>().test();
}
