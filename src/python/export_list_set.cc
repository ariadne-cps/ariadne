/***************************************************************************
 *            python/export_list_set.cc
 *
 *  21 October 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */



#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/polyhedron.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"

#include "python/typedefs.h"
#include "python/python_utilities.h"
using namespace Ariadne;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
using namespace boost::python;

inline RParallelotope plsg(const RParallelotopeListSet& s, int n) { return ::get_item(s,n); }
inline RZonotope zlsg(const RZonotopeListSet& s, int n) { return ::get_item(s,n); }

inline 
bool interiors_intersect(const RZonotopeListSet& A, const RParallelotope& B) 
{
  return interiors_intersect(A,RZonotopeListSet(RZonotope(B)));
}

template <typename R, template<typename> class BS, template<typename> class BS2>
ListSet<R,BS> 
touching_intersection(const ListSet<R,BS>& ls, const BS2<R>& bs) 
{
  ListSet<R,BS> output(ls.dimension());
    
  for (size_t i=0; i< ls.size(); i++) {
    if (!(disjoint(ls[i],bs)))
      output.adjoin(ls[i]);
  }

  return output;
}

template ListSet<Real,Rectangle> touching_intersection(const ListSet<Real,Rectangle>&, const Rectangle<Real> &);
template ListSet<Real,Parallelotope> touching_intersection(const ListSet<Real,Parallelotope>&, const Rectangle<Real> &);
template ListSet<Real,Zonotope> touching_intersection(const ListSet<Real,Zonotope>&, const Rectangle<Real> &);
template ListSet<Real,Rectangle> touching_intersection(const ListSet<Real,Rectangle>&, const Parallelotope<Real> &);
template ListSet<Real,Parallelotope> touching_intersection(const ListSet<Real,Parallelotope>&, const Parallelotope<Real> &);
template ListSet<Real,Zonotope> touching_intersection(const ListSet<Real,Zonotope>&, const Parallelotope<Real> &);
template ListSet<Real,Rectangle> touching_intersection(const ListSet<Real,Rectangle>&, const Zonotope<Real> &);
template ListSet<Real,Parallelotope> touching_intersection(const ListSet<Real,Parallelotope>&, const Zonotope<Real> &);
template ListSet<Real,Zonotope> touching_intersection(const ListSet<Real,Zonotope>&, const Zonotope<Real> &);


void export_list_set() {
  typedef bool (*RectLSBinPred) (const RRectangleListSet&, const RRectangleListSet&);
  typedef bool (*ParLSBinPred) (const RParallelotopeListSet&, const RParallelotopeListSet&);
  typedef bool (*ZonLSBinPred) (const RZonotopeListSet&, const RZonotopeListSet&);
  typedef bool (*ParLSParBinPred) (const RParallelotope&, const RParallelotopeListSet&);
  typedef bool (*LSParParBinPred) (const RParallelotopeListSet&, const RParallelotope&);
  typedef bool (*ZonLSZonBinPred) (const RZonotope&, const RZonotopeListSet&);
  typedef bool (*LSZonZonBinPred) (const RZonotopeListSet&, const RZonotope&);
  typedef bool (*LSZonParBinPred) (const RZonotopeListSet&, const RParallelotope&);
  
   typedef RRectangleListSet (*LTRlsrBinFun) (const RRectangleListSet&, const RRectangle&);
   typedef RRectangleListSet (*LTRlspBinFun) (const RRectangleListSet&, const RParallelotope&);
   typedef RRectangleListSet (*LTRlszBinFun) (const RRectangleListSet&, const RZonotope&);
   typedef RParallelotopeListSet (*LTPlsrBinFun) (const RParallelotopeListSet&, const RRectangle&);
   typedef RParallelotopeListSet (*LTPlspBinFun) (const RParallelotopeListSet&, const RParallelotope&);
   typedef RParallelotopeListSet (*LTPlszBinFun) (const RParallelotopeListSet&, const RZonotope&);
   typedef RZonotopeListSet (*LTZlsrBinFun) (const RZonotopeListSet&, const RRectangle&);
   typedef RZonotopeListSet (*LTZlspBinFun) (const RZonotopeListSet&, const RParallelotope&);
   typedef RZonotopeListSet (*LTZlszBinFun) (const RZonotopeListSet&, const RZonotope&);
  
  typedef RRectangleListSet (*RectLSBinFun) (const RRectangleListSet&, const RRectangleListSet&);
  typedef void (RRectangleListSet::* RectLSadjRectLS) (const RRectangleListSet&);
  typedef void (RParallelotopeListSet::* PltpLSadjRect) (const RRectangle&);
  typedef void (RParallelotopeListSet::* PltpLSadjRectLS) (const RRectangleListSet&);
  typedef void (RParallelotopeListSet::* PltpLSadjPltp) (const RParallelotope&);
  typedef void (RParallelotopeListSet::* PltpLSadjPltpLS) (const RParallelotopeListSet&);
  typedef void (RZonotopeListSet::* ZltzLSadjZltz) (const RZonotope&);
  typedef void (RZonotopeListSet::* ZltzLSadjZltzLZ) (const RZonotopeListSet&);


def("regular_intersection", RectLSBinFun(&regular_intersection));
  def("interiors_intersect", RectLSBinPred(&interiors_intersect));
  def("interiors_intersect", ParLSBinPred(&interiors_intersect));
  def("interiors_intersect", ZonLSBinPred(&interiors_intersect));
  //def("interiors_intersect", ParLSParBinPred(&interiors_intersect));
  //def("interiors_intersect", LSParParBinPred(&interiors_intersect));
  //def("interiors_intersect", LSZonZonBinPred(&interiors_intersect));
  //def("interiors_intersect", LSZonParBinPred(&interiors_intersect));
  def("touching_intersection", LTRlsrBinFun(&touching_intersection));
  def("touching_intersection", LTRlspBinFun(&touching_intersection));
  def("touching_intersection", LTRlszBinFun(&touching_intersection));
  def("touching_intersection", LTPlsrBinFun(&touching_intersection));
  def("touching_intersection", LTPlspBinFun(&touching_intersection));
  def("touching_intersection", LTPlszBinFun(&touching_intersection));
  def("touching_intersection", LTZlsrBinFun(&touching_intersection));
  def("touching_intersection", LTZlspBinFun(&touching_intersection));
  def("touching_intersection", LTZlszBinFun(&touching_intersection));
  def("disjoint", RectLSBinPred(&disjoint));
  def("inner_subset", RectLSBinPred(&inner_subset));
  def("subset", RectLSBinPred(&subset));

  class_<RRectangleListSet>("RectangleListSet",init<int>())
    .def(init<RRectangle>())
    .def(init<RRectangleListSet>())
    .def(init<RGridBlockListSet>())
    .def(init<RGridCellListSet>())
    .def(init<RGridMaskSet>())
    .def(init<RPartitionTreeSet>())
    .def("dimension", &RRectangleListSet::dimension)
    .def("push_back", &RRectangleListSet::push_back)
    .def("size", &RRectangleListSet::size)
    .def("empty", &RRectangleListSet::empty)
    .def("__len__", &RRectangleListSet::size)
    .def("__getitem__", &RRectangleListSet::get, return_value_policy<copy_const_reference>())
    .def("__setitem__", &RRectangleListSet::set)
    .def("__iter__", iterator<RRectangleListSet>())
    .def(self_ns::str(self))    // __self_ns::str__
  ;
  
  class_<RParallelotopeListSet>("ParallelotopeListSet",init<int>())
    .def(init<RParallelotope>())
    .def(init<RRectangleListSet>())
    .def(init<RParallelotopeListSet>())
    .def("dimension", &RParallelotopeListSet::dimension)
    .def("push_back", &RParallelotopeListSet::push_back)
    .def("adjoin", PltpLSadjPltp(&RParallelotopeListSet::adjoin))
    .def("adjoin", PltpLSadjPltpLS(&RParallelotopeListSet::adjoin))
    .def("size", &RParallelotopeListSet::size)
    .def("empty", &RParallelotopeListSet::empty)
    .def("__len__", &RParallelotopeListSet::size)
//    .def("__getitem__", &RParallelotopeListSet::get, return_value_policy<copy_const_reference>())
//    .def("__getitem__", &get_item<RParallelotopeListSet>)
    .def("__getitem__", &plsg)
    .def("__setitem__", &RParallelotopeListSet::set)
    .def("__iter__", iterator<RParallelotopeListSet>())
    .def(self_ns::str(self))    // __self_ns::str__
  ;
  
  class_<RZonotopeListSet>("ZonotopeListSet",init<int>())
    .def(init<RZonotope>())
    .def(init<RRectangleListSet>())
    .def(init<RParallelotopeListSet>())
    .def(init<RZonotopeListSet>())
    .def("dimension", &RZonotopeListSet::dimension)
    .def("push_back", &RZonotopeListSet::push_back)
    .def("adjoin", ZltzLSadjZltz(&RZonotopeListSet::adjoin))
    .def("adjoin", ZltzLSadjZltzLZ(&RZonotopeListSet::adjoin))
    .def("size", &RZonotopeListSet::size)
    .def("empty", &RZonotopeListSet::empty)
    .def("__len__", &RZonotopeListSet::size)
//    .def("__getitem__", &RParallelotopeListSet::get, return_value_policy<copy_const_reference>())
//    .def("__getitem__", &get_item<RParallelotopeListSet>)
    .def("__getitem__", &zlsg)
    .def("__setitem__", &RZonotopeListSet::set)
    .def("__iter__", iterator<RZonotopeListSet>())
    .def(self_ns::str(self))    // __self_ns::str__
  ;

}
