/***************************************************************************
 *            python/export_approximator.cc
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
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "python/float.h"

#include "geometry/box.h"
#include "geometry/zonotope.h"
#include "geometry/taylor_set.h"
#include "geometry/box_list_set.h"
#include "geometry/grid_cell_list_set.h"
#include "geometry/grid_mask_set.h"
#include "evaluation/approximator_interface.h"
#include "evaluation/standard_approximator.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class Aprx, class ES>
class ApproximatorWrapper
  : public ApproximatorInterface<Aprx,ES>,
    public wrapper< ApproximatorInterface<Aprx,ES> >
{
  typedef typename ES::real_type R;
  typedef Interval<R> I;
  typedef typename Aprx::Paving Paving;
  typedef typename Aprx::BasicSet BasicSet;
  typedef typename Aprx::CoverListSet CoverListSet;
  typedef typename Aprx::PartitionListSet PartitionListSet;
  typedef typename Aprx::PartitionTreeSet PartitionTreeSet;
  typedef ES EnclosureSet;
  typedef ListSet<ES> EnclosureSetList;
  
 public:
  ApproximatorWrapper<Aprx,ES>* clone() const { return this->get_override("clone")(); }
  EnclosureSet enclosure_set(const BasicSet&) const { return this->get_override("enclosure_set")(); }
  R radius(const EnclosureSet&) const { return this->get_override("radius")(); }
  tribool superset(const EnclosureSet&, const BasicSet&) const { return this->get_override("superset")(); }
  tribool disjoint(const EnclosureSet&, const BasicSet&) const { return this->get_override("disjoint")(); }
  BasicSet bounding_box(const EnclosureSet&) const { return this->get_override("bounding_box")(); }
  CoverListSet over_approximation(const EnclosureSet&) const { return this->get_override("over_approximation")(); }
  PartitionListSet outer_approximation(const EnclosureSet&, const Paving&) const { return this->get_override("outer_approximation")(); }

  BasicSet bounding_box(const EnclosureSetList&) const { return this->get_override("bounding_box")(); }
  CoverListSet over_approximations(const EnclosureSetList&) const { return this->get_override("ovwer_approximations")(); }
  PartitionListSet outer_approximation(const EnclosureSetList&, const Paving&) const { return this->get_override("outer_approximation")(); }

  void adjoin_over_approximation(CoverListSet&, const EnclosureSet&) const { this->get_override("adjoin_over_approximation")(); }
  void adjoin_outer_approximation(PartitionListSet&, const EnclosureSet&) const { this->get_override("adjoin_outer_approximation")(); }
  void adjoin_outer_approximation(PartitionTreeSet&, const EnclosureSet&) const { this->get_override("adjoin_outer_approximation")(); }

  void adjoin_over_approximations(CoverListSet&, const EnclosureSetList&) const { this->get_override("adjoin_over_approximations")(); }
  void adjoin_outer_approximation(PartitionListSet&, const EnclosureSetList&) const { this->get_override("adjoin_outer_approximation")(); }
  void adjoin_outer_approximation(PartitionTreeSet&, const EnclosureSetList&) const { this->get_override("adjoin_outer_approximation")(); }
  std::ostream& write(std::ostream& os) const { this->get_override("write")(); return os; }
};


template<class R>
void export_approximator() 
{
  typedef GridApproximationScheme<R> GAS;
  typedef Grid<R> Pv;
  typedef Box<R> Bx;
  typedef Zonotope<R> ZES;

  class_< ApproximatorWrapper<GAS,ZES>, boost::noncopyable >("ZonotopeApproximatorInterface",init< >());

  class_< StandardApproximator<ZES>, bases< ApproximatorInterface<GAS,ZES> > >
    approximator_class("StandardApproximator",init<>());
  approximator_class.def("enclosure_set",&StandardApproximator<ZES>::enclosure_set);
  approximator_class.def("radius",&StandardApproximator<ZES>::radius);
  approximator_class.def("bounding_box",&StandardApproximator<ZES>::bounding_box);
  approximator_class.def("over_approximation",&StandardApproximator<ZES>::over_approximation);
  approximator_class.def("outer_approximation",(GridCellListSet<R>(StandardApproximator<ZES>::*)(const ZES&,const Pv&)const)&StandardApproximator<ZES>::outer_approximation);
  approximator_class.def("adjoin_outer_approximation",(void(StandardApproximator<ZES>::*)(GridCellListSet<R>&,const ZES&)const)&StandardApproximator<ZES>::adjoin_outer_approximation);

  
}

template void export_approximator<FloatPy>();
