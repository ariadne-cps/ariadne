/***************************************************************************
 *            orbit.h
 *
 *  Copyright 2007  Pieter Collins
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
 
/*! \file orbit.h
 *  \brief Orbits of dynamic systems
 */

#ifndef ARIADNE_ORBIT_H
#define ARIADNE_ORBIT_H

#include <utility>
#include <iostream>
#include <boost/shared_ptr.hpp>

namespace Ariadne {

#ifdef DOXYGEN
//! \brief Class for storing evolution data.
template<class E> class Orbit {
  public:
    //! The type used to store a single enclosure set.
    typedef E EnclosureType;
    //! The type used to store a list of enclosure sets.
    typedef EL EnclosureListType;
    //! The initial set of the orbit.
    EnclosureType const& initial() const;
    //! The set of points reached by evolution from the given initial set over the evolution time.
    EnclosureListType const& reach() const;
    //! The set of points reached by evolution from the given initial set at the final evolution time.
    EnclosureListType const& final() const;
};
#endif


template<class ES> class Orbit;

template<class BS> class ListSet;
class GridCell;
class GridTreeSet;
class HybridGridCell;
class HybridGridTreeSet;

template<>
class Orbit<GridCell>
{
    class Data;
  public:
    typedef GridCell EnclosureType;
    typedef GridTreeSet EnclosureListType;

    Orbit(const GridCell&);
    Orbit(const GridCell&, const GridTreeSet&,
          const GridTreeSet&, const GridTreeSet&);
    GridCell const& initial() const;
    GridTreeSet const& reach() const;
    GridTreeSet const& intermediate() const;
    GridTreeSet const& final() const;
  private:
    boost::shared_ptr<Data> _data;
};

template<>
class Orbit<HybridGridCell>
{
    class Data;
  public:
    typedef HybridGridCell EnclosureType;
    typedef HybridGridTreeSet EnclosureListType;

    Orbit(const HybridGridCell&);
    Orbit(const HybridGridCell&, const HybridGridTreeSet&,
          const HybridGridTreeSet&, const HybridGridTreeSet&);
    HybridGridCell const& initial() const;
    HybridGridTreeSet const& reach() const;
    HybridGridTreeSet const& intermediate() const;
    HybridGridTreeSet const& final() const;
  private:
    boost::shared_ptr<Data> _data;
};


typedef int DiscreteState;
class ApproximateTaylorModel;
typedef ApproximateTaylorModel TaylorSetType;
typedef std::pair<DiscreteState,TaylorSetType> HybridTaylorSetType;
typedef ListSet<TaylorSetType> TaylorSetListType;
typedef ListSet<HybridTaylorSetType> HybridTaylorSetListType;

template<class ES>
class Orbit
{
    typedef ListSet<ES> ESL;
  public:
    typedef ES EnclosureType;
    typedef ListSet<ES> EnclosureListType;

    Orbit(const ES& set) : _initial(set) { }
    void adjoin_reach(const EnclosureType& set) { this->_reach.adjoin(set); }
    void adjoin_intermediate(const EnclosureType& set) { this->_intermediate.adjoin(set); }
    void adjoin_final(const EnclosureType& set) { this->_final.adjoin(set); }

    void adjoin_reach(const EnclosureListType& set) { this->_reach.adjoin(set); }
    void adjoin_intermediate(const EnclosureListType& set) { this->_intermediate.adjoin(set); }
    void adjoin_final(const EnclosureListType& set) { this->_final.adjoin(set); }

    EnclosureType const& initial() const { return this->_initial; }
    EnclosureListType const& reach() const { return this->_reach; }
    EnclosureListType const& intermediate() const { return this->_intermediate; }
    EnclosureListType const& final() const { return this->_final; }
  private:
    ES _initial;
    ESL _reach;
    ESL _intermediate;
    ESL _final;
};

template<>
class Orbit<TaylorSetType>
{
    class Data;
    typedef TaylorSetListType list_set_const_iterator;
  public:
    typedef TaylorSetType EnclosureType;
    typedef TaylorSetListType EnclosureListType;

    Orbit(const TaylorSetType&);
    void adjoin_reach(const TaylorSetType& set);
    void adjoin_intermediate(const TaylorSetType& set);
    void adjoin_final(const TaylorSetType& set);

    void adjoin_reach(const TaylorSetListType& set);
    void adjoin_intermediate(const TaylorSetListType& set);
    void adjoin_final(const TaylorSetListType& set);

    TaylorSetType const& initial() const;
    TaylorSetListType const& reach() const;
    TaylorSetListType const& intermediate() const;
    TaylorSetListType const& final() const;
  private:
    boost::shared_ptr<Data> _data;
};

template<>
class Orbit<HybridTaylorSetType>
{
    class Data;
    typedef HybridTaylorSetListType list_set_const_iterator;
  public:
    typedef HybridTaylorSetType EnclosureType;
    typedef HybridTaylorSetListType EnclosureListType;

    Orbit(const HybridTaylorSetType&);
    void adjoin_reach(const HybridTaylorSetType& set);
    void adjoin_intermediate(const HybridTaylorSetType& set);
    void adjoin_final(const HybridTaylorSetType& set);

    void adjoin_reach(const HybridTaylorSetListType& set);
    void adjoin_intermediate(const HybridTaylorSetListType& set);
    void adjoin_final(const HybridTaylorSetListType& set);

    HybridTaylorSetType const& initial() const;
    HybridTaylorSetListType const& reach() const;
    HybridTaylorSetListType const& intermediate() const;
    HybridTaylorSetListType const& final() const;
  private:
    boost::shared_ptr<Data> _data;
};

template<class ES> std::ostream& operator<<(std::ostream& os, const Orbit< ES >& orb);
template<> std::ostream& operator<<(std::ostream& os, const Orbit<TaylorSetType>& orb);


template<class ES> 
std::ostream& 
operator<<(std::ostream& os, const Orbit< ES >& orb)
{
    os << "Orbit(\n  initial=" << orb.initial()
       << "\n  intermediate=" << orb.intermediate()
       << "\n  reach=" << orb.reach()
       << "\n  final=" << orb.final()
       << ")\n";
    return os;
}


template<class G, class ES> void draw(G& graphic, const Orbit<ES>& orbit) 
{
    draw(graphic,orbit.reach()); 
    draw(graphic,orbit.initial());
    draw(graphic,orbit.final());
}

} // namespace Ariadne

#endif // ARIADNE_ORBIT_H
