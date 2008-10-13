/***************************************************************************
 *            grid_set.h
 *
 *  Copyright 2005-8  Alberto Casagrande, Pieter Collins
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
 
/*! \file grid_set.h
 *  \brief Sets based on a coordinate-aligned grid, useful for storage.
 */

#ifndef ARIADNE_GRID_SET_H
#define ARIADNE_GRID_SET_H

#include "macros.h"
#include "tribool.h"
#include "vector.h"

#include "set_interface.h"
#include "list_set.h"

namespace Ariadne {


class EuclideanSpace;

class Box;

template<class Bx> class BoxExpression;

template<class BS> class ListSet;


class Grid;
class GridCell;
class GridCellListSet;
class GridTreeSet;



/*! \ingroup Grid
 *
 *  \brief An infinite, uniform grid of rectangles in Euclidean space.
 *  
 *  \internal Maybe a Grid should be a type of Paving or Cover. 
 *  Then rather than having GridXXX classes, we can have classes such that cells of
 *  some type are mapped into concrete sets by a Paving or Cover.
 *  This should be more general, and will unify the concepts of Paving and Cover,
 *  as well as different types of covers.
 */
class Grid
{ 
 public:
  //! Destructor.
  ~Grid() { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Default constructor constructs a grid from a null pointer. Needed for some iterators.
  explicit Grid() { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Construct from a dimension and a spacing in each direction. 
  explicit Grid(uint d) { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Construct from a dimension and a spacing in each direction. 
  explicit Grid(uint d, Float l) { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Construct from a centre point and a vector of offsets.
  explicit Grid(const Vector<Float>& pt, const Vector<Float>& v) { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Copy constructor. Copies a reference to the grid data.
  Grid(const Grid& g) { ARIADNE_NOT_IMPLEMENTED; }
  
  //! The underlying dimension of the grid.
  uint dimension() const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Tests equality of two grids. Tests equality of references.
  bool operator==(const Grid& g) const { ARIADNE_NOT_IMPLEMENTED; } 
  
  //! The origin of the grid.
  Vector<Float> origin() const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! The strides between successive integer points.
  Vector<Float> lengths() const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Write to an output stream.
  friend std::ostream& operator<<(std::ostream& os, const Grid& g) { ARIADNE_NOT_IMPLEMENTED; }
};



/*! \ingroup Grid 
 *  \ingroup BasicSet
 *
 *  \brief A cell in a grid.
 *  
 *  A %GridCell is defined by mapping a subdivision cell into \f$\mathbb{R}^d\f$
 *  via a grid. 
 *
 *  A %GridCell satisfies the requirements of a RectangleExpression.
 */
class GridCell 
//  : public BoxExpression< GridCell >
{
 public:
  //! The type of denotable real number defining the vertices and cells of the grid.
  typedef Float real_type;
  //! The type of denotable point contained by the set.
  typedef Vector<Float> state_type;
  
  //! Copy constructor.
  GridCell(const GridCell& gc);
  
  //! The grid containing the cell.
  const Grid& grid() const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Convert to a box.
  operator Vector<Interval>() const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! The underlying space of the set.
  uint dimension() const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Convert the cell into a ordinary box.
  Vector<Interval> box() const { ARIADNE_NOT_IMPLEMENTED; }

  //! Split into a "left" or "right" half.
  GridCell split(bool) const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Tests subset.
  friend bool subset(const GridCell& gc1, const GridCell& gc2) { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Tests intersection.
  friend bool overlap(const GridCell& gc1, const GridCell& gc2) { ARIADNE_NOT_IMPLEMENTED; }

  //! Write to an output stream.
  friend std::ostream& operator<<(std::ostream& os, const GridCellListSet& gcls);
};





/*! \ingroup Grid 
 *  \ingroup DenotableSet
 *
 *  \brief A denotable set on a grid, defined as a list of cells.
 *  
 *  A cell list set is useful for storing small sparse sets, such as an 
 *  over-approximation to a basic set.
 *
 */
class GridCellListSet 
{
 public:
  //! The type used for iterating through the cells in the list. 
  typedef std::vector<GridCell>::iterator iterator;
  typedef std::vector<GridCell>::const_iterator const_iterator;
  
  //! A tag describing the type of set.
  //typedef denotable_set_tag set_category;
  //! The type of basic set used in the list.
  typedef GridCell basic_set_type;
  //! The type of denotable real number defining the vertices and cells of the grid.
  typedef Float real_type;
  //! The type of denotable point contained by the set.
  typedef Vector<Float> state_type;
  //! The type of basic set contained by the denotable set.
  typedef GridCell value_type;
  
  //@{
  //! \name Constructors.
  
  //! \brief Virtual destructor. */
  virtual ~GridCellListSet() { }
  //! \brief Construct an empty set based on a Grid.
  explicit GridCellListSet(const Grid& g) { ARIADNE_NOT_IMPLEMENTED; }
  //! Copy constructor.
  GridCellListSet(const GridCellListSet& gcls) { ARIADNE_NOT_IMPLEMENTED; }
  //! Convert from a GridTreeSet.
  GridCellListSet(const GridTreeSet& gts) { ARIADNE_NOT_IMPLEMENTED; }
  //@}
  
  
  //@{ 
  //! \name Conversion operators.
  
  //! Convert to a list of ordinary boxes.
  operator ListSet< Vector<Interval> >() const;
  //@}
  
  
  
  
  //@{ 
  //! \name Access operations.
  
  //! The underlying grid.
  const Grid& grid() const;
  
  //! The numeber of cells in the list.
  uint size() const;
  
  //! The \a i th cell in the list.
  GridCell operator[] (const uint i) const;
  
  //! A constant iterator to the beginning of the list.
  const_iterator begin() const;
  
  //! A constant iterator to the end of the list.
  const_iterator end() const;
  
  //! Tests if the set contains a cell. Works faster on a sorted list. 
  bool contains(GridCell& gc) const;
  
  //! Attempts to find a cell in the list. Works faster on a sorted list. 
  //! This function may not be needed.
  const_iterator find(GridCell& gc) const;
  
  //@}
  
  
  //@{ 
  //! \name Modifying operations
  
  //! Empties the set.
  void clear() { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Removes and returns the last cell in the list.
  GridCell pop() { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Sorts the cells lexicographically, removing duplicates.
  void unique_sort() { ARIADNE_NOT_IMPLEMENTED; } 
  //@}
  
  //@{ 
  //! \name Union, intersection and difference
  
  //! Append the cell \a c to the list.
  void adjoin(const GridCell& gc) { ARIADNE_NOT_IMPLEMENTED; }
  //! Append the cells of \a cls to the list.
  void adjoin(const GridCellListSet& gcls) { ARIADNE_NOT_IMPLEMENTED; }
  //! Adjoins cells contained in \a gts.
  void adjoin(const GridTreeSet& gts) { ARIADNE_NOT_IMPLEMENTED; }
  //! Restricts to cells contained in \a gts.
  void restrict(const GridTreeSet& gts) { ARIADNE_NOT_IMPLEMENTED; }
  //! Removes cells contained in \a bl.
  void remove(const GridTreeSet& gts) { ARIADNE_NOT_IMPLEMENTED; }
  //@}
  
  //@{
  //! \name SetInterface methods
  //! These operators are needed for compliance with the set interface. 
  //! They are probably not all actually needed.
  
  //! Make a dynamically-allocated copy of the set.
  virtual GridCellListSet* clone() const;
  
  //! Returns the denotable set's space.
  virtual uint dimension() const;
  
  //! True if the set is empty.
  tribool empty() const;
  
  //! Checks if a denotable set includes a point.
  virtual tribool contains(const Vector<Float>& p) const;
  
  //! Tests for superset of a box. 
  virtual tribool superset(const Vector<Interval>& r) const;
  
  //! Tests for intersection with a box.
  virtual tribool intersects(const Vector<Interval>& r) const;
  
  //! Tests for disjointness with a box.
  virtual tribool disjoint(const Vector<Interval>& r) const;
  
  //! Tests for subset of a box.
  virtual tribool subset(const Vector<Interval>& r) const;
  
  //! True if the set is bounded.
  tribool bounded() const;
  
  //! The rectangle bounding the region of the mask.
  virtual Vector<Interval> bounding_box() const;
  //@}
  
  //@{ 
  //! \name Approximation methods
  //! Adjoin an over approximation to the set \a bx, computing to depth \a d.
  void adjoin_over_approximation(const Vector<Interval>& bx, uint d) { ARIADNE_NOT_IMPLEMENTED; }
  //! Adjoin an outer approximation to the set \a s, computing to depth \a d.
  void adjoin_outer_approximation(const SetInterface& s, uint d) { ARIADNE_NOT_IMPLEMENTED; }
  //@}

  //@{ 
  //! \name Input/output operators 
  
  //! Write to an output stream.
  friend std::ostream& operator<<(std::ostream& os, const GridCellListSet& gcls) { ARIADNE_NOT_IMPLEMENTED; }
  //@}
};





/*! \ingroup Grid 
 *  \ingroup DenotableSet
 *
 *  \brief A denotable set on a grid, defined using a tree of cells.
 *  
 */
class GridTreeSet 
{
 private:
  class Iterator { };
 public:
  //! A constant iterator through the cells in the set.
  typedef Iterator const_iterator;
  typedef Iterator iterator;
  
  //! A tag describing the type of set.
  //typedef denotable_set_tag set_category;
  //! A tag describing the type of set.
  typedef GridCell basic_set_type;
  
  
  
  //@{
  //! \name Constructors
  
  //! Virtual destructor.
  virtual ~GridTreeSet() { }
  //! Construct an empty set based on a grid.
  GridTreeSet(const Grid& g);
  //! Construct an empty set of dimension \a d.
  GridTreeSet(uint d) { ARIADNE_NOT_IMPLEMENTED; }
  
  //@}
  
  //@{
  //! \name SetInterface methods
  
  //! Tests for disjointness with a box.
  virtual GridTreeSet* clone() const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! The space dimension of the set.
  virtual uint dimension() const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Checks if a denotable set includes a point.
  virtual tribool contains(const Vector<Float>& p) const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Tests for superset of a box. 
  virtual tribool superset(const Vector<Interval>& r) const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Tests for intersection with a box.
  virtual tribool intersects(const Vector<Interval>& r) const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Tests for disjointness with a box.
  virtual tribool disjoint(const Vector<Interval>& r) const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Tests for subset of a box.
  virtual tribool subset(const Vector<Interval>& r) const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! A rectangle containing the grid cell.
  virtual Vector<Interval> bounding_box() const { ARIADNE_NOT_IMPLEMENTED; }

  //! Write to an output stream.
  std::ostream& write(std::ostream&) const { ARIADNE_NOT_IMPLEMENTED; };
  //@}
  
  //@{ 
  //! \name Conversion operators
  
  //! Convert to a ListSet of box.
  operator ListSet< Vector<Interval> > () const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Convert to a GridCellListSet.
  operator GridCellListSet () const { ARIADNE_NOT_IMPLEMENTED; }
  //@}
  
  
  //@{ 
  //! \name Access methods 
  
  //! The underlying grid.
  const Grid& grid() const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! The number of cells in the tree.
  uint capacity() const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! The number of cells in the set.
  uint size() const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! The height of the top cell in the set.
  uint height() const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! The depth of the smallest cell in the set.
  uint depth() const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Constant iterator to the beginning of the cells in the set.
  const_iterator begin() const { ARIADNE_NOT_IMPLEMENTED; }
  //! Constant iterator to the end of the cells in the set.
  const_iterator end() const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Tests if the set is empty.
  tribool empty() const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! True if the set is bounded.
  tribool bounded() const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! An approximation to the volume of the set.
  Float volume() const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! The one-box neighbourhood, consisting of all cells whose closure intersects the set.
  GridTreeSet neighbourhood() const { ARIADNE_NOT_IMPLEMENTED; }
  
  //! The set of all cells which share a facet a cell in the set.
  GridTreeSet adjoining() const { ARIADNE_NOT_IMPLEMENTED; }
  
  //@}
  
  //@{ 
  //! \name Modifying operations.
  
  //! Removes all cells.
  void clear() { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Adjoins a cell to the set.
  void adjoin(const GridCell& gc) { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Adjoins a GridCellListSet to the set.
  void adjoin(const GridCellListSet& gcls) { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Adjoins a GridTreeSet to the set.
  void adjoin(const GridTreeSet& gts) { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Intersects the set with a GridCell.
  void restrict(const GridCell& gc) { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Intersects the set with a GridCellListSet.
  void restrict(const GridCellListSet& gcls) { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Intersects the set with a GridTreeSet.
  void restrict(const GridTreeSet& gts) { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Removes a cell from the set.
  void remove(const GridCell& gc) { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Removes a list of cells from the set.
  void remove(const GridCellListSet& gcls) { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Removes a GridTreeSet from the set.
  void remove(const GridTreeSet& gts) { ARIADNE_NOT_IMPLEMENTED; }
  //@}
  
  //@{ 
  //! \name Approximation-based operations
  
  //! Adjoins an over-approximation of the box \a bx.
  void adjoin_over_approximation(const Vector<Interval>& r) { ARIADNE_NOT_IMPLEMENTED; }
  
  //@}
  
  
  //@{ 
  //! \name Intersection and subset queries
  
  //! Tests subset.
  friend bool subset(GridCell&, GridTreeSet&) { ARIADNE_NOT_IMPLEMENTED; }
  //! Tests subset.
  friend bool subset(GridCellListSet&, GridTreeSet&) { ARIADNE_NOT_IMPLEMENTED; }
  //! Tests subset.
  friend bool subset(GridTreeSet&, GridTreeSet&) { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Tests intersection.
  friend bool overlap(GridCell&, GridTreeSet&) { ARIADNE_NOT_IMPLEMENTED; }
  //! Tests intersection.
  friend bool overlap(GridCellListSet&, GridTreeSet&) { ARIADNE_NOT_IMPLEMENTED; }
  //! Tests intersection.
  friend bool overlap(GridTreeSet&, GridTreeSet&) { ARIADNE_NOT_IMPLEMENTED; }
  //@}
  
  //@{ 
  //! \name Approximation
  
  //! Adjoin to \a gts an outer approximation to the set \a s, computing to depth \a d.
  void adjoin_outer_approximation(const CompactSetInterface& s, uint d) { ARIADNE_NOT_IMPLEMENTED; }
  
  //! Adjoin to \a gts an outer approximation to the set \a s, computing to depth \a d.
  void adjoin_inner_approximation(const OpenSetInterface& s, uint d) { ARIADNE_NOT_IMPLEMENTED; }
  //@}
  
  //@{ 
  //! \name Input/output operators 
  
  //! Write to an output stream.
  friend std::ostream& operator<<(std::ostream& os, const GridTreeSet& gts) { ARIADNE_NOT_IMPLEMENTED; }
  //@}
};

std::ostream& operator<<(std::ostream& os, const GridTreeSet& gts);





} // namespace Ariadne 

#endif /* ARIADNE_GRID_H */

