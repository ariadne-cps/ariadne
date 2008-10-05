#ifndef ARIADNE_GRID_H
#define ARIADNE_GRID_H

#include "tribool.h"
#include "vector.h"

namespace Ariadne {


class basic_set_tag;

class EuclideanSpace;
class SetInterface;

/*! \brief A coordinate-aligned box in Euclidean space. */
class Box;
/*! \brief %Base class for expressions exposing the box interface. */
template<class Bx> class BoxExpression;
/*! \brief A union of sets of type \a BS, implemented as an ordered list. */
template<class BS> class ListSet {};

/*! \brief A topological partition or paving of a space by sets of type \a BS. */
template<class BS> class Paving { };
/*! \brief An open cover of a space by sets of type \a BS. */
template<class BS> class Cover { };


class Grid;
class Cell;
class GridListSet;
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
  ~Grid();
  
  //! Default constructor constructs a grid from a null pointer. Needed for some iterators.
  explicit Grid();
  
  //! Construct from a dimension and a spacing in each direction. 
  explicit Grid(uint d);
  
  //! Construct from a dimension and a spacing in each direction. 
  explicit Grid(uint d, Float l);
  
  //! Construct from a centre point and a vector of offsets.
  explicit Grid(const Vector<Float>& pt, const Vector<Float>& v);
  
  //! Copy constructor. Copies a reference to the grid data.
  Grid(const Grid& g);
  
  //! The underlying dimension of the grid.
  uint dimension() const;
  
  //! Tests equality of two grids. Tests equality of references.
  bool operator==(const Grid& g) const; 
  
  //! The origin of the grid.
  Vector<Float> origin() const;
  
  //! The strides between successive integer points.
  Vector<Float> lengths() const;
  
  //! Compute an over-approximation of a box on the grid.
  //! Maybe this should be a friend... or a constructor of Cell?
  Cell over_approximation(const Vector<Interval>& bx) const;
  
  
  //! Write to an output stream.
  std::ostream& write(std::ostream& os) const;
  //! Read from an input stream.
  std::istream& read(std::istream& is);
  
  //! Write to an output stream.
  friend std::ostream& operator<<(std::ostream& os, const Grid& g);
};



/*! \ingroup Grid 
 *  \ingroup BasicSet
 *
 *  \brief A cell in a grid.
 *  
 *  A %Cell is defined by mapping a subdivision cell into \f$\mathbb{R}^d\f$
 *  via a grid. 
 *
 *  A %Cell satisfies the requirements of a RectangleExpression.
 */
class Cell 
//  : public BoxExpression< Cell >
{
 public:
  //! A tag describing the type of set.
  typedef basic_set_tag set_category;
  //! The type of denotable real number defining the vertices and cells of the grid.
  typedef Float real_type;
  //! The type of denotable point contained by the set.
  typedef Vector<Float> state_type;
  
  //! Copy constructor.
  Cell(const Cell& gc);
  
  //! The grid containing the cell.
  const Grid& grid() const;
  
  //! Convert to a box.
  operator Vector<Interval>() const;
  
  //! The underlying space of the set.
  uint dimension() const;
  
  //! Convert the cell into a ordinary box.
  Vector<Interval> box() const;

  //! Split into a "left" or "right" half.
  Cell split(bool) const;
  
  //! Write to an output stream.
  std::ostream& write(std::ostream&) const;
  
  //! Tests subset.
  friend bool subset(const Cell& gc1, const Cell& gc2);
  
  //! Tests intersection.
  friend bool overlap(const Cell& gc1, const Cell& gc2);
  
  //! Write to an output stream.
  friend std::ostream& operator<<(std::ostream& os, const Cell& gc);
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
class CellList 
{
 private:
  class Iterator;
 public:
  //! The type used for iterating through the cells in the list. 
  typedef Iterator const_iterator;
  typedef Iterator iterator;
  
  //! A tag describing the type of set.
  //typedef denotable_set_tag set_category;
  //! The type of basic set used in the list.
  typedef Cell basic_set_type;
  //! The type of denotable real number defining the vertices and cells of the grid.
  typedef Float real_type;
  //! The type of denotable point contained by the set.
  typedef Vector<Float> state_type;
  //! The type of basic set contained by the denotable set.
  typedef Cell value_type;
  
  //@{
  //! \name Constructors.
  
  //! \brief Virtual destructor. */
  virtual ~CellList() { }
  //! \brief Construct an empty set based on a Grid.
  explicit CellList(const Grid& g);
  //! Copy constructor.
  CellList(const CellList& gcls);
  //! Convert from a GridTreeSet.
  CellList(const GridTreeSet& gts);
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
  Cell operator[] (const uint i) const;
  
  //! A constant iterator to the beginning of the list.
  const_iterator begin() const;
  
  //! A constant iterator to the end of the list.
  const_iterator end() const;
  
  //! Tests if the set contains a cell. Works faster on a sorted list. 
  bool contains(Cell& gc) const;
  
  //! Attempts to find a cell in the list. Works faster on a sorted list. 
  //! This function may not be needed.
  const_iterator find(Cell& gc) const;
  
  //@}
  
  
  //@{ 
  //! \name Modifying operations
  
  //! Empties the set.
  void clear();
  
  //! Removes and returns the last cell in the list.
  Cell pop();
  
  //! Sorts the cells lexicographically, removing duplicates.
  void unique_sort(); 
  //@}
  
  //@{ 
  //! \name Union, intersection and difference
  
  //! Append the cell \a c to the list.
  void adjoin(const Cell& gc);
  //! Append the cells of \a cls to the list.
  void adjoin(const CellList& gcls);
  //! Adjoins cells contained in \a gts.
  void adjoin(const GridTreeSet& gts);
  //! Restricts to cells contained in \a gts.
  void restrict(const GridTreeSet& gts);
  //! Removes cells contained in \a bl.
  void remove(const GridTreeSet& gts);
  //@}
  
  //@{
  //! \name SetInterface methods
  //! These operators are needed for compliance with the set interface. 
  //! They are probably not all actually needed.
  
  //! Make a dynamically-allocated copy of the set.
  virtual CellList* clone() const;
  
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
  void adjoin_over_approximation(const Vector<Interval>& bx, uint d);
  //! Adjoin an outer approximation to the set \a s, computing to depth \a d.
  void adjoin_outer_approximation(const SetInterface& s, uint d);
  //@}

  //@{ 
  //! \name Input/output operators
  
  //! Write to an output stream.
  std::ostream& write(std::ostream&) const;
  //@}
  
  //@{ 
  //! \name Input/output operators 
  
  //! Write to an output stream.
  friend std::ostream& operator<<(std::ostream& os, const CellList& gcls);
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
  class Iterator;
 public:
  //! A constant iterator through the cells in the set.
  typedef Iterator const_iterator;
  typedef Iterator iterator;
  
  //! A tag describing the type of set.
  //typedef denotable_set_tag set_category;
  //! A tag describing the type of set.
  typedef Cell basic_set_type;
  
  
  
  //@{
  //! \name Constructors
  
  //! Virtual destructor.
  virtual ~GridTreeSet() { }
  //! Construct an empty set based on a grid.
  GridTreeSet(const Grid& g);
  //! Construct an empty set of dimension \a d.
  GridTreeSet(uint d);
  
  //@}
  
  //@{
  //! \name SetInterface methods
  
  //! Tests for disjointness with a box.
  virtual GridTreeSet* clone() const;
  
  //! The space dimension of the set.
  virtual uint dimension() const;
  
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
  
  //! A rectangle containing the grid cell.
  virtual Vector<Interval> bounding_box() const;
  //@}
  
  //@{ 
  //! \name Conversion operators
  
  //! Convert to a ListSet of box.
  operator ListSet< Vector<Interval> > () const;
  
  //! Convert to a CellList.
  operator CellList () const;
  //@}
  
  
  //@{ 
  //! \name Access methods 
  
  //! The underlying grid.
  const Grid& grid() const;
  
  //! The number of cells in the tree.
  uint capacity() const;
  
  //! The number of cells in the set.
  uint size() const;
  
  //! The height of the top cell in the set.
  uint height() const;
  
  //! The depth of the smallest cell in the set.
  uint depth() const;
  
  //! Constant iterator to the beginning of the cells in the set.
  const_iterator begin() const;
  //! Constant iterator to the end of the cells in the set.
  const_iterator end() const;
  
  //! Tests if the set is empty.
  tribool empty() const;
  
  //! True if the set is bounded.
  tribool bounded() const;
  
  //! An approximation to the volume of the set.
  Float volume() const;
  
  //! The one-box neighbourhood, consisting of all cells whose closure intersects the set.
  GridTreeSet neighbourhood() const;
  
  //! The set of all cells which share a facet a cell in the set.
  GridTreeSet adjoining() const;
  
  //@}
  
  //@{ 
  //! \name Modifying operations.
  
  //! Removes all cells.
  void clear();
  
  //! Adjoins a cell to the set.
  void adjoin(const Cell& gc);
  
  //! Adjoins a CellList to the set.
  void adjoin(const CellList& gcls);
  
  //! Adjoins a GridTreeSet to the set.
  void adjoin(const GridTreeSet& gts);
  
  //! Intersects the set with a Cell.
  void restrict(const Cell& gc);
  
  //! Intersects the set with a CellList.
  void restrict(const CellList& gcls);
  
  //! Intersects the set with a GridTreeSet.
  void restrict(const GridTreeSet& gts);
  
  //! Removes a cell from the set.
  void remove(const Cell& gc);
  
  //! Removes a list of cells from the set.
  void remove(const CellList& gcls);
  
  //! Removes a GridTreeSet from the set.
  void remove(const GridTreeSet& gts);
  //@}
  
  //@{ 
  //! \name Approximation-based operations
  
  //! Adjoins an over-approximation of the box \a bx.
  void adjoin_over_approximation(const Vector<Interval>& r);
  
  //@}
  
  //@{ 
  //! \name Input/output operators 
  
  //! Writes a summary of the set.
  std::string summary() const;
  
  //! Write to an output stream.
  std::ostream& write(std::ostream&) const;
  
  //@}
  
  //@{ 
  //! \name Intersection and subset queries
  
  //! Tests subset.
  friend bool subset(Cell&, GridTreeSet&);
  //! Tests subset.
  friend bool subset(CellList&, GridTreeSet&);
  //! Tests subset.
  friend bool subset(GridTreeSet&, GridTreeSet&);
  
  //! Tests intersection.
  friend bool overlap(Cell&, GridTreeSet&);
  //! Tests intersection.
  friend bool overlap(CellList&, GridTreeSet&);
  //! Tests intersection.
  friend bool overlap(GridTreeSet&, GridTreeSet&);
  //@}
  
  //@{ 
  //! \name Approximation
  
  //! Adjoin to \a gts an outer approximation to the set \a s, computing to depth \a d.
  void adjoin_outer_approximation(const SetInterface& s, uint d);
  
  //! Adjoin to \a gts an outer approximation to the set \a s, computing to depth \a d.
  void adjoin_inner_approximation(const SetInterface& s, uint d);
  //@}
  
  //@{ 
  //! \name Input/output operators 
  
  //! Write to an output stream.
  friend std::ostream& operator<<(std::ostream& os, const GridTreeSet& gts);
  //@}
};







} // namespace Ariadne 

#endif /* ARIADNE_GRID_H */

