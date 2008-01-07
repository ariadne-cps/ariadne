/***************************************************************************
 *            partition_tree_set.inline.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
namespace Ariadne {
  namespace Geometry {

    template<class R>
    class PartitionTreeIterator 
      : public boost::iterator_facade<PartitionTreeIterator<R>,
                                      PartitionTreeCell<R>,
                                      boost::forward_traversal_tag,
                                      PartitionTreeCell<R> >
    {
      friend class PartitionTree<R>;
     public:
      PartitionTreeIterator(const Box<R>& bb, 
                            const Combinatoric::SubdivisionSequence& ss, 
                            Combinatoric::BinaryTree::const_iterator i)
        : _bounding_box(bb), _subdivisions(ss), _base(i) { }
     private:
      bool equal(const PartitionTreeIterator<R>& other) const {
        return this->_base == other._base; }
      void increment() { ++_base; }
      PartitionTreeCell<R> dereference() const { 
        return PartitionTreeCell<R>(_bounding_box,_subdivisions,*_base); }
     private:
      const Box<R> _bounding_box;
      const Combinatoric::SubdivisionSequence _subdivisions;
      Combinatoric::BinaryTree::const_iterator _base;
    };

    template<class R>
    class PartitionTreeSetIterator 
      : public boost::iterator_facade<PartitionTreeSetIterator<R>,
                                      PartitionTreeCell<R>,
                                      boost::forward_traversal_tag,
                                      PartitionTreeCell<R> >
    {
      friend class PartitionTree<R>;
     public:
     private:
      bool equal(const PartitionTreeSetIterator<R>& other) const {
        return this->_base == other._base; }
      void increment() { ++_base; }
      PartitionTreeCell<R> dereference() const { 
        return PartitionTreeCell<R>(_bounding_box,_subdivisions,*_base); }
     private:
      const Box<R> _bounding_box;
      const Combinatoric::SubdivisionSequence _subdivisions;
      Combinatoric::MaskedBinaryTree::const_iterator _base;
    };

  } // namespace Geometry
} // namespace Ariadne





namespace Ariadne {

template<class R>
Geometry::PartitionScheme<R>::PartitionScheme(const Box<R>& bb, const Combinatoric::SubdivisionSequence& sc)
  : _unit_box(bb), _subdivisions(sc)
{
}


template<class R>
bool 
Geometry::PartitionScheme<R>::operator==(const PartitionScheme<R>& pg) const 
{
  return _unit_box==pg._unit_box && _subdivisions==pg._subdivisions;
}

template<class R> inline
bool 
Geometry::PartitionScheme<R>::operator!=(const PartitionScheme<R>& pg) const 
{
  return !(*this==pg); 
}

template<class R> inline
const Geometry::Box<R>& 
Geometry::PartitionScheme<R>::unit_box() const 
{
  return _unit_box; 
}

template<class R> inline
const Combinatoric::SubdivisionSequence& 
Geometry::PartitionScheme<R>::subdivisions() const 
{
  return _subdivisions; 
}

template<class R> inline
dimension_type 
Geometry::PartitionScheme<R>::dimension() const 
{
  return _subdivisions.dimension(); 
}



template<class R> inline
Geometry::PartitionTreeCell<R>::PartitionTreeCell(const Box<R>& r, const Combinatoric::SubdivisionCell& c)
  : _unit_box(r), _subdivision_cell(c)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(r,c,"PartitionTreeCell::PartitionTreeCell(Box r, SubdivisionCell c)");
}

template<class R> inline
Geometry::PartitionTreeCell<R>::PartitionTreeCell(const Box<R>& r, 
                                                  const Combinatoric::SubdivisionSequence& s, 
                                                  const Combinatoric::BinaryWord& w) 
  : _unit_box(r), _subdivision_cell(s,w) 
{ 
}

template<class R> inline
const Geometry::Box<R>& 
Geometry::PartitionTreeCell<R>::unit_box() const 
{
  return this->_unit_box; 
}

template<class R> inline
const Combinatoric::SubdivisionCell& 
Geometry::PartitionTreeCell<R>::subdivision_cell() const 
{
  return this->_subdivision_cell; 
}

template<class R> inline
dimension_type 
Geometry::PartitionTreeCell<R>::dimension() const 
{
  return this->_subdivision_cell.dimension(); 
}

template<class R> inline
tribool 
Geometry::PartitionTreeCell<R>::empty() const 
{
  return false; 
}

template<class R> inline
tribool 
Geometry::PartitionTreeCell<R>::bounded() const 
{
  return true; 
}





template<class R> inline
Geometry::PartitionTree<R>::PartitionTree(const Box<R>& r, 
                                          const Combinatoric::SubdivisionSequence& s, 
                                          const Combinatoric::BinaryTree& t)
  : _unit_box(r), _subdivision_tree(s,t) 
{
}

template<class R> inline
Geometry::PartitionTree<R>::PartitionTree(const PartitionScheme<R>& ps, const Combinatoric::BinaryTree& t)
  : _unit_box(ps.unit_box()), _subdivision_tree(ps.subdivisions(),t) 
{
}

template<class R> inline
const Geometry::Box<R>& 
Geometry::PartitionTree<R>::unit_box() const 
{
  return _unit_box; 
}

template<class R> inline
const Combinatoric::SubdivisionTree& 
Geometry::PartitionTree<R>::subdivision_tree() const 
{
  return _subdivision_tree; 
}

template<class R> inline
dimension_type 
Geometry::PartitionTree<R>::dimension() const 
{
  return _subdivision_tree.dimension(); 
}

template<class R> inline
const Combinatoric::SubdivisionSequence& 
Geometry::PartitionTree<R>::subdivisions() const 
{
  return _subdivision_tree.subdivisions(); 
}

template<class R> inline
const Combinatoric::BinaryTree& 
Geometry::PartitionTree<R>::binary_tree() const 
{
  return _subdivision_tree.binary_tree(); 
}

template<class R> inline
size_type 
Geometry::PartitionTree<R>::size() const 
{
  return _subdivision_tree.size(); 
}

template<class R> inline
Geometry::PartitionScheme<R> 
Geometry::PartitionTree<R>::scheme() const 
{
  return  PartitionScheme<R>(unit_box(),subdivisions()); 
}

template<class R> inline
typename Geometry::PartitionTree<R>::const_iterator 
Geometry::PartitionTree<R>::begin() const 
{
  return const_iterator(_unit_box,_subdivision_tree.begin()); 
}
template<class R> inline
typename Geometry::PartitionTree<R>::const_iterator 
Geometry::PartitionTree<R>::end() const 
{
  return const_iterator(_unit_box,_subdivision_tree.end()); 
}





template<class R> inline
Geometry::PartitionTreeSet<R>::PartitionTreeSet(const PartitionScheme<R>& g)
  : _unit_box(g.unit_box()), _subdivision_set(g.subdivisions()) 
{
}

template<class R> inline
Geometry::PartitionTreeSet<R>::PartitionTreeSet(const PartitionScheme<R>& g, const Combinatoric::BinaryTree& t, const BooleanArray& m)
  : _unit_box(g.unit_box()), _subdivision_set(g.subdivisions(),t,m) 
{
}

template<class R> inline
Geometry::PartitionTreeSet<R>::PartitionTreeSet(const PartitionTree<R>& t, const BooleanArray& m)
  : _unit_box(t.unit_box()), _subdivision_set(t.subdivisions(),t.binary_tree(),m)
{
}

template<class R> inline
Geometry::PartitionTreeSet<R>::PartitionTreeSet(const Box<R>& r, 
                                                const Combinatoric::SubdivisionSequence& s, 
                                                const Combinatoric::BinaryTree& t, 
                                                const BooleanArray& m)
  : _unit_box(r), _subdivision_set(s,t,m)
{
}


template<class R> inline
Geometry::Box<R> 
Geometry::PartitionTreeSet<R>::bounding_box() const 
{
  return _unit_box; 
}

template<class R> inline
const Geometry::Box<R>& 
Geometry::PartitionTreeSet<R>::unit_box() const 
{
  return _unit_box; 
}

template<class R> inline
const Combinatoric::SubdivisionTreeSet& 
Geometry::PartitionTreeSet<R>::subdivision_set() const 
{
  return _subdivision_set; 
}

template<class R> inline
dimension_type 
Geometry::PartitionTreeSet<R>::dimension() const 
{
  return _subdivision_set.dimension(); 
}

template<class R> inline
const Combinatoric::SubdivisionSequence& 
Geometry::PartitionTreeSet<R>::subdivisions() const 
{
  return _subdivision_set.subdivisions(); 
}

template<class R> inline
const Combinatoric::BinaryTree& 
Geometry::PartitionTreeSet<R>::binary_tree() const 
{
  return _subdivision_set.binary_tree(); 
}

template<class R> inline
const BooleanArray& 
Geometry::PartitionTreeSet<R>::mask() const 
{
  return _subdivision_set.mask(); 
}

template<class R> inline
size_type 
Geometry::PartitionTreeSet<R>::capacity() const 
{
  return _subdivision_set.capacity(); 
}

template<class R> inline
size_type 
Geometry::PartitionTreeSet<R>::size() const 
{
  return _subdivision_set.size(); 
}

template<class R> inline
SizeArray 
Geometry::PartitionTreeSet<R>::depths() const 
{
  return _subdivision_set.depths(); 
}

template<class R> inline
size_type 
Geometry::PartitionTreeSet<R>::depth() const 
{
  return _subdivision_set.depth(); 
}

template<class R> inline
Geometry::PartitionScheme<R> 
Geometry::PartitionTreeSet<R>::scheme() const 
{
  return PartitionScheme<R>(bounding_box(),subdivisions()); 
}

template<class R> inline
Geometry::PartitionTree<R> 
Geometry::PartitionTreeSet<R>::partition_tree() const 
{
  return PartitionTree<R>(bounding_box(),subdivisions(),binary_tree()); 
}

template<class R> inline
typename Geometry::PartitionTreeSet<R>::const_iterator 
Geometry::PartitionTreeSet<R>::begin() const 
{
  return const_iterator(_unit_box,_subdivision_set.begin()); 
}

template<class R> inline
typename Geometry::PartitionTreeSet<R>::const_iterator 
Geometry::PartitionTreeSet<R>::end() const 
{
  return const_iterator(_unit_box,_subdivision_set.end()); 
}

template<class R> inline
tribool 
Geometry::PartitionTreeSet<R>::empty() const 
{
  return this->size()==0; 
}

template<class R> inline
tribool 
Geometry::PartitionTreeSet<R>::bounded() const 
{
  return true; 
}






template<class R> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const PartitionScheme<R>& ps) 
{ 
  return ps.write(os);
}

template<class R> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const PartitionTreeCell<R>& ptc)
{ 
  return ptc.write(os);
}

template<class R> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const PartitionTree<R>& pt) 
{ 
  return pt.write(os);
}

template<class R> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const PartitionTreeSet<R>& pts)
{ 
  return pts.write(os);
}


} // namespace Ariadne

