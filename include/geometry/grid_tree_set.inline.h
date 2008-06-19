/***************************************************************************
 *            partition_tree_set.inline.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
                            const SubdivisionSequence& ss, 
                            BinaryTree::const_iterator i)
        : _bounding_box(bb), _subdivisions(ss), _base(i) { }
     private:
      bool equal(const PartitionTreeIterator<R>& other) const {
        return this->_base == other._base; }
      void increment() { ++_base; }
      PartitionTreeCell<R> dereference() const { 
        return PartitionTreeCell<R>(_bounding_box,_subdivisions,*_base); }
     private:
      const Box<R> _bounding_box;
      const SubdivisionSequence _subdivisions;
      BinaryTree::const_iterator _base;
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
      const SubdivisionSequence _subdivisions;
      MaskedBinaryTree::const_iterator _base;
    };

  } // namespace Geometry
} // namespace Ariadne





namespace Ariadne {

template<class R>
PartitionScheme<R>::PartitionScheme(const Box<R>& bb, const SubdivisionSequence& sc)
  : _unit_box(bb), _subdivisions(sc)
{
}


template<class R>
bool 
PartitionScheme<R>::operator==(const PartitionScheme<R>& pg) const 
{
  return _unit_box==pg._unit_box && _subdivisions==pg._subdivisions;
}

template<class R> inline
bool 
PartitionScheme<R>::operator!=(const PartitionScheme<R>& pg) const 
{
  return !(*this==pg); 
}

template<class R> inline
const Box<R>& 
PartitionScheme<R>::unit_box() const 
{
  return _unit_box; 
}

template<class R> inline
const SubdivisionSequence& 
PartitionScheme<R>::subdivisions() const 
{
  return _subdivisions; 
}

template<class R> inline
dimension_type 
PartitionScheme<R>::dimension() const 
{
  return _subdivisions.dimension(); 
}



template<class R> inline
PartitionTreeCell<R>::PartitionTreeCell(const Box<R>& r, const SubdivisionTreeCell& c)
  : _unit_box(r), _subdivision_cell(c)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(r,c,"PartitionTreeCell::PartitionTreeCell(Box r, SubdivisionTreeCell c)");
}

template<class R> inline
PartitionTreeCell<R>::PartitionTreeCell(const Box<R>& r, 
                                                  const SubdivisionSequence& s, 
                                                  const BinaryWord& w) 
  : _unit_box(r), _subdivision_cell(s,w) 
{ 
}

template<class R> inline
const Box<R>& 
PartitionTreeCell<R>::unit_box() const 
{
  return this->_unit_box; 
}

template<class R> inline
const SubdivisionTreeCell& 
PartitionTreeCell<R>::subdivision_cell() const 
{
  return this->_subdivision_cell; 
}

template<class R> inline
dimension_type 
PartitionTreeCell<R>::dimension() const 
{
  return this->_subdivision_cell.dimension(); 
}

template<class R> inline
tribool 
PartitionTreeCell<R>::empty() const 
{
  return false; 
}

template<class R> inline
tribool 
PartitionTreeCell<R>::bounded() const 
{
  return true; 
}





template<class R> inline
PartitionTree<R>::PartitionTree(const Box<R>& r, 
                                          const SubdivisionSequence& s, 
                                          const BinaryTree& t)
  : _unit_box(r), _subdivision_tree(s,t) 
{
}

template<class R> inline
PartitionTree<R>::PartitionTree(const PartitionScheme<R>& ps, const BinaryTree& t)
  : _unit_box(ps.unit_box()), _subdivision_tree(ps.subdivisions(),t) 
{
}

template<class R> inline
const Box<R>& 
PartitionTree<R>::unit_box() const 
{
  return _unit_box; 
}

template<class R> inline
const SubdivisionTree& 
PartitionTree<R>::subdivision_tree() const 
{
  return _subdivision_tree; 
}

template<class R> inline
dimension_type 
PartitionTree<R>::dimension() const 
{
  return _subdivision_tree.dimension(); 
}

template<class R> inline
const SubdivisionSequence& 
PartitionTree<R>::subdivisions() const 
{
  return _subdivision_tree.subdivisions(); 
}

template<class R> inline
const BinaryTree& 
PartitionTree<R>::binary_tree() const 
{
  return _subdivision_tree.binary_tree(); 
}

template<class R> inline
size_type 
PartitionTree<R>::size() const 
{
  return _subdivision_tree.size(); 
}

template<class R> inline
PartitionScheme<R> 
PartitionTree<R>::scheme() const 
{
  return  PartitionScheme<R>(unit_box(),subdivisions()); 
}

template<class R> inline
typename PartitionTree<R>::const_iterator 
PartitionTree<R>::begin() const 
{
  return const_iterator(_unit_box,_subdivision_tree.begin()); 
}
template<class R> inline
typename PartitionTree<R>::const_iterator 
PartitionTree<R>::end() const 
{
  return const_iterator(_unit_box,_subdivision_tree.end()); 
}





template<class R> inline
PartitionTreeSet<R>::PartitionTreeSet(const PartitionScheme<R>& g)
  : _unit_box(g.unit_box()), _subdivision_set(g.subdivisions()) 
{
}

template<class R> inline
PartitionTreeSet<R>::PartitionTreeSet(const PartitionScheme<R>& g, const BinaryTree& t, const BooleanArray& m)
  : _unit_box(g.unit_box()), _subdivision_set(g.subdivisions(),t,m) 
{
}

template<class R> inline
PartitionTreeSet<R>::PartitionTreeSet(const PartitionTree<R>& t, const BooleanArray& m)
  : _unit_box(t.unit_box()), _subdivision_set(t.subdivisions(),t.binary_tree(),m)
{
}

template<class R> inline
PartitionTreeSet<R>::PartitionTreeSet(const Box<R>& r, 
                                                const SubdivisionSequence& s, 
                                                const BinaryTree& t, 
                                                const BooleanArray& m)
  : _unit_box(r), _subdivision_set(s,t,m)
{
}


template<class R> inline
Box<R> 
PartitionTreeSet<R>::bounding_box() const 
{
  return _unit_box; 
}

template<class R> inline
const Box<R>& 
PartitionTreeSet<R>::unit_box() const 
{
  return _unit_box; 
}

template<class R> inline
const SubdivisionTreeSet& 
PartitionTreeSet<R>::subdivision_set() const 
{
  return _subdivision_set; 
}

template<class R> inline
dimension_type 
PartitionTreeSet<R>::dimension() const 
{
  return _subdivision_set.dimension(); 
}

template<class R> inline
const SubdivisionSequence& 
PartitionTreeSet<R>::subdivisions() const 
{
  return _subdivision_set.subdivisions(); 
}

template<class R> inline
const BinaryTree& 
PartitionTreeSet<R>::binary_tree() const 
{
  return _subdivision_set.binary_tree(); 
}

template<class R> inline
const BooleanArray& 
PartitionTreeSet<R>::mask() const 
{
  return _subdivision_set.mask(); 
}

template<class R> inline
size_type 
PartitionTreeSet<R>::capacity() const 
{
  return _subdivision_set.capacity(); 
}

template<class R> inline
size_type 
PartitionTreeSet<R>::size() const 
{
  return _subdivision_set.size(); 
}

template<class R> inline
SizeArray 
PartitionTreeSet<R>::depths() const 
{
  return _subdivision_set.depths(); 
}

template<class R> inline
size_type 
PartitionTreeSet<R>::depth() const 
{
  return _subdivision_set.depth(); 
}

template<class R> inline
PartitionScheme<R> 
PartitionTreeSet<R>::scheme() const 
{
  return PartitionScheme<R>(bounding_box(),subdivisions()); 
}

template<class R> inline
PartitionTree<R> 
PartitionTreeSet<R>::partition_tree() const 
{
  return PartitionTree<R>(bounding_box(),subdivisions(),binary_tree()); 
}

template<class R> inline
typename PartitionTreeSet<R>::const_iterator 
PartitionTreeSet<R>::begin() const 
{
  return const_iterator(_unit_box,_subdivision_set.begin()); 
}

template<class R> inline
typename PartitionTreeSet<R>::const_iterator 
PartitionTreeSet<R>::end() const 
{
  return const_iterator(_unit_box,_subdivision_set.end()); 
}

template<class R> inline
tribool 
PartitionTreeSet<R>::empty() const 
{
  return this->size()==0; 
}

template<class R> inline
tribool 
PartitionTreeSet<R>::bounded() const 
{
  return true; 
}






template<class R> inline
std::ostream& 
operator<<(std::ostream& os, const PartitionScheme<R>& ps) 
{ 
  return ps.write(os);
}

template<class R> inline
std::ostream& 
operator<<(std::ostream& os, const PartitionTreeCell<R>& ptc)
{ 
  return ptc.write(os);
}

template<class R> inline
std::ostream& 
operator<<(std::ostream& os, const PartitionTree<R>& pt) 
{ 
  return pt.write(os);
}

template<class R> inline
std::ostream& 
operator<<(std::ostream& os, const PartitionTreeSet<R>& pts)
{ 
  return pts.write(os);
}


} // namespace Ariadne

