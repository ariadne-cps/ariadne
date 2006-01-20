/***************************************************************************
 *            quadtree.h
 *
 *  Sun Jan 23 15:26:27 2005
 *  Copyright  2005  Alberto Casagrande
 *  Email casagrande@dimi.uniud.it
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
 
#ifndef _PTREE_H
#define _PTREE_H

#include <rectangle.h>
#include <ptree_node.h>

namespace Ariadne {
namespace Geometry {

template <typename R>
std::ostream& operator<<(std::ostream &os, PartitionTree<R> &t) {
  
  os << "Bounding Box=" << t._bounding_box << std::endl;
  os << "Max Depth=" << t._max_depth;
  
  for (size_t i=0; i<t.size();i++){
    os << std::endl<< "Rectangle["<<i<<"]="<< 
      t._vector[i];
  }
  
  //os << "Tree=" << t._root << std::endl;
  
  return os;
}



template <typename R>
void pov_print(PartitionTree<R> &t, std::string f_name, 
               size_t dim1, size_t dim2, size_t dim3) 
{
  std::ofstream fos;
  typedef typename R::State State;
  
  State lower(1),upper(1);
  
  f_name+=".pov";
  
  fos.open(f_name.c_str() , std::ios::out);
  
  fos << "merge {" << std::endl;
  
  for (size_t i=0; i<t.size();i++){
    lower=t[i].lower_corner();
    upper=t[i].upper_corner();
    
    //mpf_class l[3], u[3];

    fos << std::endl << "box {" << std::endl 
        << "< " << (mpf_class) lower[dim1] << " , " << (mpf_class) lower[dim2] << " , " << (mpf_class) lower[dim3]<<" >," << std::endl 
        << "< " << (mpf_class) upper [dim1] << " , " << (mpf_class) upper[dim2] << " , "  << (mpf_class) upper[dim3]<<" >" << std::endl 
        << "texture {" << std::endl 
        << "pigment { color rgb 0.4 }" << std::endl 
        << "}" << std::endl 
        << "}" << std::endl;
  }
  fos << "}" << std::endl;
  
  fos.close();
}

/*! \brief A partition tree representation for gridded sets.
 *
 * This class represents partition trees.
 */ 
template <typename R = Rectangle < State < Rational > > >
class PartitionTree {
 public:
  typedef R Rectangle;
  typedef typename R::State State;
  typedef typename R::State::Real Real;
  
  typedef typename std::vector<Rectangle>::const_iterator const_iterator;
  
 private:
  /*! \brief The ptree's bounding box. */
  Rectangle _bounding_box;
  
  /*! \brief The ptree's root. */
  PartitionTreeNode _root;
  
  /*! \brief The ptree's maximum depth */
  size_t _max_depth;
  
  /*! \brief The number of rectangles into the ptree. */
  //size_t _size;
  
  /*! \brief List of maintained rectangles. */
  std::vector<Rectangle> _vector;
  
  uint _count_rect(const PartitionTreeNode &node,
                   const Rectangle &node_box) {
    
    if (node.full()) {
      return 1;
    }
    
    if (node.empty())
      return 0;
    
    uint rect_n=0;
    
    for (size_t i=0; i< node._bits_per_space_member(); i++) {
      if ((node._space).test(i)) {
        
        Rectangle new_box=
          node_box.find_quadrant(i);
        
        rect_n+=this->_count_rect(*(node._sons[i]),
                                  new_box);
      }
    }
    
    return rect_n;
  }
  
  void _fill_rect_vector(const PartitionTreeNode &node,
                         const Rectangle &node_box, size_t &idx) {
    
    if (node.full()) {
      this->_vector[idx]=node_box;
      idx++;
      return;
    }
    
    if (node.empty())
      return;
    
    for (size_t i=0; i< node._bits_per_space_member(); i++) {
      if ((node._space).test(i)) {
        Rectangle new_box=
          node_box.find_quadrant(i);
        
        this->_fill_rect_vector(
                                *(node._sons[i]),
                                new_box, idx);
      }
    }
  }
  
  void _fill_rect_vector(const PartitionTreeNode &node,
                         const Rectangle &node_box) {
    
    uint i=0,rect_n=this->_count_rect(node,node_box);
    
    (this->_vector).resize(rect_n);
    
    this->_fill_rect_vector(node,node_box,i);
  }
  
  template <typename SetType>
  inline void _fill_ptree_check_subtrees(PartitionTreeNode &node,
                                         const SetType &A, const Rectangle &bbox, 
                                         const size_t &depth,const size_t &max_depth) {
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    for (size_t i=0; i< node._bits_per_space_member(); i++) {
      
      Rectangle B=bbox.find_quadrant(i);
      
      if (interiors_intersect(B,A)) {
        PartitionTreeNode *son= 
          new PartitionTreeNode(A.dimension());
        
        this->_fill_ptree_with(*son , A, B, 
                               depth+1 , max_depth);
        
        if (son->empty()) {
          delete son;
          node._sons[i]=NULL;
        } else {
          (node._space).set(i);
          node._sons[i]=son;
        }
        
      }
      
    }
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    
  }
  
  
  template <typename SetType>
  void _fill_ptree_with(PartitionTreeNode &node,
                        const SetType &A, 
                        const Rectangle &bbox, 
                        size_t depth,const size_t &max_depth) {
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    node._depth=depth;
    
    if (interior_subset_of_interior(bbox,A)) {
      
      /* set full */
      node.set_full();
      
#ifdef DEBUG
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
      
      return;
    }
    
    if (depth>max_depth) {
      
      /* default behaviour is to under-approximate */
      node.set_empty();
      
#ifdef DEBUG
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
      
      return;
    }
    
    this->_fill_ptree_check_subtrees(node, A, bbox, 
                                     depth, max_depth);
    
    if ((node._sons).size()==0) {
      node.set_empty();
      
#ifdef DEBUG
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
      
      return;
    }
    
    node.set_not_empty();
    
    node.eval_level_by_sons();
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
  }
  
  template <typename SetType>
  inline void _refill_ptree_check_subtrees(PartitionTreeNode &node,
                                           const SetType &A, const Rectangle &bbox, 
                                           const size_t &depth,const size_t &max_depth) {
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    for (size_t i=0; i< node._bits_per_space_member(); i++) {
      
      Rectangle B=bbox.find_quadrant(i);
      
      if ((node._space).test(i)) {
        
        if (interiors_intersect(B,A)) {
          
          this->_refill_ptree_with( *(node._sons[i]),
                                    A, B, depth+1 , max_depth);
        }
        
      } else {
        
        if (interiors_intersect(B,A)) {
          PartitionTreeNode *son= 
            new PartitionTreeNode(A.dimension());
          
          this->_refill_ptree_with(*son , A, B, 
                                   depth+1 , max_depth);
          
          if (son->empty()) {
            delete son;
            node._sons[i]=NULL;
          } else {
            (node._space).set(i);
            node._sons[i]=son;
          }
        }
        
      }
      
    }
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
  }
  
  template <typename SetType>
  void _refill_ptree_with(PartitionTreeNode &node,
                          const SetType &A, 
                          const Rectangle &bbox, 
                          size_t depth,const size_t &max_depth) {
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    if (node.full()) {
      /* if the node is full we don't need to
       * add anything */
      
#ifdef DEBUG
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
      
      return;
    }
    
    node._depth=depth;
    
    if (interior_subset_of_interior(bbox,A)) {
      
      /* set full */
      node.set_full();
      
#ifdef DEBUG
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
      
      return;
    }
    
    if (depth>max_depth) {
      
      /* default behaviour is to under-approximate */
      node.set_empty();
      
#ifdef DEBUG
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
      
      return;
    }
    
    this->_refill_ptree_check_subtrees(node, A, bbox, 
                                       depth, max_depth);
    
    if ((node._sons).size()==0) {
      node.set_empty();
      
#ifdef DEBUG
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
      
      return;
    }
    
    node.set_not_empty();
    
    node.compact_node();
    
    node.eval_level_by_sons();
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
  }
  
  void _clear_vector(){
    (this->_vector).clear();
  }
  
 public:
  
  /*! \brief A class costructor.
   *
   * Creates an PartitionTree. 
   */
  PartitionTree() {}
  
  /*! \brief A class costructor.
   *
   * Creates an PartitionTree.
   * \param bbox is the bounding box of the ptree.
   * \param m_depth is the maximum depth of the ptree. 
   */
  PartitionTree(const Rectangle &bbox, const size_t m_depth): 
    _bounding_box(bbox), _root(bbox.dimension()),
    _max_depth(m_depth){}
  
  /*! \brief A class costructor.
   *
   * Creates an PartitionTree.
   * \param A is the original Ptree.
   */
  PartitionTree(const PartitionTree<Rectangle> &A): 
    _bounding_box(A._bounding_box), _root(A._root),
    _max_depth(A._max_depth){}
  
  /*! \brief A class costructor. 
   *
   * Creates an PartitionTree and inserts an object.
   * \param A is an object that should be inserted into the ptree.
   * \param bbox is the bounding box of the ptree.
   * \param m_depth is the maximum depth of the ptree.
   */
  template <typename SetType>
  PartitionTree(const SetType &A,
                const Rectangle &bbox, const size_t m_depth):  
    _bounding_box(bbox), _root(bbox.dimension()), 
    _max_depth(m_depth) {
    
    if (A.dimension()!=bbox.dimension()) 
      throw std::domain_error("Two of the parameters have different space dimensions");
    
    if (!interior_subset_of_interior(A, this->_bounding_box))
      throw std::domain_error("The object can not be contained into the ptree.");
    
    this->_fill_ptree_with(this->_root, A, bbox, 0, m_depth);
    
  }
  
  /*! \brief A class destructor. */
  ~PartitionTree() {
    this->_clear_vector();
  }
  
  /*! \brief Return the number of rectangles into the PartitionTree. 
   *
   * \return The number of rectangles maitained by the PartitionTree.
   */
  inline size_t size() {
    if (this->_vector.size()==0) {
      this->_fill_rect_vector(this->_root,this->_bounding_box);
    }
    
    return (this->_vector).size();
  }
  
  /*! \brief Return the space dimension of the PartitionTree. 
   *
   * \return The space dimension of the PartitionTree.
   */
  inline size_t dimension() const{
    return (this->_root).dimension();
  }
  
  /*! \brief Return the space dimension of the PartitionTree. (Deprecated) */
  inline size_t dim() const{
    return (this->_root).dimension();
  }
  
  /*! \brief Accesses the i-th rectangle.
   *
   * \param index is the index of the returned rectangle.
   * \return The i-th rectangle maitained by the PartitionTree.
   */
  inline const Rectangle &operator[](size_t index){
    if (this->_vector.size()==0) {
      this->_fill_rect_vector(this->_root,this->_bounding_box);
    }
    
    return this->_vector[index];
  }
  
  /*! \brief Checks the inclusion of a state.
   *
   * \param state is a state in the space represented by the PartitionTree.
   * \return \a true if the state is contained into the rectangles maintained 
   * by the PartitionTree, \a false otherwise.
   */
  inline bool contains(const State & state) {
    if (this->_vector.size()==0) {
      this->_fill_rect_vector(this->_root,this->_bounding_box);
    }
    
    for (size_t i=0; i<this->size();i++){
      if ((this->_vector[i]).contains(state))
        return true;
    }
    
    return false;
    
  }
  
  /*! \brief Checks the interior inclusion of a state.
   *
   * \param state is a state in the space represented by the PartitionTree.
   * \return \a true if the state is contained into the interior of 
   * the rectangles maintained by the PartitionTree, \a false otherwise.
   */
  inline bool interior_contains(const State & state) const {
    
    /* TO REIMPLEMENT */
    return false;
    
  }
  
  /*! \brief Checks the emptyness.
   *
   * \return \a true if the space represented by the PartitionTree is empty,
   * \a false otherwise.
   */
  inline bool empty() const {
    return ((this->_root).empty());
  }
  
  /*! \brief Returns the begin of the maintained rectangle vector.
   *
   * \return The begin of the maintained rectangle vector.
   */
  inline const_iterator begin() {
    if (this->_vector.size()==0) {
      this->_fill_rect_vector(this->_root,this->_bounding_box);
    }
    
    return _vector().begin();
  }
  
  /*! \brief Returns the end of the maintained rectangle vector.
   *
   * \return The end of the maintained rectangle vector.
   */
  inline const_iterator end() {
    if (this->_vector.size()==0) {
      this->_fill_rect_vector(this->_root,this->_bounding_box);
    }
    
    return _vector().end();
  }
  
  /*! \brief Makes the union of two PartitionTree.
   *
   * Makes the union of the current PartitionTree and an other PartitionTree.
   * The result is stored into the current object.
   * \param \a A is an PartitionTree.
   */
  inline void inplace_union(const PartitionTree& A) {
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    if (this->dimension()!=A.dimension()) 
      throw std::domain_error("The two parameters have different space dimensions.");
    
    if (this->_bounding_box!=A._bounding_box)
      throw std::domain_error("The two parameters have different bounding boxes.");
    
    (this->_root).inplace_union(A._root);
    
    this->_clear_vector();
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
  }
  
  /*! \brief Inserts an object into the PartitionTree.
   *
   * Inserts an object into the PartitionTree. The result is stored 
   * into the current object.
   * \param \a A is an object representing a set.
   */
  template <typename SetType>
  inline void inplace_union(const SetType &A) {
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    if (A.empty()) return; 
    
    if (this->dimension()!=A.dimension()) 
      throw std::domain_error("The two parameters have different space dimensions.");
    
    if (this->_bounding_box!=A._bounding_box)
      throw std::domain_error("The two parameters have different bounding boxes.");
    
    this->_refill_ptree_with(this->_root, A, this->_bounding_box , 0, 
                             this->_max_depth );
    
    this->_clear_vector();
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
  }
  
  /*! \brief Inserts an object into the PartitionTree.
   *
   * Inserts an object into the PartitionTree. The result is stored 
   * into the current object.
   * \param \a A is an object representing a set.
   */
  template <typename SetType>
  inline void inplace_union_bbox_intersection(const SetType &A) {
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    if (A.empty()) return; 
    
    if (this->dimension()!=A.dimension()) 
      throw std::domain_error("The two parameters have different space dimensions.");
    
    this->_refill_ptree_with(this->_root, A, this->_bounding_box , 0, 
                             this->_max_depth );
    
    this->_clear_vector();
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
  }
  
  inline const Rectangle &bounding_box() const { return this->_bounding_box; }
  
  inline const PartitionTree<R> &operator=(const PartitionTree<R> &A) {
    this->_bounding_box=A._bounding_box;
    this->_root=A._root;
    this->_max_depth=A._max_depth;
    
    return *this;
  }
  
  /*! \brief Prints the ptree on a ostream */
  template <typename Rect>
  friend std::ostream& operator<<(std::ostream &os, 
                                  PartitionTree<Rect> &r);
};
  
}
}
#endif /* _PTREE_H */
