/***************************************************************************
 *            grid_set.h
 *
 *  Copyright  2008  Ivan S. Zapreev, Pieter Collins
 *            ivan.zapreev@gmail.com, pieter.collins@cwi.nl
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
 *  Foundation, Inc., 59 Templece Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file grid_set.h
 *  \brief Grid paving is used to represent sets, based on integer and dyadic coordinate cells, of a grid.
 */

#ifndef ARIADNE_GRID_SET_H
#define ARIADNE_GRID_SET_H

#include <iostream>
#include <string>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/shared_ptr.hpp>

#include "tribool.h"
#include "array.h"

#include "binary_word.h"

#include "exceptions.h"
#include "box.h"
#include "point.h"
#include "list_set.h"

#include "numeric.h"

#include "set_interface.h"
#include "vector.h"
#include "grid.h"

#include "graphics_interface.h"

using namespace std;
using namespace Ariadne;

namespace Ariadne {

/*Some type definitions*/
typedef std::vector<bool> BooleanArray;
typedef array<int> IndexArray;
typedef array<unsigned int> SizeArray;

typedef unsigned short dimension_type;

/*Some pre-declarations*/
class BinaryTreeNode;
class Grid;
class GridAbstractCell;
class GridCell;
class GridOpenCell;
class GridTreeSubset;
class GridTreeSet;

class GridTreeCursor;
class GridTreeConstIterator;
    
/*Declarations of classes in other files*/
template<class BS> class ListSet;

std::ostream& operator<<(std::ostream& output_stream, const BinaryTreeNode & binary_tree );
std::ostream& operator<<(std::ostream& os, const GridCell& theGridCell);
std::ostream& operator<<(std::ostream& os, const GridOpenCell& theGridOpenCell );
std::ostream& operator<<(std::ostream& os, const GridTreeCursor& theGridTreeCursor);
std::ostream& operator<<(std::ostream& os, const GridTreeSubset& theGridTreeSubset);
std::ostream& operator<<(std::ostream& os, const GridTreeSet& theGridTreeSet);

bool subset( const GridCell& theCellOne, const GridCell& theCellTwo, BinaryWord * pPathPrefixOne = NULL,
             BinaryWord * pPathPrefixTwo = NULL, uint * pPrimaryCellHeight = NULL );
bool subset(const GridCell& theCell, const GridTreeSubset& theSet);
bool overlap(const GridCell& theCell, const GridTreeSubset& theSet);
bool subset(const GridTreeSubset& theSet1, const GridTreeSubset& theSet2);
bool overlap(const GridTreeSubset& theSet1, const GridTreeSubset& theSet2);
    
GridTreeSet join(const GridTreeSubset& theSet1, const GridTreeSubset& theSet2);
GridTreeSet intersection(const GridTreeSubset& theSet1, const GridTreeSubset& theSet2);
GridTreeSet difference(const GridTreeSubset& theSet1, const GridTreeSubset& theSet2);
    
GridTreeSet outer_approximation(const Box& theBox, const Grid& theGrid, const uint depth);
GridTreeSet outer_approximation(const CompactSetInterface& theSet, const Grid& theGrid, const uint depth);
GridTreeSet outer_approximation(const CompactSetInterface& theSet, const uint depth);
template<class BS> GridTreeSet outer_approximation(const ListSet<BS>& theSet, const uint depth);

template<class A> void serialize(A& archive, const GridTreeSet& set, const uint version);


/*! \brief The binary-tree node operation is not allowed on a non-leaf node. */
class NotALeafNodeException : public std::logic_error {
  public: 
    NotALeafNodeException(const std::string& str) : std::logic_error(str) { }
};
    
/*! \brief The binary-tree node operation is not allowed on a leaf node. */
class IsALeafNodeException : public std::logic_error {
  public: 
    IsALeafNodeException(const std::string& str) : std::logic_error(str) { }
};
    
/*! \brief The GridTreeCursor throws this exception if we try to go beyond the binary tree. */
class NotAllowedMoveException : public std::logic_error {
  public: 
    NotAllowedMoveException(const std::string& str) : std::logic_error(str) { }
};


/*! \brief The binary tree node.
 *
 * This node is to be used in a binary tree designed for subdividing the state
 * space into enabled and disabled cells. This is required for representing
 * subsets of the state space.
 * 
 * \b Storage: We only store pointers to the left and right subtrees and the tribool
 * value indicating whether this cell is enabled/disabled or we do not know.
 */
class BinaryTreeNode {
  protected: 
    /*! \brief Defines whether the given node of the tree is on/off or we do not know*/
    tribool _isEnabled;
        
    /*! \brief The left and right subnodes of the tree. Note that,
     * \a pLeftNode == \a NULL iff \a pRightNode == \a NULL. The latter
     * is allowed iff \a _isEnabled == \a TRUE or \a _isEnabled == \a FALSE,
     * i.e. the node can be a leaf.
     */
    BinaryTreeNode* _pLeftNode;
    BinaryTreeNode* _pRightNode;

    /*! \brief This method splits the enabled subtrees of the tree rooted to
     * \a pCurrentNode in such a way that the tree depth becomes \a depth.
     * If the initial tree depth is greater than \a depth then nothing is done.
     */
    void mince_node(BinaryTreeNode* pCurrentNode, const uint depth);
            
    /*! \brief This method recombined the sub tree nodes rooted to \a pCurrentNode.
     * Note that, the two leaf nodes with the same parent are removed if they have
     * the same value of isEnabled fields.
     */
    void recombine_node(BinaryTreeNode * pCurrentNode);
            
    /*! \brief This method is used for recursive restoration of the binary tree from
     *  the used data, i.e. \a theTree and \a theEnabledCells. It is used in the
     *  constructor \a BinaryTreeNode( const BooleanArray& , const BooleanArray& )
     */
    void restore_node( BinaryTreeNode * pCurrentNode, uint & arr_index, uint & leaf_counter,
                       const BooleanArray& theTree, const BooleanArray& theEnabledCells);
            
    /*! \brief This method is used in constructors for the node initialization */
    void init( tribool isEnabled, BinaryTreeNode* pLeftNode, BinaryTreeNode* pRightNode );
            
  public:
    //@{
    //! \name Constructors
            
    /*! \brief Construct a tree node. */
    explicit BinaryTreeNode(const tribool _isEnabled = false );
            
    /*! \brief The copy constructor.
     * The the node and all it's sub nodes are copied.
     */
    explicit BinaryTreeNode(const BinaryTreeNode& theTreeNode);
            
    /*! \brief Constructs a binary tree from the boolean arrays. \a theTree defines
     * the tree structure, \a theEnabledCells defines the tree nodes.
     * IVAN S. ZAPREEV:
     * WARNING: We assume that every node has either no children or both of them!
     * NOTE: The dinary data of the tree is organized in the following way:
     *          1      The \a theTree contains Depth first search lay out of the tree,
     *         / \     where 1 stands for the non-leaf node and 0 for a leaf node, we 
     *        /   \    always visit the left su-node first, e.g. the tree on the left
     *       2     5   is encodes the array:     [1, 1, 0, 0, 0]
     *      / \        where the corresponding    ^  ^  ^  ^  ^
     *     /   \       tree nodes are:            1  2  3  4  5
     *    3     4      \a theEnabledCells contains true/false values for the leaf nodes
     * of the tree. Their order is the same as in \a theTree, e.g. here it is: 3,4,5
     */
    explicit BinaryTreeNode( const BooleanArray& theTree, const BooleanArray& theEnabledCells );
            
    //@}
            
    ~BinaryTreeNode();
            
    //@{
    //! \name Properties
            
    /*! \brief Returns true if the node is marked as enabled, otherwise false */
    bool is_enabled() const;

    /*! \brief Returns true if some of the leaf nodes in the tree rooted to this node are enabled, otherwise false */
    bool has_enabled() const;

    /*! \brief Returns true if all leaf nodes in the tree rooted to this node are enabled, otherwise false */
    bool all_enabled() const;

    /*! \brief This method returns true if the given path defines a node in the tree that is either enabled
     *  or is in a "virtual" subtree of some enabled node (note that enabled nodes can only be leafs).
     * Note that: \a path is treated as if it is rooted to this node, we assume that the path starts
     * from position \a position of \a path. The parameter \a position is used for recursive calls only.
     */
    bool is_enabled( const BinaryWord & path, const uint position = 0) const;
            
    /*! \brief Returns true if the node is marked as disabled, otherwise false */
    bool is_disabled() const;
            
    /*! \brief Returns true if the node is a leaf (pLeftNode == NULL && pRightNode == NULL) otherwise false */
    bool is_leaf() const;
            
    /*! \brief Return the left sub-node */
    BinaryTreeNode * left_node() const;
            
    /*! \brief Return the right sub-node */
    BinaryTreeNode * right_node() const;

    /*! \brief Returns the depth of the sub-tree rooted to the given node, i.e. the depth of it's deepest node */
    uint depth() const;
    
    /*! \brief Allows to compare to binaty tree nodes */
    bool operator==(const BinaryTreeNode & otherNode ) const;
            
    static bool is_equal_nodes( const BinaryTreeNode * pFirstNode, const BinaryTreeNode * pSecondNode );
            
    //@}
            
    //@{
    //! \name Leaf Operations

    /*! \brief This method makes the node to become a leaf node with the enabled value : \a is_enabled
     * NOTE: the leat and the right sub-trees (are deallocated.
     * WARNING: this method MUST NOT be called on a non-leaf node!!!
     */
    void make_leaf( tribool is_enabled );

    /*! \brief Marks the leaf node as enabled, otherwise they through \a NotALeafNodeEsception */
    void set_enabled();
            
    /*! \brief Marks the leaf node as disabled, otherwise they through \a NotALeafNodeEsception */
    void set_disabled();
            
    /*! \brief Marks the node as neither enabled nor disabled, is only applicable to non-leaf nodes.
     * When applied to a leaf node, throws IsALeafNodeException.
     */
    void set_unknown();
            
    /*! \brief Splits the leaf node, i.e. adds two subnodes with the _isEnabled field value inherited from the parent node.
     * The parent's _isEnabled is set to intermediate, because the subsequent operation on the subtree might enable/disable
     * subnodes and we do not want to keep track of these changes. If the node is not a leaf then nothing is done.
     */
    void split();
            
    /*! \brief This method splits the enabled subtrees of the tree rooted to
     * this node in such a way that the tree depth becomes \a depth.
     * If the initial tree depth is greater than \a depth then nothing is done.
     */
    void mince(const uint depth);
            
    //@}
            
    //@{
    //! \name
    
    /*! \brief Allows to assign one binary tree node to the other, this is done by copying
     * the nodes value and the sub-trees. The sub-trees of this node are deallocated.
     */
    BinaryTreeNode& operator=(const BinaryTreeNode & otherNode );
    
    /*! \brief Copy all the data (including the sub-nodes) from the node pointed by \a pOtherNode into the given node.
     * Note that here we will create copies of the sub nodes, and NOT just copy pointers to them!
     */
    void copy_from( const BinaryTreeNode * pOtherNode );
            
    /*! \brief This method recombined the sub tree nodes. Note that, the two leaf nodes with
     * the same parent are removed if they have the same value of isEnabled fields.
     */
    void recombine();
            
    /*! \brief Stores the binary tree in a form of two arrays, their structure is the same as needed for
     *   the BinaryTreeNode( const BooleanArray& , const BooleanArray&  ) constructor
     */
    void tree_to_binary_words( BinaryWord & tree, BinaryWord & leaves ) const;
            
    /*! \brief Stores the binary tree node as a string*/
    string node_to_string() const;
            
    /*! \brief Finds(creates) the leaf node defined by the \a path and marks it as enabled.
     * If some prefix of the \a path references an enabled node then nothing is done.
     */
    void add_enabled( const BinaryWord & path );

    /*! \brief This method adjoins the enabled nodes of \a subTree to this tree.
     * Note that, the position of the root node of \a subTree within this
     * tree is defined by path.
     */
    void add_enabled( const BinaryTreeNode * pOtherSubTree, const BinaryWord & path );

    /*! \brief This method merges \a pFromTreeRoot into \a pToTreeRoot.
     *  The enabled nodes of the former tree are added to the latter tree.
     *  If in the latter tree there is an enabled leaf node and in the former
     *  tree the corresponding node is not a leaf, then this node of the latter
     *  tree stays intact. If on the other hand we have a non-leaf node in
     *  \a pToTreeRoot and we are adding to it an enabled node of \a pFromTreeRoot
     *  then we just "substitute" the former one with the latter.
     *  NOTE: 1. This function is recursive. 2. No pointers are copied between 
     *  \a pToTreeRoot and \a pFromTreeRoot.
     */
    void add_enabled( BinaryTreeNode* pToTreeRoot, const BinaryTreeNode* pFromTreeRoot );
            
    /*! \brief Starting in the \a pNode node as at the root, this method counts 
     *  the number of enabled leaf nodes in the subtree
     *  rooted at pNode.
     */
    static size_t count_enabled_leaf_nodes( const BinaryTreeNode* pNode );
    
    /*! \brief Starting in the \a pRootTreeNode node as at the root, this method finds(creates)
     *  the leaf node defined by the \a path and marks it as enabled. If some prefix of the \a path
     *  references an enabled node then nothing is done.
     *  NOTE: This is a recursive method on the position (\a position) in the binary path (\a path).
     *  Therefore, the initial evaluate of this method should be done with \a position == 0;
     */
    static void add_enabled( BinaryTreeNode* pRootTreeNode, const BinaryWord& path, const uint position = 0 );
    
    /*! \brief Creates a binary tree of the height rootNodePath.size(), puts the subtree oldRootNode
     * into the node defined by the path \a rootNodePath, returns the root node of the extended tree.
     */
    static BinaryTreeNode * prepend_tree( const BinaryWord & rootNodePath, BinaryTreeNode * oldRootNode);
    
    /*! \brief This method restricts \a pThisNode to \a pOtherNode.
     * In essance we do the inplace AND on the tree node pThisNode.
     * Note that, this method is recursive.
     */
    static void restrict( BinaryTreeNode * pThisNode, const BinaryTreeNode * pOtherNode );
    
    /*! \brief This method removed enabled nodes of \a pOtherNode from \a pThisNode.
     * Note that, this method is recursive.
     */
    static void remove( BinaryTreeNode * pThisNode, const BinaryTreeNode * pOtherNode );
    
    /*! \brief checks if two trees overlap in a set-theory sence.
     * I.e. we assume that pRootNodeOne and pRootNodeTwo correspond to the same (virtual) root
     * and then we see if the enabled leaf node of one tree contain enabled leaf nodes of
     * another tree as their (virtual) children.
     */
    static bool overlap( const BinaryTreeNode * pRootNodeOne, const BinaryTreeNode * pRootNodeTwo );
    
    /*! \brief checks if the tree pRootNodeOne is a subset of the tree pRootNodeTwo, in a set-theory sence.
     * I.e. we assume that pRootNodeOne and pRootNodeTwo correspond to the same (virtual) root
     * and then we see if every enabled leaf node of pRootNodeOne is contained in the enabled leaf nodes of
     * pRootNodeTwo, or it is covered by the enabled leaf nodes of pRootNodeTwo. When we write, contained and
     * covered then we mean: is a subnode in the (virtual) tree and all it's subnodes in the (virtual) tree.
     */
    static bool subset( const BinaryTreeNode * pRootNodeOne, const BinaryTreeNode * pRootNodeTwo );
    
    //@}
};

/*! \brief An abstract cell of a grid paving. This class is the base of the GridCell - a regular cell on the Grid
 *  and the GridOpenCell - an open cell on a Grid. Here we only store common data and operations
 * 
 * This class contains the geometric data of a cell, as well as the combinatoric data needed to easily adjoin it to a tree.
 * It does not contain the tree structure, so cannot be used in cursors/iterators.
 * NOTE: The "Primary root cell" is the cell that the \a GridTreeSet can be rooted to
 * (the \a GridTreeSubset can be rooted to a non primary cell).
 */
class GridAbstractCell {
  protected:
    /*! \brief The intrinsic grid of the cell. Note that grids are internally passed by reference. */
    Grid _theGrid;

    /*! \brief The level of the primary root cell relative to the zero level. */
    uint _theHeight;

    /*! \brief The word describing the path in a binary tree. This path defines the cell in the Grid. */
    BinaryWord _theWord;

    /*! \brief The geometric box represented by the cell. */
    Box _theBox;
            
    /*! \brief The box is given as an explicit parameter. This should only be used
     *  by friend classes which have already computed the box.
     */
    GridAbstractCell(const Grid& theGrid, const uint theHeight, const BinaryWord& theWord, const Box& theBox);
    
    /*! \brief The copy constructor for the \a GridAbstractCell */
    GridAbstractCell(const GridAbstractCell& theGridCell);
            
    friend class GridTreeCursor;

    /*! \brief having \a theHeight the height of the primary cell, with \a leftBottomCorner and \a rightTopCorner
     * defining the primary cell corners for \a theHeight-1, we recompute \a leftBottomCorner and \a rightTopCorner
     * for the current level \a theHeight.
     */
    static inline void primary_cell_at_height( const uint theHeight, int & leftBottomCorner, int & rightTopCorner );

    /*! \brief This function allows to compare to cells it is used by the operator== and operator< methods of this class
     *  The value of \a comparator should be either \a COMPARE_EQUAL or \a COMPARE_LESS
     *  The function checks that both cells are on the same grid and then alignes their primary cells.
     *  The latter is done by extending the binary word of the cell with the lowest primary cell with
     *  the corresponding prefix. When the words are alligned, wi simply use the == and < operators of
     *  the \a BinaryWord class.
     *  NOTE: Since the function is only based on comparison of the Grids, the words definings the cells and the
     *  primary cell heights we can use it in sub-classes for comparing cells of the same types
     */
    static bool compare_abstract_grid_cells(const GridAbstractCell * pCellLeft, const GridAbstractCell &cellRight, const uint comparator );
    
  public:

    /*! \brief The constant for the grid cells comparison */
    static const uint COMPARE_EQUAL = 0;
    static const uint COMPARE_LESS = 1;
    
    /*! \brief The underlying grid. */
    const Grid& grid() const;
            
    /*! \brief The height of the primary cell to which this cell is rooted. */
    uint height() const;
            
    /*! \brief The depth in the grid at which the cell lies. */
    int depth() const;
            
    /*! \brief The word describing the path in a binary tree from the primary cell of height (this.height()) to this cell. */
    const BinaryWord& word() const;
            
    /*! \brief The geometric box represented by the cell. */
    const Box& box() const;

    /*! \brief The dimension of the cell. Needed by the BoxExpression interface. */
    dimension_type dimension() const;
    
    /*! \brief The upper and lower bound in the \a i<sup>th</sup> coordinate. */
    Interval operator[](dimension_type i) const;
    
    /*! \brief Allows to assign one GridCell to another */
    GridAbstractCell& operator=( const GridAbstractCell & otherCell );

    /*! \brief The primary-cell box on some grid is defined by the height above the zero level cells.
     *     The coordinates of the primary cells change uniformly. For example (3-D):
     *       Heigth | (left bottom) | (right top)
     *       ---------------------------------------
     *          0   | (0, 0, 0)     | (1, 1, 1)
     *          1   | (-1, -1, -1)  | (1, 1, 1)
     *          2   | (-1, -1, -1)  | (3, 3, 3)
     *          3   | (-5, -5, -5)  | (3, 3, 3)
     *          4   | (-5, -5, -5)  | (11, 11, 11)
     *     The cell size is always double the size of the previous-level cell and one
     *     of the corners stays the same as well. If Li, Ri are the left-, right-corned
     *     coordinates at livel i then we have the following recursive definition:
     *       L0 = 0, R0 = 0
     *       Li = 2*L(i-1) - R(i-1) (if i is odd)
     *       Li = L(i-1)            (if i is even)
     *       Ri = R(i-1)            (if i is odd)
     *       Ri = 2*R(i-1) - L(i-1) (if i is even)
     *
     *  WARNING: The returned box must be interptered relative to some grid.
     *  In other words, we assume some lattice and the relation to the original
     *  space is not taken into account.
     *
     * Here
     * \a theHeight defines the height of the primary cell above the zero level and
     * \a dimensions is the number of dimension in the considered space
     */
    static Vector<Interval> primary_cell_lattice_box( const uint theHeight, const dimension_type dimensions );
            
    /*! \brief Takes the lattice box \a theLatticeBox related to some (unknown) grid and computes the
     *   smallest primary cell on that grid that will contain this box. 
     *   The method returns the hight of that primary cell.
     */
    static uint smallest_enclosing_primary_cell_height( const Vector<Interval>& theLatticeBox );
    
    /*! /brief Computes the height of the primary cell that encloses
     *  (may be exactly) the given box on the given grid.
     */
    static uint smallest_enclosing_primary_cell_height( const Box & theBox, const Grid& theGrid );

    /*! \brief This method returns the path from the \a topPCellHeight to the \a bottomPCellHeight
     *  in the \a dimensions dimensional space. We assume that \a topPCellHeight >= \a bottomPCellHeight
     *  if not, then we return an empty binary word.
     */
    static BinaryWord primary_cell_path( const uint dimensions, const uint topPCellHeight, const uint bottomPCellHeight);

    /*! \brief Apply grid data \a theGrid to \a theLatticeBox in order to compute
     * the box dimensions in the original space
     */
    static Box lattice_box_to_space(const Vector<Interval> & theLatticeBox, const Grid& theGrid );
};

/*! \brief A cell of a grid paving. Note that this class represents the closed cell.
 *  It is uniquely defined by the path from the primary cell and is exactly on cell in a grid.
 *  This is different from the GridOpenCell which is open and the path from the root cell defines
 *  its lower left sub cell wheres the complete open cell is a union of 4 cells.
 * 
 * This class does not contain the tree structure, so cannot be used in cursors/iterators.
 */
class GridCell: public GridAbstractCell {
  protected:
    
    friend class GridTreeCursor;

    /*! \brief Tests if theCellOne is a subset of theCellTwo, the paths to the cells,
     *  aligned to some common primary cell, which height will be referenced by
     *  \a pPrimaryCellHeight, are returned as \a pathPrefixOne and \a pathPrefixTwo
     */
    friend bool subset( const GridCell& theCellOne, const GridCell& theCellTwo, BinaryWord * pPathPrefixOne,
                            BinaryWord * pPathPrefixTwo, uint * pPrimaryCellHeight );

  public:
    /*! \brief Default constructor. Needed for some containers and iterators. */
    GridCell();

    /*! \brief Construct a cell based on \a theGrid, the primary cell of level \a theHeight,
     * and the path to the SubPaving's root cell which is accessible from the primary cell
     * via the path \a _theWord. Note that, \a theHeight is the height relative to the Grid,
     * but not to the original space!
     */
    GridCell(const Grid& theGrid, const uint theHeight, const BinaryWord& theWord);
	
	/*! \brief Allows to split the given cell into two sub-cells. When isRight == true
	 * then we return the right sub-cell, otherwise the left one */
	GridCell split(bool isRight) const;

    /*! \brief The equality operator. */
    bool operator==(const GridCell& otherCell) const;
    
    /*! \brief Allows to assign one GridCell to another */
    GridCell& operator=( const GridCell & otherCell );
    
    /*! \brief A total order on cells on the same grid, by height and word prefix. */
    bool operator<(const GridCell& otherCell) const;
    
    /*! \brief Allows to convert the given GridCell into an open grid cell (GridOpenCell)*/
    GridOpenCell interior() const;
    
    /*! \brief this method computes the box corresponding to this cell in the grid lattice.*/
    static Vector<Interval> compute_lattice_box( const uint dimensions, const uint theHeight, const BinaryWord& theWord );
    
    /*! \brief this method computes the box in the original space based on the \a theGrid,
     *  and a cell which is obtained by traversing the path given by \a _theWord from the
     * primary cell located at the heigth \a theHeight above the zero level cells
     * (relative to the grid).
     */
    static Box compute_box(const Grid& theGrid, const uint theHeight, const BinaryWord& theWord);

    /*! /brief Creates an over approximation for the \a theBox on \a theGrid. \a theBox
     * is in the original space coordinates. We compute the over approximation as the
     * smallest primary cell on the Grid, such that it contains \a theBox (after it's
     * mapping on \a theGrid )
     */
    static GridCell smallest_enclosing_primary_cell( const Box & theBox, const Grid& theGrid );

    /*! \brief This function allows to compare to cells it is used by the operator== and operator< methods of this class
     *  The value of \a comparator should be either \a COMPARE_EQUAL or \a COMPARE_LESS
     *  The function checks that both cells are on the same grid and then alignes their primary cells.
     *  The latter is done by extending the binary word of the cell with the lowest primary cell with
     *  the corresponding prefix. When the words are alligned, wi simply use the == and < operators of
     *  the \a BinaryWord class.
     */
    static bool compare_grid_cells(const GridCell * pCellLeft, const GridCell &cellRight, const uint comparator );
    
    /*! \brief Computes the neighboring cell to the given cell (defined by \a theGrid, \a theHeight, \a theWord)
     *  in the given dimension \a dim. The resulting cell will the of the same size as the given one.
     */
    static GridCell neighboringCell( const Grid& theGrid, const uint theHeight, const BinaryWord& theWord, const uint dim );
};

/*! \brief An open cell of a grid paving. This cell is open and the path from the primary cell
 *  defines its lower left sub cell wheres the complete open cell is a union of 4 cells.
 * 
 * This class does not contain the tree structure, so cannot be used in cursors/iterators.
 * NOTE: The open cell is defined by its base cell which is a simple GridCell corresponding
 * to the left-most quadrant cell of the open cell. In this respect, height() and word() return
 * values defining the base cell and box() returns the entire box related to the open cell.
 * This box should be treated as an open set, i.e. its borders should be excluded.
 */
class GridOpenCell: public GridAbstractCell {
  protected:
    
    friend class GridCell;
    
    /*! \brief The box is given as an explicit parameter. This should only be used
     *  by friend classes which have already computed the box.
     */
    GridOpenCell(const Grid& theGrid, const uint theHeight, const BinaryWord& theWord, const Box& theBox);
    
    /*! \brief This method allows to enumerate the GridCells that are located in the positive
     *  axis direction from the base cell. The base cell is defined by \a theHeight and
     *  \a theBaseCellWord. When called \a cellPosition should be an empty word. This is a
     *  technical parameter used in the method's recursivce calls. \a theResultSet is the vector
     *  to which all the neighboring cell will be added, i.e. it is a return parameter.
     *  NOTE: The cell defined by \a theHeight, \a theBaseCellWord will also be in \a theResultSet.
     */  
    void neighboring_cells( const uint theHeight, const BinaryWord& theBaseCellWord,
                            BinaryWord& cellPosition, GridTreeSet& theResultSet ) const;
    
    /*! \brief This method allows to compute the neighboring (to the right) cell of
     *  the base cell given by \a theHeight and \a theBaseCellWord. Here \a cellPosition
     *  is the position of the neighboring cell with resect to the base cell in the
     *  theGrid.dimensions() dimensional space. Note that if \a cellPosition consists
     *  of theGrid.dimensions() zeroes then we adjoin the base cell itself.
     */
    static GridCell neighboring_cell( const Grid& theGrid, const uint theHeight,
                                      const BinaryWord& theBaseCellWord, BinaryWord& cellPosition );
    
    /*! \brief This method allows to find the smallest open cell that contains \a theBox.
     *  The search is started from \a theOpenCell that is an open cell covering \a theBox.
     *  Note: This method is recursive and it assumes that \a theOpenCell
     *  and \a theBox are on the same Grid! The latter is not checked, but only assumed.
     *  In the open cell based on \a theOpenCell does not cover \a theBox this method returns NULL.
     */
    static GridOpenCell * smallest_open_subcell( const GridOpenCell &theOpenCell, const Box & theBox );
    
    /*! \brief Takes the cell \a theCell of the set where it belongs, \a theSet, and see if it has
     *  neighbors in this set. If yes then it creates the open cell lying inside theSet such that
     *  it covers the borders. All such open cells are added to the vector \a result. \a cellPosition
     *  is a technical parameter that has to be set to an empty word. Also, this method add the interior
     *  of theCell to the vector \a result.
     */
    static void cover_cell_and_borders( const GridCell& theCell, const GridTreeSet& theSet,
                                        BinaryWord& cellPosition, std::vector<GridOpenCell>& result );
    
  public:
    /*! \brief Default constructor. Needed for some containers and iterators. */
    GridOpenCell();

    /*! \brief Construct a cell based on \a theGrid, the primary cell of level \a theHeight,
     * and the path to the SubPaving's root cell which is accessible from the primary cell
     * via the path \a _theWord. Note that, \a theHeight is the height relative to the Grid,
     * but not to the original space!
     */
    GridOpenCell(const Grid& theGrid, const uint theHeight, const BinaryWord& theWord);
	
	/*! \brief Allows to split the given cell into two sub-cells. When isRight == true
	 * then we return the right sub-cell, if false then the left one, otherwise the middle one */
	GridOpenCell split(tribool isRight) const;

    /*! \brief The equality operator. */
    bool operator==(const GridOpenCell& otherCell) const;
    
    /*! \brief Allows to assign one GridCell to another */
    GridOpenCell& operator=( const GridOpenCell & otherCell );
    
    /*! \brief A total order on cells on the same grid, by height and word prefix. */
    bool operator<(const GridOpenCell& otherCell) const;
    
    /*! \brief Computes all the cells that constitute the GridOpenCell in the form of the GridTreeSet.*/
    GridTreeSet closure() const;
    
    /*! \brief Computes the intersection of two GridOpenCell as a list of open cells whoes union
     * gives the intersection. This is done becase the intersection can no always be represented
     * as one open cell. If the intersection is empty then it returns an empty vector.
     * WARNING: This operation is expensive because we do not always get a single cell as a result.
     * Moreover, currently intersecting an open cell with itself or a subcell sharing a common border
     * with the cell will result in more that one open cell. Still the operation is robust, it always
     * results in a set of open cells whoes union is the exact intersection of the given open cells.
     */
    static std::vector<GridOpenCell> intersection( const GridOpenCell & theLeftOpenCell, const GridOpenCell & theRightOpenCell );
    
    /*! \brief Tests if the two open cells overlap */
    static bool overlap( const GridOpenCell & theLeftOpenCell, const GridOpenCell & theRightOpenCell );
    
    /*! \brief this method computes the box in the original space based on the \a theGrid.
     * This box should be treated as an open set. I.e. the borders of the box must be excluded.
     * Note tha, \a theHeight and \a theWord define the left bottom cell (base cell) of the open cell.
     */
    static Box compute_box(const Grid& theGrid, const uint theHeight, const BinaryWord& theWord);

    /*! \brief This function allows to compare to cells it is used by the operator== and operator< methods of this class
     *  The value of \a comparator should be either \a COMPARE_EQUAL or \a COMPARE_LESS
     *  The function checks that both cells are on the same grid and then alignes their primary cells.
     *  The latter is done by extending the binary word of the cell with the lowest primary cell with
     *  the corresponding prefix. When the words are alligned, wi simply use the == and < operators of
     *  the \a BinaryWord class.
     */
    static bool compare_grid_cells(const GridOpenCell * pCellLeft, const GridOpenCell &cellRight, const uint comparator );
    
    /*! \brief This method computes the smallest open cell that contains the given box \a theBox
     *  (in the original space), taking into account the given grid \a theGrid.
     */
    static GridOpenCell outer_approximation( const Box & theBox, const Grid& theGrid );
};

/*! \brief This class represents a subpaving of a paving. Note that, the subtree enclosed into
 * this class is just a pointer to the node in the tree of some paving. This class is not
 * responsible for deallocation of that original tree.
 */
class GridTreeSubset {
  protected:
            
    friend class GridTreeCursor;

    friend GridTreeSet outer_approximation( const CompactSetInterface& theSet, const Grid& theGrid, const uint depth );
    friend GridTreeSet inner_approximation( const OpenSetInterface& theSet, const Grid& theGrid, const uint depth );
            
    /*! \brief The pointer to the root node of the subpaving tree.
     * Note that, this is not necessarily the root node of the corresponding paving tree.
     */
    BinaryTreeNode * _pRootTreeNode;
            
    /*! \brief The paving cell corresponding to the root node of the SubPaving.*/
    GridCell _theGridCell;

    /*! \brief this function takes the interval width and computes how many binary subdivisions
     * one has to make in order to have sub-intervals of the width <= \a theMaxWidth
     */
    uint compute_number_subdiv( Float theWidth, const Float theMaxWidth) const;
    
    /*! \brief This method checks whether the set defined by \a pCurrentNode is a superset
     *  of \a theBox, in case when it is known that the cell corresponding to the root of
     *  pCurrentNode [and defined by \a theGrid, the primary cell (\a theHeight) to which
     *  this tree is (virtually) rooted via the path theWord] encloses \a theBox.
     *  This is a recursive procedure and it returns true only if there are no disabled
     *  cells in \a pCurrentNode that intersect with theBox.
     */
    static tribool covers( const BinaryTreeNode* pCurrentNode, const Grid& theGrid,
                                const uint theHeight, BinaryWord &theWord, const Box& theBox );
    
    /*! \brief This method checks whether the set defined by \a pCurrentNode is a subset
     *  of \a theBox. The set of \a pCurrentNode is defined by \a theGrid, the primary
     *  cell (\a theHeight) to which this tree is (virtually) rooted via the path theWord.
     *  This is a recursive procedure and it returns true only if all enabled sub-cells of
     *  \a pCurrentNode are sub-sets of \a theBox.
     */
    static tribool subset( const BinaryTreeNode* pCurrentNode, const Grid& theGrid,
                           const uint theHeight, BinaryWord &theWord, const Box& theBox );

    /*! \brief This method checks whether \a theBox is disjoint from the set defined by
     *  \a pCurrentNode, \a theGrid, the primary cell (\a theHeight) to which this
     *  tree is (virtually) rooted via the path theWord. This is done using the recursive
     *  procedure by checking the cells of the tree that intersect with the box and going
     *  down to the leaves. When reaching a leaf node that is enabled we conclude that we
     *  have an intersection. If there are no such nodes then there is no intersection,
     *  and the sets are disjoint.
     */
    static tribool disjoint( const BinaryTreeNode* pCurrentNode, const Grid& theGrid,
                             const uint theHeight, BinaryWord &theWord, const Box& theBox );
    
    /*! \brief This method checks whether \a theBox overlaps the set defined by
     *  \a pCurrentNode, \a theGrid, the primary cell (\a theHeight) to which this
     *  tree is (virtually) rooted via the path theWord. This is done using the recursive
     *  procedure by checking the cells of the tree that overlap the box and going
     *  down to the leaves. When reaching a leaf node that is enabled we conclude that we
     *  have an overlap. If there are no such nodes then there is no overlap. Note that this 
     *  method cannot be implemented using disjoint since disjoint tests if the closures 
     *  are disjoint, and overlaps tests if the interiors are not disjoint.
     */
    static tribool overlaps( const BinaryTreeNode* pCurrentNode, const Grid& theGrid,
                             const uint theHeight, BinaryWord &theWord, const Box& theBox );
    
  public:
    
    /*! \brief A short name for the constant iterator */
    typedef GridTreeConstIterator const_iterator;
            
    //@{
    //! \name Constructors
            
    /*! \brief The new root node can only be constructed from the existing tree node.
     * Here, \a pRootTreeNode is not copied, we simply store its pointer inside this class.
     * Note that, \a pRootTreeNode should correspond to the sub-paving root node. Thus,
     * \a theHeight defines the height of the primary root cell of the GridTreeSet
     * (Remember that every GridTreeSubset is just a reference to a subtree of a GridTreeSet).
     * \a theWord defines the path to the \a pRootTreeNode node from the primary root cell
     * of the corresponding GridTreeSet.
     */
    GridTreeSubset( const Grid& theGrid, const uint theHeight, const BinaryWord& theWord, BinaryTreeNode * pRootTreeNode );
    
    /*! \brief A copy constructor that only copies the pointer to the root of the binary tree and the cell */
    GridTreeSubset( const GridTreeSubset &otherSubset);
            
    //@}

    /*! Virtual destructor. The destructor needs to be virtual since GridTreeSet is a subclass 
     *  with different memory management. */
    virtual ~GridTreeSubset();
        
    //@{
    //! \name Properties
             
    /*! \brief True if the set is empty. */
    bool empty() const; 

    /*! \brief The number of activated cells in the set. */
    size_t size() const; 

    /*! \brief The dimension of the set. */
    dimension_type dimension() const;

    /*! \brief Returns a constant reference to the underlying grid. */
    const Grid& grid() const;
            
    /*! \brief Returns the const pointer to the root BinaryTreeNode of the SubPaving*/
    const BinaryTreeNode * binary_tree() const;
            
    /*! Recalculate the depth of the tree rooted at \a _pRootTreeNode */
    uint depth() const;

    /*! The measure (area, volume) of the set in Euclidean space. */
    double measure() const;

    /*! \brief Returns the \a GridCell corresponding to the ROOT NODE of this \a GridTreeSubset
     * WARNING: It is NOT the primary cell of the paving! 
     */
    GridCell cell() const;
            
    /*! \brief Computes a bounding box for a grid set. */
    Box bounding_box() const;
            
    /*! \brief Allows to test if the two subpavings are "equal". The method returns true if
     * the grida are equal and the binary trees are equal. Note that, only in case both 
     * GridTreeSubset objects are recombines, this method is guaranteed to tell you that
     * the two GridTreeSubset represent equal sets.
     */
    bool operator==(const GridTreeSubset& anotherGridTreeSubset) const;
            
    //@}

    //@{
    //! \name Subdivisions
            
    /*! \brief Subdivides the tree up to the depth specified by the parameter.
     * Note that, we start from the root of the sub-paving and thus the subdivision
     * is done relative to it but not to the root of the original paving.
     */
    void mince( const uint theNewDepth );

    /*! \brief Subdivide the paving until the smallest depth such that the leaf
     * cells size is <= \a theMaxCellWidth. Note that, the disabled cells are
     * not subdivided.
     */
    void subdivide( Float theMaxCellWidth );
            
    /*! \brief Recombines the subdivisions, for instance if all subcells of a cell are
     * enabled/disabled then they are put together.
     */
    void recombine();
            
    //@}
    
    //@{
    //! \name Geometric Predicates

    /*! \brief Tests if a cell is a subset of a set. */
    friend bool subset( const GridCell& theCell, const GridTreeSubset& theSet );
    
    /*! \brief Tests if a cell overlaps (as an open set) a paving set. */
    friend bool overlap( const GridCell& theCell, const GridTreeSubset& theSet );

    /*! \brief Tests if a grid set \a theSet1 is a subset of \a theSet2. */
    friend bool subset( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 );
    
    /*! \brief Tests if two grid paving sets overlap (i.e. intersect as open sets.)
     *  If at least one of the GridTreeSubsets represents an empty set, then the result is false.
     */
    friend bool overlap( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 );

    /* PIETER: We may want the two predicates below to return tribool, reflecting the
     * case that the box is an approximation. Unfortunately, this is very difficult
     * for subset.
     */

    /*! \brief Tests if a grid set is a subset of a box. */
    tribool subset( const Box& theBox ) const;
    
    /*! \brief Tests if a grid set is a superset of a box. */
    tribool superset( const Box& theBox ) const;
    
    /*! \brief Tests if (the closure of) a grid set is disjoint from a box. */
    tribool disjoint( const Box& theBox  ) const;
    
    /*! \brief Tests if a grid set overlaps a box. */
    tribool overlaps( const Box& theBox ) const;
            
    //@}
            
    //@{
    //! \name Iterators

    /*! \brief A constant iterator through the enabled leaf nodes of the subpaving. */
    const_iterator begin() const;
            
    /*! \brief A constant iterator to the end of the enabled leaf nodes of the subpaving. */
    const_iterator end() const;

    //@}

    //@{
    //! \name Conversions 

    /*! \brief An assignment operator that only copies the pointer to the root of the binary tree and the cell */
    GridTreeSubset& operator=( const GridTreeSubset &otherSubset);
            
    /*! \brief Convert to a list of ordinary boxes, unrelated to the grid. */
    operator ListSet<Box>() const;
            
    //@}
};

/*! \brief The GridTreeSet class that represents a set of cells with mixed integer and dyadic coordinates.
 * The cells can be enabled or disabled (on/off), indicating whether they belong to the paving or not.
 * It is possible to have cells that are neither on nor off, indicating that they have enabled and
 * disabled sub cells.
 * 
 * A trivial cell (level 0) of the paving corresponds to a unit hypercube, of the n-dimensional state
 * space. Subdivision of any 0-level paving is based on dyadic numbers and thus is a binary-style
 * partitioning of the cell. All cells with purely integer coordinates are enclosed in a (virtual)
 * binary tree. We illustrate this for a two dimensional case:
 *  1. Take the unit cell [0, 0]*[1, 1],
 *  2. Cells [-1, 0]*[0, 1] and [0, 0]*[1, 1] are rooted to [-1, 0]*[1, 1],
 *  3. Cells [-1, -1]*[0, 0] and [0, -1]*[1, 0] are rooted to [-1, -1]*[1, 0],
 *  4. Cells [-1, -1]*[1, 0] and [-1, 0]*[1, 1] are rooted to [-1, -1]*[1, 1].
 */
class GridTreeSet : public GridTreeSubset {        
  protected:

    friend class GridTreeCursor;

    /*! \brief This method takes the height of the primary cell
     *  \a otherPavingPCellHeight and if it is:
     *    (a) higher then for this paving, pre-pends \a _pRootTreeNode.
     *    (b) lower then for this paving, locates it in this paging's tree.
     *    (c) equal to the hight of this paving's primary cell, does nothing.
     *  This method returns the node corresponding to the primary cell of
     *  height \a otherPavingPCellHeight in the (updated) paving.
     *  If \a stop_on_enabled is set to true then for the case (b) if we meet
     *  a leaf node on the path from the primary node of the paving to the primary
     *  node of the cell, we stop locating the node corresponding to the primary
     *  cell of height \a otherPavingPCellHeight and set \a has_stopped to true.
     *  If \a stop_on_disabled is set to true then for the case (b) if we meet
     *  a leaf node on the path from the primary node of the paving to the primary
     *  node of the cell, we stop locating the node corresponding to the primary
     *  cell of height \a otherPavingPCellHeight and set \a has_stopped to true.
     */
    BinaryTreeNode* align_with_cell( const uint otherPavingPCellHeight, const bool stop_on_enabled, const bool stop_on_disabled, bool & has_stopped );

    /*! \brief This method adjoins the outer approximation of \a theSet (computed on the fly) to this paving.
     *  We use the primary cell (enclosed in this paving) of height \a primary_cell_hight and represented 
     *  by the paving's binary node \a pBinaryTreeNode. When adding the outer approximation, we compute it
     *  up to the level of accuracy given by \a max_mince_depth. This parameter defines, how many subdivisions
     *  of the binary tree we should make to get the proper cells for outer approximating \a theSet.
     *  This method is recursive, the parameter \a pPath defines the path to the current node pBinaryTreeNode
     *  from the root node in recursive calls, thus the initial evaluate for this method must be done with an empty word.
     */
    static void _adjoin_outer_approximation( const Grid & theGrid, BinaryTreeNode * pBinaryTreeNode, const uint primary_cell_height,
                                             const uint max_mince_depth, const CompactSetInterface& theSet, BinaryWord * pPath );

    /*! \brief This method adjoins the inner approximation of \a theSet (computed on the fly) to this paving.
     *  We use the primary cell (enclosed in this paving) of height \a primary_cell_height and represented 
     *  by the paving's binary node \a pBinaryTreeNode. When adding the inner approximation, we compute it
     *  up to the level of accuracy given by \a max_mince_depth. This parameter defines, how many subdivisions
     *  of the binary tree we should make to get the proper cells for inner approximating \a theSet.
     *  This method is recursive, the parameter \a pPath defines the path to the current node pBinaryTreeNode
     *  from the root node in recursive calls, thus the initial evaluate for this method must be done with an empty word.
     */
    static void _adjoin_inner_approximation( const Grid & theGrid, BinaryTreeNode * pBinaryTreeNode, const uint primary_cell_height,
                                             const uint max_mince_depth, const OpenSetInterface& theSet, BinaryWord * pPath );

    /*! \brief This method adjoins the lower approximation of \a theSet (computed on the fly) to this paving.
     *  We use the primary cell (enclosed in this paving) of height \a primary_cell_hight and represented 
     *  by the paving's binary node \a pBinaryTreeNode. When adding the lower approximation, we compute it
     *  up to the level of accuracy given by \a max_mince_depth. This parameter defines, how many subdivisions
     *  of the binary tree we should make to get the proper cells for lower approximating \a theSet.
     *  This method is recursive, the parameter \a pPath defines the path to the current node pBinaryTreeNode
     *  from the root node in recursive calls, thus the initial evaluate for this method must be done with an empty word.
     *  The approximation method does not recombine cells, as knowing that both children intersect a set is more
     *  information than knowing that the parent does.
     */
    static void _adjoin_lower_approximation( const Grid & theGrid, BinaryTreeNode * pBinaryTreeNode, const uint primary_cell_height,
                                             const uint max_mince_depth, const OvertSetInterface& theSet, BinaryWord * pPath );

    /*! \brief This method adjoins the lower approximation of \a theSet (computed on the fly) to this paving.
     *  It is specialised for open sets, for which we have the superset() operator. If a set is a superset of
     *  a cell, then we know it overlaps the cell and all its children.
     */
    static void _adjoin_lower_approximation( const Grid & theGrid, BinaryTreeNode * pBinaryTreeNode, const uint primary_cell_height,
                                             const uint max_mince_depth, const OpenSetInterface& theSet, BinaryWord * pPath );

    /*! \brief This method is uset to do restriction of this set to the set given by
     *  \a theOtherSubPaving Note that, here we require that the height of the primary
     *  root cell of this set is >= the height of the primary root cell of \a theOtherSubPaving.
     */
    void restrict_to_lower( const GridTreeSubset& theOtherSubPaving );
            
    /*! \brief This method is uset to remove \a theOtherSubPaving from this set.
     *  Note that, here we require that the height of the primary root cell of
     *  this set is >= the height of the primary root cell of \a theOtherSubPaving.
     */
    void remove_from_lower( const GridTreeSubset& theOtherSubPaving );

    /*! \brief This method changes the primary cell of this GridTreeSet.
     *  We only can increase the height of the primary cell, this is why
     *  if toPCellHeight <= this->cell().height(), then nothing is done.
     */
    void up_to_primary_cell( const uint toPCellHeight );

  public:
    //@{
    //! \name Constructors
            
    /*! \brief Create a %GridTreeSet based on zero dimensions. 
     *  This constructor is needed to use the Boost Serialization library.
     */
    GridTreeSet( );
            
    /*! \brief The new root node can only be constructed from the existing tree node.
     *  Here, \a pRootTreeNode is not copied, we simply store its pointer inside this class.
     *  Note that, \a pRootTreeNode should correspond to the root node. \a theHeight defines
     *  the height of the primary root cell corresponding to the \a pRootTreeNode node.
     */
    GridTreeSet( const Grid& theGrid, const uint theHeight, BinaryTreeNode * pRootTreeNode );
            
    /*! \brief Construct a grid tree set from a single cell.
     */
    GridTreeSet( const GridCell & theGridCell );

    /*! \brief The copy constructor that actually copies all the data,
     *  including the paving tree. I.e. the new copy of the tree is created.
     */
    GridTreeSet( const GridTreeSet & theGridTreeSet );

    /*! A simple constructor that creates the [0, 1]*...*[0, 1] cell in the
     *  \a theDimension - dimensional space. Here we assume that we have a non scaling
     *  grid with no shift of the coordinates. I.e. Grid._data._origin = {0, ..., 0}
     *  and Grid._data._lengths = {1, ..., 1}. If enable == true then the cell is enabled
     */
    explicit GridTreeSet( const uint theDimension, const bool enable = false );

    /*! \brief Construct an empty tree. The \a theBoundingBox is used to define the lattice
     *  block (in theGrid) that will correspond to the root of the paving tree \a pRootTreeNode.
     */
    explicit GridTreeSet( const Grid& theGrid, const bool enable = false  );

    /*! \brief Construct an empty tree. The \a theLatticeBox is used to define the lattice
     *  block (in theGrid) that will correspond to the root of the paving tree \a pRootTreeNode.
     */
    explicit GridTreeSet( const Grid& theGrid, const Box & theLatticeBox );
            
    /*! \brief Construct the paving based on the block's coordinates, defined by: \a theLeftLowerPoint
     *  and \a theRightUpperPoint. These are the coordinates in the lattice defined by theGrid.
     *  The primary cell, enclosing the given block of cells is computed automatically and the
     *  binary tree is rooted to that cell. The \a theEnabledCells array defines the enabled/disabled
     *  0-level cells of the block (lexicographic order).
     */
    explicit GridTreeSet( const Grid& theGrid, const IndexArray theLeftLowerPoint,
                          const IndexArray theRightUpperPoint, const BooleanArray& theEnabledCells );

    /*! \brief Creates a new paving from the user data. \a theTree is an array representation of the binary
     * tree structure, \a theEnabledCells tells whether a node is or is not a leaf, \a theHeight gives the
     * height of the primary cell which is assumed to correspond to the root node of \a theTree.
     */
    explicit GridTreeSet( const Grid& theGrid, uint theHeight, const BooleanArray& theTree, const BooleanArray& theEnabledCells );

    //@}

    //@{
    //! \name Cloning/Copying/Assignment
            
    /*! \brief The copy assignment operator, which copies all the data 
     *  including the paving tree if necessary.
     */
    GridTreeSet& operator=( const GridTreeSet & theGridTreeSet );

    /*! \brief Return a new dynamically-allocated copy of the %GridTreeSet.
     *  In this case, all the data is copied.
     */
    GridTreeSet* clone() const;

    //@}
            
    /*! \brief Destructor, removes all the dynamically allocated data, any */
    /* GridTreeSubset referencing this %GridTreeSet becomes invalid. */
    virtual ~GridTreeSet();
            
       
            
    //@{
    //! \name Geometric Operations
    /* PIETER: You may prefer to make inplace operations may return
     *a reference to *this, allowing chaining of operations.
     */

    /*! \brief Clears the set (makes empty set on same grid). */
    void clear( );
            
    /*! \brief Adjoin (make inplace union with) a single cell. */
    void adjoin( const GridCell& theCell );
            
    /*! \brief Remove a single cell. */
    void remove( const GridCell& theCell );

    /*! \brief Adjoin (make inplace union with) another grid paving set. */
    void adjoin( const GridTreeSubset& theOtherSubPaving );
            
    /*! \brief Restrict to (make inplace intersection with) another grid paving set. */
    void restrict( const GridTreeSubset& theOtherSubPaving );
            
    /*! \brief Remove cells in another grid paving set. */
    void remove( const GridTreeSubset& theOtherSubPaving );
            
    /*! \brief Restrict to cells rooted to the primary cell with the height (at most) \a theHeight. */
    void restrict_to_height( const uint theHeight );
 
    //@}
            
    //@{
    //! \name Geometric Operations
           
    /*! \brief Join (make union of) two grid paving sets. */
    friend GridTreeSet join( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 );
    
    /*! \brief The intersection of two grid paving sets. Points only lying on the
     *  intersection of the boundaries of the two sets are not included in the result.
     */
    friend GridTreeSet intersection( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 );
            
    /*! \brief The difference of two grid paving sets. (Results in theSet1 minus theSet2) */
    friend GridTreeSet difference( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 );
            
    //@}

    //@{
    //! \name Geometric Approximation

    /*! \brief Adjoin an over approximation to box, computing to the given depth.
     *  \pre The box must have nonempty interior.
     */
    void adjoin_over_approximation( const Box& theBox, const uint depth );
    
    /*! \brief Adjoin an outer approximation to a given set, computing to the given depth.
     *  This method computes an outer approximation for the set \a theSet on the grid \a theGrid.
     *  Note that, the depth is the total number of subdivisions (in all dimensions) of the unit
     *  cell of the grid. This method does the followig:
     * 1. Computes the smallest Primary cell enclosing \a theSet
     * 2. Allocates the paving for this cell
     * 3. Minces the paving to the level: depth + \<the primary cell height\>
     * 4. Iterates through the enabled leaf nodes of the paving (all the nodes are initially enabled)
     * 5. Disables the cells that are disjoint with the \a theSet
     */
    void adjoin_outer_approximation( const CompactSetInterface& theSet, const uint depth );
            
    /*! \brief Adjoin a lower approximation to a given set, computing to the given height and depth. 
     *   A lower approximation comprises all cells intersecting a given set.
     */
    void adjoin_lower_approximation( const OvertSetInterface& theSet, const uint height, const uint depth );

    /*! \brief Adjoin a lower approximation to a given set restricted to the given bounding box, computing to the given depth. 
     *   A lower approximation comprises all cells intersecting a given set.
     */
    void adjoin_lower_approximation( const OvertSetInterface& theSet, const Box& bounding_box, const uint depth );

    /*! \brief Adjoin a lower approximation to a given set, computing to the given depth. 
     *   A lower approximation comprises all cells intersecting a given set.
     */
    void adjoin_lower_approximation( const LocatedSetInterface& theSet, const uint depth );

    /*! \brief Adjoin an inner approximation to a given set, computing to the given height and depth. 
     *   An inner approximation comprises all cells that are sub-cells of the given set.
     */
    void adjoin_inner_approximation( const OpenSetInterface& theSet, const uint height, const uint depth );

    /*! \brief Adjoin an inner approximation to a given set restricted to the given bounding box, computing to the given depth. 
     *   An inner approximation comprises all cells that are sub-cells of the given set.
     */
    void adjoin_inner_approximation( const OpenSetInterface& theSet, const Box& bounding_box, const uint depth );

    //@}

    //@{
    //! \name Input/output routines.
    //@}

};

/*! \brief This class represents a cursor/iterator that can be used to traverse a subtree. */
class GridTreeCursor {
  private:
    //The size with which the stack size will be incremented
    static const uint STACK_SIZE_INCREMENT = 256;
            
    //The index of the top stack element (within the stack array)
    //WARNING: the initial value should be -1 to indicate that
    //there are no elements in the stack
    int _currentStackIndex;
            
    /* The nodes traversed to the current location. */
    array<BinaryTreeNode*> _theStack;
            
    /* The subpaving to cursor on */
    const GridTreeSubset * _pSubPaving;
            
    GridCell _theCurrentGridCell;

    /*! \brief Push the node into the atack */
    void push( BinaryTreeNode* pLatestNode );
            
    /*! \brief Pop the node from the atack, return NULL is the stack is empty */
    BinaryTreeNode* pop( );
            
    /*! \brief Check is the stack contains just one element, i.e. we are at the root */
    bool is_at_the_root() const;
            
    /*! \brief this method is supposed to update the _theCurrentGridCell value, when the cursos moves
     * if left_or_right == false then we go left, if left_or_right == true then right
     * if left_or_right == indeterminate then we are going one level up.
     */
    void updateTheCurrentGridCell( tribool left_or_right );

    friend std::ostream& operator<<(std::ostream& os, const GridTreeCursor& theGridTreeCursor);
            
  public:
    /*! \brief Default constructor constructs an invalid cursor. */
    GridTreeCursor();
            
    /*! \brief The constructor that accepts the subpaving to cursor on */
    GridTreeCursor(const GridTreeSubset * pSubPaving);
            
    /*! \brief The simple copy constructor, this constructor copies the pointer to
     *  GridTreeSubset and thus the internal stack information remains valid.
     */
    GridTreeCursor(const GridTreeCursor & otherCursor);
    
    /*! \brief This is an assignment operator that acts similar to the copy
     *  constructor. Here we copy the pointer to underlying GridTreeSubset.
     */
    GridTreeCursor & operator=(const GridTreeCursor & otherCursor);

    /*! This destructor does not deallocate the enclosed sub paving. */
    ~GridTreeCursor();

    /*! \brief Test if the current node is enabled. */
    bool is_enabled() const;
            
    /*! \brief Test if the current node is disabled*/
    bool is_disabled() const;
            
    /*! \brief Test if the current node is a leaf. */
    bool is_leaf() const;
            
    /*! \brief Test if the current node is the root of the subtree. 
     *  Returns true if the stack has only one element. */
    bool is_root() const;

    /*! \brief Returns true if the current node is the left child of the parent node */
    bool is_left_child() const;

    /*! \brief Returns true if the current node is the right child of the parent node */
    bool is_right_child() const;

    //@{
    //! \name Leaf Operations
            
    /*! \brief Markes the leaf node as enabled, otherwise they through \a NotALeafNodeEsception */
    void set_enabled() const;
            
    /*! \brief Markes the leaf node as disabled, otherwise they through \a NotALeafNodeEsception */
    void set_disabled() const;
            
    //@}

    /*! \brief Move to the parent node. Throws an NotAllowedMoveException if the current node is the root. 
     *   Returns a reference to itself. */
    GridTreeCursor& move_up();
            
    /*! \brief Move to the left child node. Throws an NotAllowedMoveException if the current node is a leaf. */
    GridTreeCursor& move_left();
            
    /*! \brief Move to the right child node. Throws an NotAllowedMoveException if the current node is a leaf. */
    GridTreeCursor& move_right();
            
    /*! \brief Move to a child node. Throws an NotAllowedMoveException if the current node is a leaf.
     * if left_or_right == false then we go left, if left_or_right == true then right
     */
    GridTreeCursor& move(bool left_or_right);

    /*! \brief Convert to a GridCell. */
    const GridCell& cell() const;
            
    /*! \brief Allows to test if the two cursors are equal, this is determined by
     * the fact that they point to the same binary-tree node.
     * NOTE: if _currentStackIndex < 0 for at least one of the cursors then the
     * result is always false.
     */
    bool operator==(const GridTreeCursor& anotherGridTreeCursor) const;

    /*! \brief The dereferencing operator which returns a
     * reference to the GridTreeSubset for the current node.
     */
    GridTreeSubset operator*();
            
    /*! \brief The dereferencing operator which returns a constant
     * reference to the GridTreeSubset for the current node.
     */
    const GridTreeSubset operator*() const;
};
    
/*! \brief This class allows to iterate through the enabled leaf nodes of GridTreeSubset.
 * The return objects for this iterator are constant GridCells.
 */
class GridTreeConstIterator : public boost::iterator_facade< GridTreeConstIterator, GridCell const, boost::forward_traversal_tag > {
  private:
    /*! \brief When set to true indicates that this is the "end iterator" */
    bool _is_in_end_state;
            
    /*! \brief the cursor object that is created based on the given sub paving*/
    GridTreeCursor _pGridTreeCursor;
            
    friend class boost::iterator_core_access;
            
    //@{
    //! \name Iterator Specific
            
    void increment();
            
    /*! \brief Returns true if:
     * both iterators are in the "end iterator" state
     * both iterators are NOT in the "end iterator" state
     * and they point to the same node of the same sub paving
     */
    bool equal( GridTreeConstIterator const & theOtherIterator) const;
            
    GridCell const& dereference() const;
            
    //@}

    //@{
    //! \name Local
            
    /*! \brief Allows to navigate to the first (\a firstLast==true ),
     * last (\a firstLast==false) enabled leaf of the sub paving
     * Returns true if the node was successfully found. If nothing is
     * found then the cursor should be in the root node again.
     */
    bool navigate_to(bool firstLast);
            
    /*! \brief A recursive search function, that looks for the next enabled leaf in the tree.
     *  The search is performed from left to right. (Is used for forward iteration)
     */
    void find_next_enabled_leaf();
            
    //@}
            
  public:
    //@{
    //! \name Constructors

    /*! \brief Default constructor constructs an invalid iterator. 
     *  \internal Note that a default constructor is required for compliance with the 
     *  STL Trivial Iterator concept. */
    GridTreeConstIterator();

    /*! \brief The constructor that accepts the subpacing \a pSubPaving to iterate on
     * The paramerter \a firstLastNone indicatges whether we want to position the iterator
     * on the first enabled leaf node (firstLastNone == true) or the last one (firstLastNone == false)
     * or we are constructing the "end iterator" that does not point anywhere.
     */
    explicit GridTreeConstIterator( const GridTreeSubset * pSubPaving, const tribool firstLastNone );
            
    /*! \brief The copy constructor that copies the current state by copying the
     *   underlying GridTreeCursor. The latter is copied by it's copy constructor.
     */
    GridTreeConstIterator( const GridTreeConstIterator& theGridPavingIter );
    
    //@}

    /*! \brief An assignment operator that copies the current state by copying the
     *  underlying GridTreeCursor. The latter is copied by it's assignment operator.
     */
    GridTreeConstIterator & operator=( const GridTreeConstIterator& theGridPavingIter );

            
    /*! This destructor only deallocates the GridTreeCursor, note that the
     * latter one does not deallocate the enclosed sub paving.
     */
    ~GridTreeConstIterator();

    //@{
    //! \name Get the cursor of the iterator
                
    //This cursor is only needed to get access to enable/disable node functionality
    GridTreeCursor const& cursor() const;
                
    //@}
};
    
/****************************************************************************************************/
/***************************************Inline functions*********************************************/
/****************************************************************************************************/

/****************************************BinaryTreeNode**********************************************/
    
inline void BinaryTreeNode::init( tribool isEnabled, BinaryTreeNode* pLeftNode, BinaryTreeNode* pRightNode ){
    _isEnabled = isEnabled;
    _pLeftNode = pLeftNode;
    _pRightNode = pRightNode;
}

inline BinaryTreeNode::BinaryTreeNode(const tribool isEnabled){
    init( isEnabled, NULL, NULL );
}
    
inline BinaryTreeNode::BinaryTreeNode(const BinaryTreeNode& theTreeNode){
    if( (theTreeNode._pLeftNode) != NULL && ( theTreeNode._pRightNode != NULL ) ){
        init( theTreeNode._isEnabled, new BinaryTreeNode( *theTreeNode._pLeftNode ), new BinaryTreeNode( *theTreeNode._pRightNode ) );
    }else{
        //NOTE: We do not allow for nodes where one leaf is NULL and another is not
        init( theTreeNode._isEnabled, NULL, NULL );
    }
}
    
inline BinaryTreeNode::BinaryTreeNode( const BooleanArray& theTree, const BooleanArray& theEnabledCells ) {
    //Make default initialization
    init( false, NULL, NULL ) ;
        
    //If the tree is not empry and there are enabled leafs then do the thing
    if( ( theTree.size() ) > 0 && ( theEnabledCells.size() > 0 ) ) {
        uint arr_index = 0, leaf_counter = 0;
        restore_node( this, arr_index, leaf_counter, theTree, theEnabledCells );
    } else {
        //Otherwise settle with one disabled node
        this->set_disabled();
    }
}

inline BinaryTreeNode::~BinaryTreeNode(){
    if( _pLeftNode != NULL ) {
        delete _pLeftNode;
    }
    if( _pRightNode != NULL ) {
        delete _pRightNode;
    }
}
    
inline bool BinaryTreeNode::is_enabled() const{
    return definitely(_isEnabled);
}

inline bool BinaryTreeNode::is_disabled() const{
    return ! possibly(_isEnabled) ;
}

inline bool BinaryTreeNode::is_leaf() const{
    return (_pLeftNode == NULL) && (_pRightNode == NULL);
}
    
inline BinaryTreeNode * BinaryTreeNode::left_node() const {
    return _pLeftNode;
}
    
inline BinaryTreeNode * BinaryTreeNode::right_node() const {
    return _pRightNode;
}
    
inline void BinaryTreeNode::set_enabled() {
    if ( is_leaf() ) {
        _isEnabled = true;
    } else {
        throw NotALeafNodeException(ARIADNE_PRETTY_FUNCTION);
    }
}

inline void BinaryTreeNode::copy_from( const BinaryTreeNode * pOtherNode ){
    if( pOtherNode != NULL ){
        _isEnabled = pOtherNode->_isEnabled;
        if( _pLeftNode != NULL){ delete _pLeftNode; _pLeftNode = NULL; }
        if( _pRightNode != NULL){ delete _pRightNode; _pRightNode = NULL; }
        if( pOtherNode->_pLeftNode != NULL ){ _pLeftNode = new BinaryTreeNode( * (pOtherNode->_pLeftNode) ); }
        if( pOtherNode->_pRightNode != NULL ){ _pRightNode = new BinaryTreeNode( * (pOtherNode->_pRightNode) ); }
    }
}

inline void BinaryTreeNode::set_disabled() {
    if ( is_leaf() ) {
        _isEnabled = false;
    } else {
        throw NotALeafNodeException(ARIADNE_PRETTY_FUNCTION);
    }
}
    
inline void BinaryTreeNode::set_unknown() {
    if ( ! is_leaf() ) {
        _isEnabled = indeterminate;
    } else {
        throw IsALeafNodeException(ARIADNE_PRETTY_FUNCTION);
    }
}

inline void BinaryTreeNode::make_leaf(tribool is_enabled ){
    _isEnabled = is_enabled;
    if( _pLeftNode != NULL ) { delete _pLeftNode; _pLeftNode= NULL; }
    if( _pRightNode != NULL ) { delete _pRightNode; _pRightNode= NULL; }
}
    
inline void BinaryTreeNode::split() {
    if ( is_leaf() ) {
        _pLeftNode  = new BinaryTreeNode(_isEnabled);
        _pRightNode = new BinaryTreeNode(_isEnabled);
        set_unknown();
    }
}
    
inline void BinaryTreeNode::add_enabled( const BinaryWord& path ){
    add_enabled( this, path, 0 );
}
    
inline void BinaryTreeNode::mince(const uint depth) {
    mince_node(this, depth);
}
    
inline void BinaryTreeNode::recombine() {
    recombine_node(this);
}
    
inline string BinaryTreeNode::node_to_string() const {
    stringstream tmp_stream;
    tmp_stream << "BinaryTreeNode( isLeaf = " << is_leaf() << ", isEnabled = " << is_enabled() << ", isDisabled = " << is_disabled() << " )";
    return tmp_stream.str();
}

inline BinaryTreeNode& BinaryTreeNode::operator=( const BinaryTreeNode & otherNode ) {
    //Copy the node value
    _isEnabled = otherNode._isEnabled;
    
    //Deallocate memory for the children, if any
    if( _pLeftNode != NULL ){ delete _pLeftNode; _pLeftNode = NULL; }
    if( _pRightNode != NULL ){ delete _pRightNode; _pRightNode = NULL; }
    
    //Copy the children trees from the otherNode
    if( otherNode._pLeftNode != NULL ) { _pLeftNode = new BinaryTreeNode( * ( otherNode._pLeftNode ) ); }
    if( otherNode._pRightNode != NULL ) { _pRightNode = new BinaryTreeNode( * ( otherNode._pRightNode ) ); }

    return *this;
}

/********************************************GridTreeCursor***************************************/
    
inline GridTreeCursor::GridTreeCursor(  ) :
    _currentStackIndex(-1), _pSubPaving(0), _theCurrentGridCell( Grid(), 0, BinaryWord() ) {
}

inline GridTreeCursor::GridTreeCursor(const GridTreeCursor & otherCursor) :
    _currentStackIndex(otherCursor._currentStackIndex), _theStack(otherCursor._theStack), 
    _pSubPaving(otherCursor._pSubPaving), _theCurrentGridCell(otherCursor._theCurrentGridCell){
}

inline GridTreeCursor& GridTreeCursor::operator=(const GridTreeCursor & otherCursor) {
    _currentStackIndex = otherCursor._currentStackIndex;
    _theStack = otherCursor._theStack; 
    _pSubPaving = otherCursor._pSubPaving;
    _theCurrentGridCell = otherCursor._theCurrentGridCell;
    
    return *this;
}

inline GridTreeCursor::GridTreeCursor(const GridTreeSubset * pSubPaving) : 
    _currentStackIndex(-1), _pSubPaving(pSubPaving), _theCurrentGridCell( pSubPaving->cell() ) {
    //Remember that GridTreeSubset contains GridCell
        
    //Add the current node to the stack, since we are in it

    //IVAN S ZAPREEV:
    //NOTE: There is no need in allocating _theStack elements,
    //it is all done in the push method
    push(pSubPaving->_pRootTreeNode);
}

inline GridTreeCursor::~GridTreeCursor() {
    //IVAN S ZAPREEV:
    //WARNING: The subpaving should not be deallocated here!
    //There are no other data feels that we need to deallocate
    _pSubPaving = NULL;
}
    
inline void GridTreeCursor::push( BinaryTreeNode* pLatestNode ){
    //If we are out of free space, increase the array's capacity
    if( static_cast<uint>(_currentStackIndex +1) == ( _theStack.size()  ) ){
        _theStack.resize( _theStack.size() + STACK_SIZE_INCREMENT );
    }
        
    //The index for the tree node to be added
    _currentStackIndex += 1;
        
    //Put the tree node pointer into the stack
    _theStack[_currentStackIndex] = pLatestNode;

}
    
inline BinaryTreeNode* GridTreeCursor::pop( ){
    BinaryTreeNode* pLastNode = NULL;
        
    //If there are non-root nodes in the stack
    if( _currentStackIndex > 0 ){
        //Return the stack element at _currentStackIndex,
        //then decrement the current element's index.
        return _theStack[ _currentStackIndex-- ];
    }
        
    return pLastNode;
}
    
inline bool GridTreeCursor::is_at_the_root() const {
    //If we are pointing at the zero cell of the array then it
    //means that we are in the root of the tree
    return (_currentStackIndex == 0);
}

inline bool GridTreeCursor::is_enabled() const {
    return _theStack[ _currentStackIndex ]->is_enabled();
}
    
inline bool GridTreeCursor::is_disabled() const {
    return _theStack[ _currentStackIndex ]->is_disabled();
}

inline bool GridTreeCursor::is_leaf() const {
    return _theStack[ _currentStackIndex ]->is_leaf();
}
    
inline bool GridTreeCursor::is_root() const {
    return  is_at_the_root();
}

inline void GridTreeCursor::set_enabled() const {
    return _theStack[ _currentStackIndex ]->set_enabled();
}
    
inline void GridTreeCursor::set_disabled() const {
    return _theStack[ _currentStackIndex ]->set_disabled();
}

inline bool GridTreeCursor::is_left_child() const{
    //If there is a parent node and the given node is it's left child
    return ( _currentStackIndex > 0 ) && ( _theStack[ _currentStackIndex - 1 ]->left_node() == _theStack[ _currentStackIndex ] );
}
    
inline bool GridTreeCursor::is_right_child() const{
    //If there is a parent node and the given node is it's right child
    return ( _currentStackIndex > 0 ) && ( _theStack[ _currentStackIndex - 1 ]->right_node() == _theStack[ _currentStackIndex ] );
}

inline GridTreeCursor& GridTreeCursor::move_up() {
    if( ! is_root() ){
        //Remove the current node from the stack
        pop();
        //Recompute the PavingGridCell
        updateTheCurrentGridCell( indeterminate );
        //Return the object back
        return ( * this);
    } else {
        throw NotAllowedMoveException(ARIADNE_PRETTY_FUNCTION);
    }
}
    
inline GridTreeCursor& GridTreeCursor::move_left() {
    return move(false);
}
    
inline GridTreeCursor& GridTreeCursor::move_right() {
    return move(true);
}
    
inline void GridTreeCursor::updateTheCurrentGridCell( tribool left_or_right ){
    if( ! indeterminate(left_or_right) ){
        //Here left_or_right is either true or false, also true defines moving to the right branch
        _theCurrentGridCell._theWord.push_back( definitely( left_or_right ) );
    } else {
        //If left_or_right is indeterminate, this means that we go "up"
        //in the tree, so we remove the last bit of the path.
        _theCurrentGridCell._theWord.pop_back();
    }
    _theCurrentGridCell = GridCell( _theCurrentGridCell._theGrid, _theCurrentGridCell._theHeight, _theCurrentGridCell._theWord );
}

inline GridTreeCursor& GridTreeCursor::move(bool left_or_right) {
    BinaryTreeNode* pNextNode;
    if( ! is_leaf() ){
        //If we are not in the leaf node then we can go down
        if( left_or_right ){ //true moves us to the right
            pNextNode = _theStack[ _currentStackIndex ]->right_node();
        } else { //false moves us to the left
            pNextNode = _theStack[ _currentStackIndex ]->left_node();
        }
        //Put the node into the stack
        push(pNextNode);
        //Recompute the PavingGridCell
        updateTheCurrentGridCell( left_or_right );
        //Return the object back
        return ( * this);
    } else {
        throw NotAllowedMoveException(ARIADNE_PRETTY_FUNCTION);
    }
}

inline const GridCell& GridTreeCursor::cell() const {
    return _theCurrentGridCell;
}
    
inline bool GridTreeCursor::operator==(const GridTreeCursor& anotherGridTreeCursor) const {
    bool areEqual = false;
    if( (this->_currentStackIndex >=0) && (anotherGridTreeCursor._currentStackIndex >=0 ) ){
        areEqual = (this->_theStack[this->_currentStackIndex] == anotherGridTreeCursor._theStack[anotherGridTreeCursor._currentStackIndex]);
    }
    return areEqual;
}
    
inline GridTreeSubset GridTreeCursor::operator*() {
    //IVAN S ZAPREEV:
    //NOTE: The first three parameters define the location of the _theStack[ _currentStackIndex ]
    //node with respect to the primary cell of the GridTreeSet.
    return GridTreeSubset( _theCurrentGridCell._theGrid,
                           _theCurrentGridCell._theHeight,
                           _theCurrentGridCell._theWord,
                           _theStack[ _currentStackIndex ] );
}
    
inline const GridTreeSubset GridTreeCursor::operator*() const {
    return (const GridTreeSubset & ) *(* this);
}

/****************************************GridTreeConstIterator************************************/
    
inline GridTreeConstIterator::GridTreeConstIterator( const GridTreeConstIterator& theGridPavingIter ) : 
    _is_in_end_state(theGridPavingIter._is_in_end_state), _pGridTreeCursor(theGridPavingIter._pGridTreeCursor) {
    //Copy everything
}

inline GridTreeConstIterator& GridTreeConstIterator::operator=( const GridTreeConstIterator& theGridPavingIter ){
    _is_in_end_state = theGridPavingIter._is_in_end_state;
    _pGridTreeCursor = theGridPavingIter._pGridTreeCursor;
    
    return *this;
}

inline GridTreeConstIterator::GridTreeConstIterator( const GridTreeSubset * pSubPaving, const tribool firstLastNone ):
    _pGridTreeCursor(pSubPaving) {
    if( ! indeterminate( firstLastNone ) ){
        //If the first/last enabled node is not found, it means that there are no elements
        //to iterate on, then we switch to the "end iterator" state
        _is_in_end_state = ! navigate_to( definitely( firstLastNone ) );
    } else {
        //In this case if we do nothing, the cursor in the iterator will point to the root node
        //Since this can be the only node in the tre we should add a marker that indicates that
        //the iterator is at the end state.
        _is_in_end_state = true;
    }
}
    
inline void GridTreeConstIterator::increment() {
    //If we are not done iterating
    if( ! _is_in_end_state){
        //We are at some enabled-leaf node and we want to find the next one.
        //The next node is somewhere on the right in the tree.
        find_next_enabled_leaf();
    }
}
    
inline GridTreeConstIterator::~GridTreeConstIterator() {
    //IVAN S ZAPREEV:
    //WARNING: There are no memory allocations in this class that have to be cleaned
}
    
inline bool GridTreeConstIterator::equal( GridTreeConstIterator const & theOtherIterator) const {
    //Check if both iterators are in the "end iterator" state
    bool result = theOtherIterator._is_in_end_state && this->_is_in_end_state;
        
    //If not then check if the cursors are equal (i.e. if they point to the same binary-tree node of the same (subpaving)
    if( ! result ){
        if( theOtherIterator._is_in_end_state || this->_is_in_end_state ){
            //If at one is in the end state and the other is not then the answer is FALSE
            result = false;
        } else {
            result = theOtherIterator._pGridTreeCursor == this->_pGridTreeCursor;
        }
    }
        
    return result;
}
    
inline GridCell const& GridTreeConstIterator::dereference() const {
    return _pGridTreeCursor.cell();
}
    
inline GridTreeCursor const& GridTreeConstIterator::cursor() const {
    return _pGridTreeCursor;
}
 
/*****************************************GridAbstractCell*******************************************/

inline GridAbstractCell::GridAbstractCell(const GridAbstractCell& theGridCell):
    _theGrid(theGridCell._theGrid), _theHeight(theGridCell._theHeight),
    _theWord(theGridCell._theWord), _theBox(theGridCell._theBox) {
}

inline GridAbstractCell::GridAbstractCell(const Grid& theGrid, const uint theHeight,
                                            const BinaryWord& theWord, const Box& theBox):
    _theGrid(theGrid), _theHeight(theHeight),
    _theWord(theWord), _theBox(theBox) {
}
    
inline void GridAbstractCell::primary_cell_at_height( const uint theHeight, int & leftBottomCorner, int & rightTopCorner ) {
    if ( theHeight % 2 == 1 ) {
        leftBottomCorner = 2*leftBottomCorner - rightTopCorner;
    } else {
        rightTopCorner   = 2*rightTopCorner - leftBottomCorner;
    }
}
    
inline const Grid& GridAbstractCell::grid() const {
    return _theGrid;
}
        
inline uint GridAbstractCell::height() const {
    return _theHeight;
}
    
inline int GridAbstractCell::depth() const {
    return _theWord.size() - _theHeight;
}
    
inline const BinaryWord& GridAbstractCell::word() const {
    return _theWord;
}
    
inline const Box& GridAbstractCell::box() const {
    return _theBox;
}

inline dimension_type GridAbstractCell::dimension() const {
    return _theGrid.dimension();
}

inline Interval GridAbstractCell::operator[](dimension_type i) const {
    return _theBox[i];
}

/*****************************************GridAbstractCell*******************************************/

inline GridAbstractCell& GridAbstractCell::operator=(const GridAbstractCell & otherCell ) {
    _theGrid = otherCell._theGrid;
    _theHeight = otherCell._theHeight;
    _theWord = otherCell._theWord;
    _theBox = otherCell._theBox;
    
    return *this;
}

/*********************************************GridCell***********************************************/

inline GridCell::GridCell() : GridAbstractCell( Grid(), 0, BinaryWord(), Box() ) {
}

inline GridCell::GridCell(const Grid& theGrid, const uint theHeight, const BinaryWord& theWord) :
                    GridAbstractCell( theGrid, theHeight, theWord, compute_box(theGrid, theHeight, theWord) ) {
}

inline GridCell& GridCell::operator=(const GridCell & otherCell ) {
    GridAbstractCell::operator=(otherCell);
    return *this;
}

inline bool GridCell::operator==( const GridCell& other ) const {
    return compare_grid_cells( this, other, COMPARE_EQUAL );
}      

inline bool GridCell::operator<(const GridCell& other) const {
    return compare_grid_cells( this, other, COMPARE_LESS );
}

inline bool GridCell::compare_grid_cells(const GridCell * pCellLeft, const GridCell &cellRight, const uint comparator ) {
    return GridAbstractCell::compare_abstract_grid_cells( pCellLeft, cellRight, comparator );
}

//The box is not related to the Grid, i.e. it is in the original space, whereas the binary tree is related to the Lattice of the grid:
inline Box GridCell::compute_box(const Grid& theGrid, const uint theHeight, const BinaryWord& theWord) {
    //Compute the lattice box and then map it to the original space using the grid data
    return lattice_box_to_space( GridCell::compute_lattice_box( theGrid.dimension(), theHeight, theWord ), theGrid );
}

/*********************************************GridOpenCell***********************************************/

inline GridOpenCell::GridOpenCell(const Grid& theGrid, const uint theHeight, const BinaryWord& theWord, const Box& theBox) :
                            GridAbstractCell(theGrid, theHeight, theWord, theBox) {
}

inline GridOpenCell::GridOpenCell() : GridAbstractCell( Grid(), 0, BinaryWord(), compute_box(Grid(), 0, BinaryWord()) ) {
}

inline GridOpenCell::GridOpenCell(const Grid& theGrid, const uint theHeight, const BinaryWord& theWord) :
                    GridAbstractCell( theGrid, theHeight, theWord, compute_box(theGrid, theHeight, theWord) ) {
}

inline GridOpenCell& GridOpenCell::operator=(const GridOpenCell & otherCell ) {
    GridAbstractCell::operator=(otherCell);
    return *this;
}

inline bool GridOpenCell::operator==( const GridOpenCell& other ) const {
    return compare_grid_cells( this, other, COMPARE_EQUAL );
}      

inline bool GridOpenCell::operator<(const GridOpenCell& other) const {
    return compare_grid_cells( this, other, COMPARE_LESS );
}

inline bool GridOpenCell::compare_grid_cells(const GridOpenCell * pCellLeft, const GridOpenCell &cellRight, const uint comparator ) {
    return GridAbstractCell::compare_abstract_grid_cells( pCellLeft, cellRight, comparator );
}

inline bool GridOpenCell::overlap( const GridOpenCell & theLeftOpenCell, const GridOpenCell & theRightOpenCell ) {
    return definitely( theLeftOpenCell.box().overlaps(  theRightOpenCell.box() ) );
}

/********************************************GridTreeSubset******************************************/
    
inline GridTreeSubset::GridTreeSubset( const Grid& theGrid, const uint theHeight,
                                       const BinaryWord& theWord, BinaryTreeNode * pRootTreeNode ) :
                                       _pRootTreeNode(pRootTreeNode), _theGridCell(theGrid, theHeight, theWord) {
}

inline GridTreeSubset::GridTreeSubset( const GridTreeSubset &otherSubset ) : _pRootTreeNode(otherSubset._pRootTreeNode),
                                                                             _theGridCell(otherSubset._theGridCell) {
}

inline GridTreeSubset::~GridTreeSubset() {
    //IVAN S ZAPREEV:
    //WARNING: This method should have no implementation what so ever
    //All the synamically allocatged data should be destroyed from the
    //corresponding Paving object
}

inline uint GridTreeSubset::compute_number_subdiv( Float theWidth, const Float theMaxWidth) const{
    //Compute the minimum number of subdivisions N as:
    //   minimum N : ( theWidth / 2^{N} ) <= theMaxWidth
    //This is equivalent to (because all values are positive)
    //   minimum N : theWidth / theMaxWidth <= 2^{N}
    //log is a continuous-increasing dunction, thus the latter is 
    //equivalent to (since 0 < theMaxWidth, theWidth < infinite )
    //   minimum N : log( theWidth / theMaxWidth ) <= N
    //In computations then we need to take
    //   N = ceil( log( theWidth / theMaxWidth ) )
        
    //IVAN S ZAPREEV:
    //NOTE: We take log_approx, log_approx because if we use _up and _down versions to maximize the answer,
    //then it might become very inaccurate, and in case require more subdivisions than needed.
    //NOTE: We need base 2 logarithm, we get it by the rule: log_a(b)/log_a(c) = log_c(a)
    //PIETER COLLINS:
    //NOTE: Float now uses approximate operators by default
    uint result = 0;
    if ( theWidth > theMaxWidth ){
        //result = (uint) ceil( div_approx( log_approx( div_approx( theWidth, theMaxWidth ) ) , log_approx( R(2.0) ) ) );
        result = (uint) ceil( Ariadne::div( Ariadne::log( Ariadne::div( theWidth, theMaxWidth ) ) , Ariadne::log( 2.0 ) ) );
    }
    return result;
}

inline bool GridTreeSubset::empty() const {
    return BinaryTreeNode::count_enabled_leaf_nodes( this->binary_tree() ) == 0; 
}
    
inline size_t GridTreeSubset::size() const {
    return BinaryTreeNode::count_enabled_leaf_nodes( this->binary_tree() ); 
}
    
inline dimension_type GridTreeSubset::dimension( ) const {
    return grid().dimension();
}


inline const Grid& GridTreeSubset::grid() const {
    return this->_theGridCell.grid();
}
    
inline GridCell GridTreeSubset::cell() const {
    return _theGridCell;
}
    
inline Box GridTreeSubset::bounding_box() const {    
    if(this->empty()) return Box(this->dimension());
    
    GridTreeSet::const_iterator iter=this->begin();
    Box bbox = iter->box();
     
    for( ; iter!=this->end(); ++iter) {
        Box cell = iter->box();
        for(uint i = 0; i < cell.dimension(); ++i) {
            if(cell[i].lower() < bbox[i].lower()) bbox[i].set_lower(cell[i].lower());
            if(cell[i].upper() > bbox[i].upper()) bbox[i].set_upper(cell[i].upper());
        }
    }
    
    return bbox;

}
    
inline void GridTreeSubset::mince( const uint theNewDepth ) {
    _pRootTreeNode->mince( theNewDepth );
}

inline void GridTreeSubset::recombine() {
    _pRootTreeNode->recombine();
}

inline uint GridTreeSubset::depth() const {
    return _pRootTreeNode->depth();
}

inline bool GridTreeSubset::operator==(const GridTreeSubset& anotherGridTreeSubset) const {
    return ( this->_theGridCell == anotherGridTreeSubset._theGridCell ) &&
        ( ( * this->_pRootTreeNode ) == ( * anotherGridTreeSubset._pRootTreeNode ) );
}
    
inline GridTreeSubset::const_iterator GridTreeSubset::begin() const {
    return GridTreeSubset::const_iterator(this, true);
}

inline GridTreeSubset::const_iterator GridTreeSubset::end() const {
    return GridTreeSubset::const_iterator(this, indeterminate);
}
    
inline const BinaryTreeNode * GridTreeSubset::binary_tree() const {
    return _pRootTreeNode;
}

inline tribool GridTreeSubset::superset( const Box& theBox ) const {
    //Simply check if theBox is covered by the set and then make sure that
    //all tree cells that are not disjoint from theBox are enabled

    ARIADNE_ASSERT( theBox.dimension() == cell().dimension() );

    tribool isASubSet = theBox.subset( cell().box() );
    if( ! isASubSet ) {
        //If the box is not covered by the cell corresponding to the root node
        //of the set, then clearly theBox is not a subset of this set.
        return false;
    } else {
        //Otherwise, is theBox is possibly a subset then we try to see furhter 
        BinaryWord pathCopy( cell().word() );
        return GridTreeSubset::covers( binary_tree(), grid(), cell().height(), pathCopy, theBox );    
    }
}

inline tribool GridTreeSubset::subset( const Box& theBox ) const {
    //Check that the box corresponding to the root node of the set
    //is not disjoint from theBox. If it is then the set is not a
    //subset of theBox otherwise we need to traverse the tree and check
    //if all it's enabled nodes give boxes that are subsets of theBox.
    
    ARIADNE_ASSERT( theBox.dimension() == cell().dimension() );

    BinaryWord pathCopy( cell().word() );
    
    return GridTreeSubset::subset( binary_tree(), grid(), cell().height(), pathCopy, theBox );    
}

inline tribool GridTreeSubset::disjoint( const Box& theBox ) const {
    //Simply check if the box does not intersect with the set

    ARIADNE_ASSERT( theBox.dimension() == cell().dimension() );

    BinaryWord pathCopy( cell().word() );
    
    return GridTreeSubset::disjoint( binary_tree(), grid(), cell().height(), pathCopy, theBox );    
}

inline tribool GridTreeSubset::overlaps( const Box& theBox ) const {
    //Check if the box of the root cell overlaps with theBox,
    //if not then theBox does not intersect with the cell,
    //otherwise we need to find at least one enabled node
    //in the binary tree, such that it's box overlaps theBox.

    ARIADNE_ASSERT( theBox.dimension() == cell().dimension() );

    BinaryWord pathCopy( cell().word() );
    
    return GridTreeSubset::overlaps( binary_tree(), grid(), cell().height(), pathCopy, theBox );    
}

inline GridTreeSubset& GridTreeSubset::operator=( const GridTreeSubset &otherSubset) {
    _pRootTreeNode = otherSubset._pRootTreeNode;
    _theGridCell = otherSubset._theGridCell;
    
    return *this;
}

/*********************************************GridTreeSet*********************************************/
    
inline void GridTreeSet::adjoin( const GridCell& theCell ) {
    ARIADNE_ASSERT( this->grid() == theCell.grid() );
        
    bool has_stopped = false;
    //Align the paving and the cell
    BinaryTreeNode* pBinaryTreeNode = align_with_cell( theCell.height(), true, false, has_stopped );
        
    //If we are not trying to adjoin something into an enabled sub cell of the paving
    if( ! has_stopped ){
        //Add the enabled cell to the binary tree.
        pBinaryTreeNode->add_enabled( theCell.word() );
    }
}
    
inline void GridTreeSet::adjoin( const GridTreeSubset& theOtherSubPaving ) {
    ARIADNE_ASSERT( this->grid() == theOtherSubPaving.cell().grid() );
        
    bool has_stopped = false;
    //Align the paving and the cell
    BinaryTreeNode* pBinaryTreeNode = align_with_cell( theOtherSubPaving.cell().height(), true, false, has_stopped );
        
    //If we are not trying to adjoin something into an enabled sub cell of the paving
    if( ! has_stopped ){
        //Now, the pBinaryTreeNode of this paving corresponds to the primary cell common with theOtherSubPaving.
        //The theOtherSubPaving's root node is defined the path theOtherSubPaving.word() which starts in the
        //node corresponding to the common primary cell.
        pBinaryTreeNode->add_enabled( theOtherSubPaving.binary_tree(), theOtherSubPaving.cell().word() );
    }
}

/**************************************FRIENDS OF BinaryTreeNode***************************************/

/*! \brief Stream insertion operator, prints out two binary arrays, one is the tree structure
 *  and the other is the true/false (enabled/disabled) values for the leaf nodes
 */    
inline std::ostream& operator<<(std::ostream& output_stream, const BinaryTreeNode & binary_tree ) {
    BinaryWord tree, leaves;
    binary_tree.tree_to_binary_words( tree, leaves );
    return output_stream << "BinaryTreeNode( Tree: " << tree << ", Leaves: " << leaves << ")";
}

/****************************************FRIENDS OF GridOpenCell*******************************************/

inline std::ostream& operator<<(std::ostream& os, const GridOpenCell& theGridOpenCell ) {
    //Write the grid data to the string stream
    return os << "GridOpenCell( " << theGridOpenCell.grid() <<
        ", Primary cell height: " << theGridOpenCell.height() << 
        ", Base cell path: " << theGridOpenCell.word() <<
        ", Box (closure): " << theGridOpenCell.box() << " )";
}

/****************************************FRIENDS OF GridCell*******************************************/

/*! \brief Stream insertion operator for the GridCell. */
inline std::ostream& operator<<(std::ostream& os, const GridCell& gridPavingCell){
    //Write the grid data to the string stream
    return os << "GridCell( " << gridPavingCell.grid() <<
        ", Primary cell height: " << gridPavingCell.height() << 
        ", Path to the root: " << gridPavingCell.word() <<
        ", Box: " << gridPavingCell.box() << " )";
}

/*************************************FRIENDS OF GridTreeCursor*****************************************/

inline std::ostream& operator<<(std::ostream& os, const GridTreeCursor& theGridTreeCursor){
    const int curr_stack_idx = theGridTreeCursor._currentStackIndex;
    os << "GridTreeCursor( " << theGridTreeCursor._pSubPaving <<
        ", Curr. stack index: " << curr_stack_idx  <<
        ", Stack data: [ ";
    for(int i = 0; i <= curr_stack_idx; i++){
        os << theGridTreeCursor._theStack[i]->node_to_string() << ( ( i < curr_stack_idx ) ? "" : ", ");
    }
    return os<<" ], " << theGridTreeCursor._theCurrentGridCell << " )";
}

/*************************************FRIENDS OF GridTreeSubset*****************************************/

inline std::ostream& operator<<(std::ostream& os, const GridTreeSubset& theGridTreeSubset) {
    return os << "GridTreeSubset( Primary cell: " << theGridTreeSubset.cell() << ", " << (*theGridTreeSubset.binary_tree()) <<" )";
}
    
/***************************************FRIENDS OF GridTreeSet******************************************/
    
inline std::ostream& operator<<(std::ostream& os, const GridTreeSet& theGridTreeSet) {
    const GridTreeSubset& theGridTreeSubset = theGridTreeSet;
    return os << "GridTreeSet( " << theGridTreeSubset << " )";
}
    
inline GridTreeSet outer_approximation(const CompactSetInterface& theSet, const uint depth) {
    Grid theGrid(theSet.dimension());
    return outer_approximation(theSet,theGrid,depth);
}
    
inline GridTreeSet outer_approximation(const CompactSetInterface& theSet, const Grid& theGrid, const uint depth) {
    GridTreeSet result(theGrid);
    result.adjoin_outer_approximation(theSet,depth);
    return result;
}


template<class BS>
GridTreeSet outer_approximation(const ListSet<BS>& theSet, const Grid& theGrid, const uint depth) {
    ARIADNE_ASSERT_MSG(theSet.dimension()==theGrid.dimension(),"theSet="<<theSet<<", theGrid="<<theGrid);
    GridTreeSet result(theGrid);
    for(typename ListSet<BS>::const_iterator iter=theSet.begin(); iter!=theSet.end(); ++iter) {
        result.adjoin_outer_approximation(*iter,depth);
    }
    result.recombine();
    return result;
}


template<class A> void serialize(A& archive, Ariadne::GridTreeSet& set, const unsigned int version) {
    ARIADNE_NOT_IMPLEMENTED;
}

void draw(GraphicsInterface& theGraphic, const GridCell& theGridCell); 
    
void draw(GraphicsInterface& theGraphic, const GridTreeSet& theGridTreeSet);

void draw(GraphicsInterface& theGraphic, const CompactSetInterface& theSet);


} // namespace Ariadne

#endif /* ARIADNE_GRID_SET_H */

