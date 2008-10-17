/***************************************************************************
 *            grid_paving.h
 *
 *  Copyright  2008  Ivan S. Zapreev, Pieter Collins
 *  ivan.zapreev@gmail.com, Pieter.Collins@cwi.nl
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

/*! \file grid_paving.h
 *  \brief Grid paving is used to represent sets, based on integer and dyadic coordinate cells, of a grid.
 */

#ifndef GRID_PAVING_H
#define GRID_PAVING_H

#include <iostream>
#include <string>

#include "boost/iterator/iterator_facade.hpp"

#include "tribool.h"
#include "array.h"

#include "binary_word.h"

#include "grid_set.h"
#include "exceptions.h"
#include "box.h"
#include "point.h"

#include "vector.h"


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
	class GridPavingCell;
	class GridPavingCursor;
	class GridPavingConstIterator;
	class GridSubPaving;
	class GridSubPavingIterator;
	class GridPaving;
	
	std::ostream& operator<<(std::ostream& os, const GridPavingCell& gridPavingCell);
	bool subset(const GridPavingCell& theCell, const GridSubPaving& theSet);
	bool overlap(const GridPavingCell& theCell, const GridSubPaving& theSet);
	bool subset(const GridSubPaving& theSet1, const GridSubPaving& theSet2);
	bool overlap(const GridSubPaving& theSet1, const GridSubPaving& theSet2);
	bool subset(const Box& theBox, const GridSubPaving& theSet);
	bool disjoint(const Box& theBox, const GridSubPaving& theSet);
	bool intersects(const Box& theBox, const GridSubPaving& theSet);
	
	GridPaving join(const GridSubPaving& theSet1, const GridSubPaving& theSet2);
	GridPaving intersection(const GridSubPaving& theSet1, const GridSubPaving& theSet2);
	GridPaving difference(const GridSubPaving& theSet1, const GridSubPaving& theSet2);
	
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
	
	/*! \brief The GridPavingCursor throws this exception if we try to go beyond the binary tree. */
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
	class BinaryTreeNode{
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
			void restore_node( BinaryTreeNode * theCurrentNode, uint & arr_index, uint & leaf_counter,
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
			//! \name Operations
			
			/*! \brief This method recombined the sub tree nodes. Note that, the two leaf nodes with
			 * the same parent are removed if they have the same value of isEnabled fields.
			 */
			void recombine();
			
			/*! \brief Creates a string representing the binary tree rooted to the given node */
			string tree_to_string() const;

			/*! \brief Creates a string representing the given node */
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
			
			/*! \brief Starting in the \a pRootTreeNode node as at the root, this method finds(creates)
			 *  the leaf node defined by the \a path and marks it as enabled. If some prefix of the \a path
			 *  references an enabled node then nothing is done.
			 *  NOTE: This is a recursive method on the position (\a position) in the binary path (\path).
			 *  Therefore, the initial call of this method should be done with \a position == 0;
			 */
			static void add_enabled( BinaryTreeNode* pRootTreeNode, const BinaryWord& path, const uint position = 0 );

			/*! \brief Creates a binary tree of the height rootNodePath.size(), puts the subtree oldRootNode
			 * into the node defined by the path \a rootNodePath, returns the root node of the extended tree.
			 */
			static BinaryTreeNode * prepend_tree( const BinaryWord & rootNodePath, BinaryTreeNode * oldRootNode);
			
			//@}
	};

	/*! \brief A cell of a grid paving.
	 * 
	 * This class contains the geometric data of a cell, as well as the combinatoric data needed to easily adjoin it to a tree.
	 * It does not contain the tree structure, so cannot be used in cursors/iterators.
	 * NOTE: The "Primary root cell" is the cell that the \a GridPaving can be rooted to
	 * (the \a GridSubPaving can be rooted to a non primary cell).
	 */
	class GridPavingCell
        // : public BoxExpression< GridPavingCell > 
        {
		protected:
			/*! \brief The intrinsic grid of the cell. Note that grids are internally passed by reference. */
			Grid _theGrid;

			/*! \brief The level of the primary root cell relative to the zero level. */
			uint _theHeight;

			/*! \brief The word describing the path in a binary tree to the given cell from the primary root cell. */
			BinaryWord _theWord;

			/*! \brief The geometric box represented by the cell. */
			Box _theBox;
			
			/*! \brief The box is given as an explicit parameter. This should only be used
			 *  by friend classes which have already computed the box.
			 */
			GridPavingCell(const Grid& theGrid, const uint theHeight, const BinaryWord& theWord, const Box& theBox);
			
			friend class GridPavingCursor;

			/*! \brief having \a theHeight the hight of the primary cell, with \a leftBottomCorner and \a rightTopCorner
			 * defining the primary cell corners for \a theHeight-1, we recompute \a leftBottomCorner and \a rightTopCorner
			 * for the current level \a theHeight.
			 */
			static inline void primary_cell_at_height( const uint theHeight, int & leftBottomCorner, int & rightTopCorner );

		public:
			/*! \brief The copy constructor for the \a GridPavingCell */
			GridPavingCell(const GridPavingCell& theGridPavingCell);

			/*! \brief Construct a cell based on \a theGrid, the primary cell of level \a theHeight,
			 * and the path to the SubPaving's root cell which is accessible from the primary cell
			 * via the path \a _theWord.
			 */
			GridPavingCell(const Grid& theGrid, const uint theHeight, const BinaryWord& theWord);
			
			/*! \brief The underlying grid. */
			const Grid& grid() const;
			
			/*! \brief The height of the tree above the zero level. */
			const uint& height() const;
			
			/*! \brief The word describing the path in a binary tree to the cell. */
			const BinaryWord& word() const;
			
			/*! \brief The geometric box represented by the cell. */
			const Box& box() const;

			/*! \brief The dimension of the cell. Needed by the BoxExpression interface. */
			dimension_type dimension() const;

			/*! \brief The upper and lower bound in the \a i<sup>th</sup> coordinate. */
			Interval operator[](dimension_type i) const;
			
			/*! \brief Serialize the current object state, this is mostly needed for testing */
			string to_string() const;

			/*! \brief this method computes the box in the original space based on the \a theGrid,
			 *  and a cell which is obtained by traversing the path given by \a _theWord from the
			 * primary cell located at the heigth \a theHeight above the zero level cells
			 * (relative to the grid).
			 */
			static Box compute_box(const Grid& theGrid, const uint theHeight, const BinaryWord& theWord);
			
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
			 * \a theHeight defines the higth of the primary cell above the zero level and
			 * \a dimensions is the number of dimension in the considered space
			 */
			static Box primary_cell( const uint theHeight, const dimension_type dimensions );
			
			/*! \brief Takes a box \a theBox related to some grid and computes the smallest primary cell on
			 *   that grid that will contain this box. Here \a dimensions defines the space dimension and
			 *   the method returns the hight of that primary cell.
			 */
			static uint smallest_primary_cell( const dimension_type dimensions, const Box& theBox );

			/*! \brief Apply grid data \a theGrid to \a theGridBox in order to compute
			 * the box dimensions in the original space
			 */
			static Box grid_box_to_space(const Box & theGridBox, const Grid& theGrid );

			/*! \brief This method returns the path from the \a topPCellHeight to the \a bottomPCellHeight
			 *  in the \a dimensions dimensional space. We assume that \a topPCellHeight >= \a bottomPCellHeight
			 *  if not, then we return an empty binary word.
			 */
			static BinaryWord primary_cell_path( const uint dimensions, const uint topPCellHeight, const uint bottomPCellHeight);
	};

	/*! \brief This class represents a subpaving of a paving. Note that, the subtree enclosed into
	 * this class is just a pointer to the node in the tree of some paving. This class is not
	 * responsible for deallocation of that original tree.
	 */
	class GridSubPaving {
		protected:
			
			friend class GridPavingCursor;

			template<class Set> friend GridPaving outer_approximation( const Grid& theGrid, const Set& theSet, const uint depth );
			template<class Set> friend GridPaving inner_approximation( const Grid& theGrid, const Set& theSet, const uint depth );
			
			/*! \brief The pointer to the root node of the subpaving tree.
			 * Note that, this is not necessarily the root node of the corresponding paving tree.
			 */
			BinaryTreeNode * _pRootTreeNode;
			
			/*! \brief The paving cell corresponding to the root node of the SubPaving.*/
			GridPavingCell _theGridPavingCell;

			/*! \brief this function takes the interval width and computes how many binary subdivisions
			 * one has to make in order to have sub-intervals of the width <= \a theMaxWidth
			 */
			uint compute_number_subdiv( Float theWidth, const Float theMaxWidth) const;
		
		public:
			/*! \brief A short name for the constant iterator */
			typedef GridPavingConstIterator const_iterator;
			
			//@{
			//! \name Constructors
			
			/*! \brief The new root node can only be constructed from the existing tree node.
			 * Here, \a pRootTreeNode is not copied, we simply store its pointer inside this class.
			 * Note that, \a pRootTreeNode should correspond to the sub-paving root node. Thus,
			 * \a theHeight defines the higth of the primary root cell of the GridPaving
			 * (Remember that every GridSubPaving is just a reference to a subtree of a GridPaving).
			 * \a theWord defines the path to the \a pRootTreeNode node from the primary root cell
			 * of the corresponding GridPaving.
			 */
			GridSubPaving( const Grid& theGrid, const uint theHeight, const BinaryWord& theWord, BinaryTreeNode * pRootTreeNode );
			
			//@}
			
                        /*! Virtual destructor. The destructor needs to be virtual since GridPaving is a subclass 
                         *  with different memory management. */
                        virtual ~GridSubPaving();
			
			//@{
			//! \name Properties
			
			/*! \brief Returns the const pointer to the root BinaryTreeNode of the SubPaving*/
			const BinaryTreeNode * binary_tree() const;
			
			/*! Recalculate the depth of the tree rooted at \a _pRootTreeNode */
			uint depth() const;
			
			/*! \brief Serialize the current object state, this is mostly needed for testing */
			virtual string to_string() const;

			/*! \brief Returns the \a GridPavingCell corresponding to the ROOT NODE of this \a GridSubPaving
			 * WARNING: It is NOT the primary cell of the paving! 
			 */
			GridPavingCell cell() const;
			
			/*! \brief Allows to test if the two subpavings are equal, this is determined by
			 * the fact that they are rooted to the same binary tree. In other words we check
			 * that the values of \a _pRootTreeNode are the same
			*/
			bool operator==(const GridSubPaving& anotherGridSubPaving) const;
			
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
			//! \name Iterators

			/*! \brief A constant iterator through the enabled leaf nodes of the subpaving. */
			const_iterator begin() const;
			
			/*! \brief A constant iterator to the end of the enabled leaf nodes of the subpaving. */
			const_iterator end() const;

			//@}

			//@{
			//! \name Conversions 
			/*  PIETER: These conversions are not strictly necessary, but may be useful. */
			
			/*! \brief Convert to a list of ordinary boxes, unrelated to the grid. */
			operator ListSet<Box>() const;
			
			/*! \brief Convert to a list of cells. */
			operator GridCellListSet() const;
			
			//@}
	};

	/* \brief The GridPaving class that represents a set of cells with mixed integer and dyadic coordinates.
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
	class GridPaving : public GridSubPaving {		
		protected:

			friend class GridPavingCursor;

			/*! \brief This method takes the primary cell of \a theCell and if it is:
			 *   (a) higher then for this paving, pre-pends \a _pRootTreeNodeis.
			 *   (b) lower then for this paving, locates it in this paging's tree.
			 *   (c) equal to the hight of this paving's primary cell, does nothing.
			 *  This method returns the node corresponding to the primary cell of
			 *   \a theCell in the (updated) paving.
			 *  If \a stop_on_enabled is set tro true then for the case (b) if we meet
			 *  a leaf node on the path from the primary onde of the paving to the primary
			 *  node of the cell, we stop locating the primary cell of \a theCell and set
			 *  \a has_stopped to true.
			 */
			BinaryTreeNode* align_with_cell( const GridPavingCell& theCell, const bool stop_on_enabled, bool & has_stopped );

			/*! \brief This method adjoins the outer approximation of \a theSet (computed on the fly) to this paving.
			 *  We use the primary cell (enclosed in this paving) of height \a primary_cell_hight and represented 
			 *  by the paving's binary node \a pBinaryTreeNode. When adding the outer approximation, we compute it
			 *  up to the level of accuracy given by \a max_mince_depth. This parameter defines, how many subdivisions
			 *  of the vinary tree we should make to get the proper cells for outer approximating \a theSet.
			 *  This method is recursive, the parameter \a pPath defines the path to the current node pBinaryTreeNode
			 *  from the root node in recursive calls, thus the initial call for this method must be done with an empty word.
			 */
			template<class Set> void adjoin_outer_approximation( BinaryTreeNode * pBinaryTreeNode, const uint primary_cell_height,
										const uint max_mince_depth, const Set& theSet, BinaryWord * pPath );

		public:
			//@{
			//! \name Constructors
			
			/*! \brief The new root node can only be constructed from the existing tree node.
			 * Here, \a pRootTreeNode is not copied, we simply store its pointer inside this class.
			 * Note that, \a pRootTreeNode should correspond to the root node. \a theHeight defines
			 * the higth of the primary root cell corresponding to the \a pRootTreeNode node.
			 */
			GridPaving( const Grid& theGrid, const uint theHeight, BinaryTreeNode * pRootTreeNode );
			
			/*! \brief The copy constructor that actually copies all the data,
			 * including the paving tree. I.e. the new copy of the tree is created
			 */
			GridPaving( const GridPaving & theGridPaving );

			/*! A simple constructor that creates the [0, 1]*...*[0, 1] cell in the
			 * \a theDimension - dimensional space. Here we assume that we have a non scaling
			 * grid with no shift of the coordinates. I.e. Grid._data._origin = {0, ..., 0}
			 * and Grid._data._lengths = {1, ..., 1}. If enable == true then the cell is enabled
			 */
			explicit GridPaving( const uint theDimension, const bool enable = false );

			/*! \brief Construct an empty tree. The \a theBoundingBox is used to define the lattice
			 *  block (in theGrid) that will correspond to the root of the paving tree \a pRootTreeNode.
			 */
			explicit GridPaving( const Grid& theGrid, const bool enable = false  );

			/*! \brief Construct an empty tree. The \a theBoundingBox is used to define the lattice
			 *  block (in theGrid) that will correspond to the root of the paving tree \a pRootTreeNode.
			 */
			explicit GridPaving( const Grid& theGrid, const Box & theBoundingBox );
			
			/*! \brief Construct the paving based on the block's coordinates, defined by: \a theLeftLowerPoint
			 * and \a theRightUpperPoint. These are the coordinates in the lattice defined by theGrid.
			 * The primary cell, enclosing the given block of cells is computed automatically and the
			 * binary tree is rooted to that cell. The \a theEnabledCells array defines the enabled/disabled
			 * 0-level cells of the block (lexicographic order).
			 */
			explicit GridPaving( const Grid& theGrid, const IndexArray theLeftLowerPoint,
						const IndexArray theRightUpperPoint, const BooleanArray& theEnabledCells );

                        /*! \brief Creates a new paving from the user data. \a theTree is an array representation of the binary
			 * tree structure, \a theEnabledCells tells whether a node is or is not a leaf, \a theHeight gives the
			 * higth of the primary cell which is assumed to correspond to the root node of \a theTree.
			 */
			explicit GridPaving( const Grid& theGrid, uint theHeight, const BooleanArray& theTree, const BooleanArray& theEnabledCells );

			//@}

			//@{
			//! \name Cloning
			
			/*! \brief Return a new dynamically-allocated copy of the GridPaving.
			 *  In this case, all the data is copied.
			 */
                         GridPaving* clone() const;

			//@}
			
			/*! \brief Destructor, removes all the dynamically allocated data, any */
			/* GridSubPaving referencing this GridPaving becomes invalid. */
			virtual ~GridPaving();
			
			//@{
			//! \name Geometric Predicates

			/*! \brief Tests if a cell is a subset of a set. */
			friend bool subset( const GridPavingCell& theCell, const GridSubPaving& theSet );
	
			/*! \brief Tests if a cell is overlaps (intersects as an open set) a paving set. */
			friend bool overlap( const GridPavingCell& theCell, const GridSubPaving& theSet );

			/*! \brief Tests if a grid paving set is a subset of another. */
			friend bool subset( const GridSubPaving& theSet1, const GridSubPaving& theSet2 );
	
			/*! \brief Tests if two grid paving sets overlap (i.e. intersect as open sets.) */
			friend bool overlap( const GridSubPaving& theSet1, const GridSubPaving& theSet2 );

			/* PIETER: We may want the two predicates below to return tribool, reflecting the case that the box is an
			* approximation. Unfortunately, this is very difficult for subset. */
			/*! \brief Tests if a box is a subset of a set. */
			friend bool subset( const Box& theBox, const GridSubPaving& theSet );
	
			/*! \brief Tests if a box is disjoint from (the closure of) a grid set. */
			friend bool disjoint( const Box& theBox, const GridSubPaving& theSet );
	
			/*! \brief Tests if a box intersects thea set. */
			friend bool intersects( const Box& theBox, const GridSubPaving& theSet );

			//@}

			//@{
			//! \name Geometric Operations
			/* PIETER: You may prefer to make inplace operations may return
			 *a reference to *this, allowing chaining of operations.
			 */

			/*! \brief Adjoin (make inplace union with) a single cell. */
			void adjoin( const GridPavingCell& theCell );
			
			/*! \brief Remove a single cell. */
			void remove( const GridPavingCell& theCell );

			/*! \brief Adjoin (make inplace union with) another grid paving set. */
			void adjoin( const GridSubPaving& theOtherSubPaving );
			
			/*! \brief Restrict to (make inplace intersection with) another grid paving set. */
			void restrict( const GridSubPaving& theOtherSubPaving );
			
			/*! \brief Remove cells in another grid paving set. */
			void remove( const GridSubPaving& theOtherSubPaving );

			/*! \brief Join (make union of) two grid paving sets. */
			friend GridPaving join( const GridSubPaving& theSet1, const GridSubPaving& theSet2 );
	
			/*! \brief The intersection of two grid paving sets. Points only lying on the
			 *  intersection of the boundaries of the two sets are not included in the result.
			 */
			friend GridPaving intersection( const GridSubPaving& theSet1, const GridSubPaving& theSet2 );
	
			/*! \brief The difference of two grid paving sets. */
			friend GridPaving difference( const GridSubPaving& theSet1, const GridSubPaving& theSet2 );

			//@}

			//@{
			//! \name Geometric Approximation

			/*! \brief Adjoin an outer approximation to a given set, computing to the given depth.
			 *! This method computes an outer approximation for the set \a theSet on the grid \a theGrid.
			 *  Note that, the depth is the total number of subdivisions (in all dimensions) of the unit
			 *  cell of the grid. This method does the followig:
			 * 1. Computes the smallest Primary cell enclosing \a theSet
			 * 2. Allocates the paving for this cell
			 * 3. Minces the paving to the level: depth + <the primary cell hight>
			 * 4. Iterates through the enabled leaf nodes of the paving (all the nodes are initially enabled)
			 * 5. Disables the cells that are disjoint with the \a theSet
			 */
			template<class Set> void adjoin_outer_approximation( const Set& theSet, const uint depth );
			
			/*! \brief Adjoin an inner approximation to a given set, computing to the given depth. */
			template<class Set> void adjoin_inner_approximation( const Set& theSet, const uint depth );
			
			/*! \brief Adjoin a lower approximation to a given set, computing to the given depth. 
			*   A lower approximation comprises all cells intersecting a given set. */
			template<class Set> void adjoin_lower_approximation( const Set& theSet, const uint depth );

			//@}
	};

	/*! \brief This class represents a cursor/iterator that can be used to traverse a subtree. */
	/* PIETER: It might be better to have the GridPavingCell as a data field rather than a subclass. */
	class GridPavingCursor {
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
			const GridSubPaving * _pSubPaving;
			
			GridPavingCell _theCurrentGridCell;

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
			
		public:
			/*! \brief The constructor that accepts the subpaving to cursor on */
			GridPavingCursor(const GridSubPaving * pSubPaving);
			
			/*! \brief The simple copy copnstructor */
			GridPavingCursor(const GridPavingCursor & pSubPaving);
			
			/*! This destructor does not deallocate the enclosed sub paving. */
			~GridPavingCursor();

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
			GridPavingCursor& move_up();
			
			/*! \brief Move to the left child node. Throws an NotAllowedMoveException if the current node is a leaf. */
			GridPavingCursor& move_left();
			
			/*! \brief Move to the right child node. Throws an NotAllowedMoveException if the current node is a leaf. */
			GridPavingCursor& move_right();
			
			/*! \brief Move to a child node. Throws an NotAllowedMoveException if the current node is a leaf.
			 * if left_or_right == false then we go left, if left_or_right == true then right
			 */
			GridPavingCursor& move(bool left_or_right);

			/*! \brief Covert to a GridPavingCell. */
			const GridPavingCell& cell() const;
			
			/*! \brief Allows to test if the two cursors are equal, this is determined by
			 * the fact that they store the same nodes in \a _theStack[_currentStackIndex].
			 * In other words we check that they point to the same binary-tree node.
			 * NOTE: if _currentStackIndex < 0 for at least one of the cursors then the
			 * result is always false.
			*/
			bool operator==(const GridPavingCursor& anotherGridPavingCursor) const;

			/*! \brief The dereferencing operator which returns a
			 * reference to the GridSubPaving for the current node.
			 */
			GridSubPaving operator*();
			
			/*! \brief The dereferencing operator which returns a constant
			 * reference to the GridSubPaving for the current node.
			 */
			const GridSubPaving operator*() const;
			
			/*! \brief Serialize the current cursor state, this is mostly needed for testing */
			string to_string() const;
	};
	
	/*! \brief This class allows to iterate through the enabled leaf nodes of GridSubPaving.
	 * The return objects for this iterator are constant GridPavingCell
	 */
	class GridPavingConstIterator : public boost::iterator_facade< GridPavingConstIterator, GridPavingCell const, boost::forward_traversal_tag > {
		private:
			/*! \brief When set to true indicates that this is the "end iterator" */
			bool _is_in_end_state;
			
			/*! \brief We do not want to have the default constructor here */
			GridPavingConstIterator();
			
			/*! \brief the cursor object that is created based on the given sub paving*/
			GridPavingCursor _pGridPavingCursor;
			
			friend class boost::iterator_core_access;
			
			//@{
			//! \name Iterator Specific
			
			void increment();
			
			/*! \brief Returns true if:
			 * both iterators are in the "end iterator" state
			 * both iterators are NOT in the "end iterator" state
			 * and they point to the same node of the same sub paving
			 */
			bool equal( GridPavingConstIterator const & theOtherIterator) const;
			
			GridPavingCell const& dereference() const;
			
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
			
			/*! \brief The constructor that accepts the subpacing \a pSubPaving to iterate on
			 * The paramerter \a firstLastNone indicatges whether we want to position the iterator
			 * on the first enabled leaf node (firstLastNone == true) or the last one (firstLastNone == false)
			 * or we are constructing the "end iterator" that does not point anywhere.
			 */
			explicit GridPavingConstIterator( const GridSubPaving * pSubPaving, const tribool firstLastNone );
			
			/*! \brief The copy constructor */
			GridPavingConstIterator( const GridPavingConstIterator& gridPavingIter );
			
			//@}
			
			/*! This destructor only deallocates the GridPavingCursor, note that the
			 * latter one does not deallocate the enclosed sub paving.
			 */
			~GridPavingConstIterator();

			//@{
			//! \name Get the cursor of the iterator
				
			//This cursor is only needed to get access to enable/disable node functionality
			GridPavingCursor const& cursor() const;
				
			//@}
	};

	
} // namespace Ariadne

#include "grid_paving.inline.h"

#endif /* GRID_PAVING */

