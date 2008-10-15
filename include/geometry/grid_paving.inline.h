/***************************************************************************
 *            grid_paving.inline.h
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

#include "numeric/declarations.h"
#include "numeric/rational.h"
#include "numeric/float.h"
#include "numeric/interval.h"

namespace Ariadne {

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
			int arr_index = 0, leaf_counter = 0;
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
			throw NotALeafNodeException(__PRETTY_FUNCTION__);
		}
	}
	
	inline void BinaryTreeNode::set_disabled() {
		if ( is_leaf() ) {
			_isEnabled = false;
		} else {
			throw NotALeafNodeException(__PRETTY_FUNCTION__);
		}
	}
	
	inline void BinaryTreeNode::set_unknown() {
		if ( ! is_leaf() ) {
			_isEnabled = indeterminate;
		} else {
			throw IsALeafNodeException(__PRETTY_FUNCTION__);
		}
	}

	inline void BinaryTreeNode::make_leaf(tribool is_enabled ){
		//WARNING: No checks for NON-NULL pointer values, for performance reasons
		delete _pLeftNode; _pLeftNode= NULL;
		delete _pRightNode; _pRightNode= NULL;
		_isEnabled = is_enabled;
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
		string result = "";
		result += "isLeaf = ";
		result += is_leaf() ? "1" : "0";
		result += ", isEnabled = ";
		result += is_enabled() ? "1" : "0";
		result += ", isDisabled = ";
		result += is_disabled() ? "1" : "0";
		return  result;
	}

/********************************************GridPavingCursor***************************************/

	template<class R> inline GridPavingCursor<R>::GridPavingCursor(const GridPavingCursor<R> & theSubPaving) :
						_theStack(theSubPaving._theStack), _currentStackIndex(theSubPaving._currentStackIndex),
						_pSubPaving(theSubPaving._pSubPaving), _theCurrentGridCell(theSubPaving._theCurrentGridCell){
		//Copy everything
	}

	template<class R> inline GridPavingCursor<R>::GridPavingCursor(const GridSubPaving<R> * pSubPaving) : 
			_currentStackIndex(-1), _pSubPaving(pSubPaving), _theCurrentGridCell( pSubPaving->cell() ) {
		//Remember that GridSubPaving<R> contains GridPavingCell<R>
		
		//Add the current node to the stack, since we are in it

		//IVAN S ZAPREEV:
		//NOTE: There is no need in allocating _theStack elements,
		//it is all done in the push method
		push(pSubPaving->_pRootTreeNode);
	}

	template<class R> inline GridPavingCursor<R>::~GridPavingCursor() {
		//IVAN S ZAPREEV:
		//WARNING: The subpaving should not be deallocated here!
		//There are no other data feels that we need to deallocate
		_pSubPaving = NULL;
	}
	
	template<class R> inline void GridPavingCursor<R>::push( BinaryTreeNode* pLatestNode ){
		//If we are out of free space, increase the array's capacity
		if( _currentStackIndex == ( _theStack.size() -1 ) ){
			_theStack.resize( _theStack.size() + STACK_SIZE_INCREMENT );
		}
		
		//The index for the tree node to be added
		_currentStackIndex += 1;
		
		//Put the tree node pointer into the stack
		_theStack[_currentStackIndex] = pLatestNode;

	}
	
	template<class R> inline BinaryTreeNode* GridPavingCursor<R>::pop( ){
		BinaryTreeNode* pLastNode = NULL;
		
		//If there are non-root nodes in the stack
		if( _currentStackIndex > 0 ){
			//Return the stack element at _currentStackIndex,
			//then decrement the current element's index.
			return _theStack[ _currentStackIndex-- ];
		}
		
		return pLastNode;
	}
	
	template<class R> inline bool GridPavingCursor<R>::is_at_the_root() const {
		//If we are pointing at the zero cell of the array then it
		//means that we are in the root of the tree
		return (_currentStackIndex == 0);
	}

	template<class R> inline bool GridPavingCursor<R>::is_enabled() const {
		return _theStack[ _currentStackIndex ]->is_enabled();
	}
	
	template<class R> inline bool GridPavingCursor<R>::is_disabled() const {
		return _theStack[ _currentStackIndex ]->is_disabled();
	}

	template<class R> inline bool GridPavingCursor<R>::is_leaf() const {
		return _theStack[ _currentStackIndex ]->is_leaf();
	}
	
	template<class R> inline bool GridPavingCursor<R>::is_root() const {
		return  is_at_the_root();
	}

	template<class R> inline void GridPavingCursor<R>::set_enabled() const {
		return _theStack[ _currentStackIndex ]->set_enabled();
	}
	
	template<class R> inline void GridPavingCursor<R>::set_disabled() const {
		return _theStack[ _currentStackIndex ]->set_disabled();
	}

	template<class R> inline bool GridPavingCursor<R>::is_left_child() const{
		//If there is a parent node and the given node is it's left child
		return ( _currentStackIndex > 0 ) && ( _theStack[ _currentStackIndex - 1 ]->left_node() == _theStack[ _currentStackIndex ] );
	}
	
	template<class R> inline bool GridPavingCursor<R>::is_right_child() const{
		//If there is a parent node and the given node is it's right child
		return ( _currentStackIndex > 0 ) && ( _theStack[ _currentStackIndex - 1 ]->right_node() == _theStack[ _currentStackIndex ] );
	}

	template<class R> inline GridPavingCursor<R>& GridPavingCursor<R>::move_up() {
		if( ! is_root() ){
			//Remove the current node from the stack
			pop();
			//Recompute the PavingGridCell
			updateTheCurrentGridCell( indeterminate );
			//Return the object back
			return ( * this);
		} else {
			throw NotAllowedMoveException(__PRETTY_FUNCTION__);
		}
	}
	
	template<class R> inline GridPavingCursor<R>& GridPavingCursor<R>::move_left() {
		return move(false);
	}
	
	template<class R> inline GridPavingCursor<R>& GridPavingCursor<R>::move_right() {
		return move(true);
	}
	
	template<class R> inline void GridPavingCursor<R>::updateTheCurrentGridCell( tribool left_or_right ){
		if( ! indeterminate(left_or_right) ){
			//Here left_or_right is either true or false, also true defines moving to the right branch
			_theCurrentGridCell._theWord.push_back( definitely( left_or_right ) );
		} else {
			//If left_or_right is indeterminate, this means that we go "up"
			//in the tree, so we remove the last bit of the path.
			_theCurrentGridCell._theWord.pop_back();
		}
		_theCurrentGridCell = GridPavingCell<R>( _theCurrentGridCell._theGrid, _theCurrentGridCell._theHeight, _theCurrentGridCell._theWord );
	}

	template<class R> inline GridPavingCursor<R>& GridPavingCursor<R>::move(bool left_or_right) {
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
			throw NotAllowedMoveException(__PRETTY_FUNCTION__);
		}
	}

	template<class R> inline const GridPavingCell<R>& GridPavingCursor<R>::cell() const {
		return _theCurrentGridCell;
	}

	template<class R> inline bool GridPavingCursor<R>::operator==(const GridPavingCursor<R>& anotherGridPavingCursor) const {
		bool areEqual = false;
		if( (this->_currentStackIndex >=0) && (anotherGridPavingCursor._currentStackIndex >=0 ) ){
			areEqual = (this->_theStack[this->_currentStackIndex] == anotherGridPavingCursor._theStack[anotherGridPavingCursor._currentStackIndex]);
		}
		return areEqual;
	}

	template<class R> inline GridSubPaving<R> GridPavingCursor<R>::operator*() {
		//IVAN S ZAPREEV:
		//NOTE: The first three parameters define the location of the _theStack[ _currentStackIndex ]
		//node with respect to the primary cell of the GridPaving.
		return GridSubPaving<R>( _theCurrentGridCell._theGrid,
				_theCurrentGridCell._theHeight,
				_theCurrentGridCell._theWord,
				_theStack[ _currentStackIndex ] );
	}
	
	template<class R> inline const GridSubPaving<R> GridPavingCursor<R>::operator*() const {
		return (const GridSubPaving<R> & ) *(* this);
	}

/****************************************GridPavingConstIterator************************************/
	
	template<class R> inline GridPavingConstIterator<R>::GridPavingConstIterator() {
		//IVAN S ZAPREEV:
		// WARNING: This constructor should never be called, this is why it is private
		
		throw NotImplemented(__PRETTY_FUNCTION__);
	}

	template<class R> inline GridPavingConstIterator<R>::GridPavingConstIterator( const GridPavingConstIterator<R>& theGridPavingIter ) : 
			_is_in_end_state(theGridPavingIter._is_in_end_state), _pGridPavingCursor(theGridPavingIter._pGridPavingCursor) {
		//Copy everything
	}

	template<class R> inline GridPavingConstIterator<R>::GridPavingConstIterator( const GridSubPaving<R> * pSubPaving, const tribool firstLastNone ):
											_pGridPavingCursor(pSubPaving) {
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
	
	template<class R> void GridPavingConstIterator<R>::increment() {
		//If we are not done iterating
		if( ! _is_in_end_state){
			//We are at some enabled-leaf node and we want to find the next one.
			//The next node is somewhere on the right in the tree.
			find_next_enabled_leaf();
		}
	}
	
	template<class R> inline GridPavingConstIterator<R>::~GridPavingConstIterator() {
		//IVAN S ZAPREEV:
		//WARNING: There are no memory allocations in this class that have to be cleaned
	}
	
	template<class R> inline bool GridPavingConstIterator<R>::equal( GridPavingConstIterator<R> const & theOtherIterator) const {
		//Check if both iterators are in the "end iterator" state
		bool result = theOtherIterator._is_in_end_state && this->_is_in_end_state;
		
		//If not then check if the cursors are equal (i.e. if they point to the same binary-tree node of the same (subpaving)
		if( ! result ){
			if( theOtherIterator._is_in_end_state || this->_is_in_end_state ){
				//If at one is in the end state and the other is not then the answer is FALSE
				result = false;
			} else {
				result = theOtherIterator._pGridPavingCursor == this->_pGridPavingCursor;
			}
		}
		
		return result;
	}
	
	template<class R> inline GridPavingCell<R> const& GridPavingConstIterator<R>::dereference() const {
		return _pGridPavingCursor.cell();
	}
	
	template<class R> inline GridPavingCursor<R> const& GridPavingConstIterator<R>::cursor() const {
		return _pGridPavingCursor;
	}
 
/*********************************************GridPavingCell***********************************************/
	
	template<class R> inline GridPavingCell<R>::GridPavingCell(const GridPavingCell<R>& theGridPavingCell):
									_theGrid(theGridPavingCell._theGrid),
									_theHeight(theGridPavingCell._theHeight),
									_theWord(theGridPavingCell._theWord),
									_theBox(theGridPavingCell._theBox){
	}
	
	template<class R> inline GridPavingCell<R>::GridPavingCell(const Grid<R>& theGrid, const uint theHeight, const BinaryWord& theWord) :
								_theGrid(theGrid), _theHeight(theHeight), _theWord(theWord),
								_theBox(compute_box(theGrid, theHeight, theWord)) {
		
	}

	template<class R> inline GridPavingCell<R>::GridPavingCell(const Grid<R>& theGrid, const uint theHeight,
								const BinaryWord& theWord, const Box<R>& theBox):
								_theGrid(theGrid), _theHeight(theHeight),
								_theWord(theWord), _theBox(theBox) {
	}
	
	template<class R> inline void GridPavingCell<R>::primary_cell_at_height( const uint theHeight, int & leftBottomCorner, int & rightTopCorner ) {
		if ( theHeight % 2 == 1 ) {
			leftBottomCorner = 2*leftBottomCorner - rightTopCorner;
		} else {
			rightTopCorner   = 2*rightTopCorner - leftBottomCorner;
		}
	}
	
	template<class R> inline const Grid<R>& GridPavingCell<R>::grid() const {
		return _theGrid;
	}
		
	template<class R> inline const uint& GridPavingCell<R>::height() const {
		return _theHeight;
	}
	
	template<class R> inline const BinaryWord& GridPavingCell<R>::word() const {
		return _theWord;
	}
	
	template<class R> inline const Box<R>& GridPavingCell<R>::box() const {
		return _theBox;
	}

	template<class R> inline dimension_type GridPavingCell<R>::dimension() const {
		return _theGrid.dimension();
	}

	template<class R> inline Interval<R> GridPavingCell<R>::operator[](dimension_type i) const {
		return _theBox.interval(i);
	}
	
/********************************************GridSubPaving******************************************/
	
	template<class R> inline GridSubPaving<R>::GridSubPaving( const Grid<R>& theGrid, const uint theHeight,
		const BinaryWord& theWord, BinaryTreeNode * pRootTreeNode ): _theGridPavingCell(theGrid, theHeight, theWord) {
		_pRootTreeNode = pRootTreeNode;
	}

	template<class R> inline GridSubPaving<R>::~GridSubPaving() {
		//IVAN S ZAPREEV:
		//WARNING: This method should have no implementation what so ever
		//All the synamically allocatged data should be destroyed from the
		//corresponding Paving object
	}

	template<class R> inline uint GridSubPaving<R>::compute_number_subdiv( R theWidth, const R theMaxWidth) const{
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
		uint result = 0;
		if ( theWidth > theMaxWidth ){
			result = (uint) ceil( div_approx( log_approx( div_approx( theWidth, theMaxWidth ) ) , log_approx( R(2.0) ) ) );
		}
		return result;
	}

	template<class R> inline GridPavingCell<R> GridSubPaving<R>::cell() const {
		return _theGridPavingCell;
	}

	template<class R> inline void GridSubPaving<R>::mince( const uint theNewDepth ) {
		_pRootTreeNode->mince( theNewDepth );
	}

	template<class R> inline void GridSubPaving<R>::recombine() {
		_pRootTreeNode->recombine();
	}

	template<class R> inline uint GridSubPaving<R>::depth() const {
		return _pRootTreeNode->depth();
	}

	template<class R> inline bool GridSubPaving<R>::operator==(const GridSubPaving<R>& anotherGridSubPaving) const{
		return this->_pRootTreeNode == anotherGridSubPaving._pRootTreeNode;
	}
	
	template<class R> inline typename GridSubPaving<R>::const_iterator GridSubPaving<R>::begin() const {
		return typename GridSubPaving<R>::const_iterator(this, true);
	}

	template<class R> inline typename GridSubPaving<R>::const_iterator GridSubPaving<R>::end() const {
		return typename GridSubPaving<R>::const_iterator(this, indeterminate);
	}

	/*! \brief Stream insertion operator. */
	template<class R> inline std::ostream& operator<<(std::ostream& os, const GridPavingCell<R>& gridPavingCell){
		return os << gridPavingCell.to_string();
	}

	template<class R> inline const BinaryTreeNode * GridSubPaving<R>::binary_tree() const {
		return _pRootTreeNode;
	}


/*********************************************GridPaving*********************************************/
	
	template<class R> inline GridPaving<R>::GridPaving( const Grid<R>& theGrid, const bool enable  ) :
								GridSubPaving<R>( theGrid, 0, BinaryWord(), new BinaryTreeNode( enable ) ){
	}

	template<class R> inline GridPaving<R>::GridPaving( const Grid<R>& theGrid, const uint theHeight, BinaryTreeNode * pRootTreeNode ) : 
								GridSubPaving<R>( theGrid, theHeight, BinaryWord(), pRootTreeNode ){
	}

	template<class R> inline GridPaving<R>::GridPaving( const GridPaving<R> & theGridPaving ) :
							GridSubPaving<R>( theGridPaving._theGridPavingCell.grid(), theGridPaving._theGridPavingCell.height(),
									theGridPaving._theGridPavingCell.word(), new BinaryTreeNode( *theGridPaving._pRootTreeNode )) {
		//Call the super constructor: Create an exact copy of the tree, copy the bounding box
	}

	template<class R> inline GridPaving<R>::GridPaving( const uint theDimension, const bool enable ) :
							GridSubPaving<R>( Grid<R>( theDimension, R(1.0) ), 0, BinaryWord(), new BinaryTreeNode( enable )) {
		//We want a [0,1]x...[0,1] cell in N dimensional space with no sxaling or shift of coordinates:
		//1. Create a new non scaling grid with no shift of the coordinates
		//2. The higth of the primary cell is zero, since is is [0,1]x...[0,1] itself
		//3. The binary word that describes the path from the primary cell to the root
		//   of the tree is empty, because any paving always has a primary cell as a root
		//4. A new disabled binary tree node, gives us the root for the paving tree
	}

	template<class R> inline GridPaving<R>::GridPaving(const Grid<R>& theGrid, const Box<R> & theBoundingBox ) :
								GridSubPaving<R>( theGrid, GridPavingCell<R>::smallest_primary_cell( theGrid.dimension(), theBoundingBox ),
											BinaryWord(), new BinaryTreeNode( false ) ) {
		//1. The main point here is that we have to compute the smallest primary cell that contains theBoundingBox
		//2. This cell is defined by it's height and becomes the root of the GridPaving
		//3. Point 2. implies that the word to the root of GridSubPaving should be set to
		//   empty and we have only one disabled node in the binary tree
	}
	
	template<class R> inline GridPaving<R>::GridPaving( const Grid<R>& theGrid, uint theHeight, const BooleanArray& theTree, const BooleanArray& theEnabledCells ) :
								GridSubPaving<R>( theGrid, theHeight, BinaryWord(), new BinaryTreeNode( theTree, theEnabledCells ) ) {
		//Use the super class constructor and the binary tree constructed from the arrays: theTree and theEnabledCells
	}

	template<class R> inline GridPaving<R>* GridPaving<R>::clone() const {
		return new GridPaving<R>( *this );
	}
	
	template<class R> inline GridPaving<R>::~GridPaving() {
		if( GridSubPaving<R>::_pRootTreeNode != NULL){
			delete GridSubPaving<R>::_pRootTreeNode;
			GridSubPaving<R>::_pRootTreeNode = NULL;
		}
	}
	
	template<class R> inline void GridPaving<R>::adjoin( const GridPavingCell<R>& theCell ) {
		bool has_stopped = false;
		//Align the paving and the cell
		BinaryTreeNode* pBinaryTreeNode = align_with_cell( theCell, true, has_stopped );
		
		//If we are not trying to adjoin something into an enabled sub cell of the paving
		if( ! has_stopped ){
			//Add the enabled cell to the binary tree.
			pBinaryTreeNode->add_enabled( theCell.word() );
		}
	}
	
	template<class R> inline void GridPaving<R>::adjoin( const GridSubPaving<R>& theOtherSubPaving ) {
		bool has_stopped = false;
		//Align the paving and the cell
		BinaryTreeNode* pBinaryTreeNode = align_with_cell( theOtherSubPaving.cell(), true, has_stopped );
		
		//If we are not trying to adjoin something into an enabled sub cell of the paving
		if( ! has_stopped ){
			//Now, the pBinaryTreeNode of this paving corresponds to the primary cell common with theOtherSubPaving.
			//The theOtherSubPaving's root node is defined the path theOtherSubPaving.word() which starts in the
			//node corresponding to the common primary cell.
			pBinaryTreeNode->add_enabled( theOtherSubPaving.binary_tree(), theOtherSubPaving.cell().word() );
		}
	}
	
	template<class R> template<class Set> inline void GridPaving<R>::adjoin_inner_approximation( const Set& theSet, const uint depth ) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}
	
	template<class R> template<class Set> inline void GridPaving<R>::adjoin_lower_approximation( const Set& theSet, const uint depth ) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}

/*************************************FRIENDS OF BinaryTreeNode*************************************/

	inline std::ostream& operator<<(std::ostream & theOstream, const BinaryTreeNode & theBinaryTreeRoot ){
		return theOstream << theBinaryTreeRoot.tree_to_string();
	}

/*************************************FRIENDS OF GridPaving*****************************************/


} // namespace Ariadne

