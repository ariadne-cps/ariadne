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

#include "numeric.h"

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

	inline GridPavingCursor::GridPavingCursor(const GridPavingCursor & theSubPaving) :
						_currentStackIndex(theSubPaving._currentStackIndex), _theStack(theSubPaving._theStack), 
						_pSubPaving(theSubPaving._pSubPaving), _theCurrentGridCell(theSubPaving._theCurrentGridCell){
		//Copy everything
	}

	inline GridPavingCursor::GridPavingCursor(const GridSubPaving * pSubPaving) : 
			_currentStackIndex(-1), _pSubPaving(pSubPaving), _theCurrentGridCell( pSubPaving->cell() ) {
		//Remember that GridSubPaving contains GridPavingCell
		
		//Add the current node to the stack, since we are in it

		//IVAN S ZAPREEV:
		//NOTE: There is no need in allocating _theStack elements,
		//it is all done in the push method
		push(pSubPaving->_pRootTreeNode);
	}

	inline GridPavingCursor::~GridPavingCursor() {
		//IVAN S ZAPREEV:
		//WARNING: The subpaving should not be deallocated here!
		//There are no other data feels that we need to deallocate
		_pSubPaving = NULL;
	}
	
	inline void GridPavingCursor::push( BinaryTreeNode* pLatestNode ){
		//If we are out of free space, increase the array's capacity
          if( static_cast<uint>(_currentStackIndex +1) == ( _theStack.size()  ) ){
			_theStack.resize( _theStack.size() + STACK_SIZE_INCREMENT );
		}
		
		//The index for the tree node to be added
		_currentStackIndex += 1;
		
		//Put the tree node pointer into the stack
		_theStack[_currentStackIndex] = pLatestNode;

	}
	
	inline BinaryTreeNode* GridPavingCursor::pop( ){
		BinaryTreeNode* pLastNode = NULL;
		
		//If there are non-root nodes in the stack
		if( _currentStackIndex > 0 ){
			//Return the stack element at _currentStackIndex,
			//then decrement the current element's index.
			return _theStack[ _currentStackIndex-- ];
		}
		
		return pLastNode;
	}
	
	inline bool GridPavingCursor::is_at_the_root() const {
		//If we are pointing at the zero cell of the array then it
		//means that we are in the root of the tree
		return (_currentStackIndex == 0);
	}

	inline bool GridPavingCursor::is_enabled() const {
		return _theStack[ _currentStackIndex ]->is_enabled();
	}
	
	inline bool GridPavingCursor::is_disabled() const {
		return _theStack[ _currentStackIndex ]->is_disabled();
	}

	inline bool GridPavingCursor::is_leaf() const {
		return _theStack[ _currentStackIndex ]->is_leaf();
	}
	
	inline bool GridPavingCursor::is_root() const {
		return  is_at_the_root();
	}

	inline void GridPavingCursor::set_enabled() const {
		return _theStack[ _currentStackIndex ]->set_enabled();
	}
	
	inline void GridPavingCursor::set_disabled() const {
		return _theStack[ _currentStackIndex ]->set_disabled();
	}

	inline bool GridPavingCursor::is_left_child() const{
		//If there is a parent node and the given node is it's left child
		return ( _currentStackIndex > 0 ) && ( _theStack[ _currentStackIndex - 1 ]->left_node() == _theStack[ _currentStackIndex ] );
	}
	
	inline bool GridPavingCursor::is_right_child() const{
		//If there is a parent node and the given node is it's right child
		return ( _currentStackIndex > 0 ) && ( _theStack[ _currentStackIndex - 1 ]->right_node() == _theStack[ _currentStackIndex ] );
	}

	inline GridPavingCursor& GridPavingCursor::move_up() {
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
	
	inline GridPavingCursor& GridPavingCursor::move_left() {
		return move(false);
	}
	
	inline GridPavingCursor& GridPavingCursor::move_right() {
		return move(true);
	}
	
	inline void GridPavingCursor::updateTheCurrentGridCell( tribool left_or_right ){
		if( ! indeterminate(left_or_right) ){
			//Here left_or_right is either true or false, also true defines moving to the right branch
			_theCurrentGridCell._theWord.push_back( definitely( left_or_right ) );
		} else {
			//If left_or_right is indeterminate, this means that we go "up"
			//in the tree, so we remove the last bit of the path.
			_theCurrentGridCell._theWord.pop_back();
		}
		_theCurrentGridCell = GridPavingCell( _theCurrentGridCell._theGrid, _theCurrentGridCell._theHeight, _theCurrentGridCell._theWord );
	}

	inline GridPavingCursor& GridPavingCursor::move(bool left_or_right) {
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

	inline const GridPavingCell& GridPavingCursor::cell() const {
		return _theCurrentGridCell;
	}

	inline bool GridPavingCursor::operator==(const GridPavingCursor& anotherGridPavingCursor) const {
		bool areEqual = false;
		if( (this->_currentStackIndex >=0) && (anotherGridPavingCursor._currentStackIndex >=0 ) ){
			areEqual = (this->_theStack[this->_currentStackIndex] == anotherGridPavingCursor._theStack[anotherGridPavingCursor._currentStackIndex]);
		}
		return areEqual;
	}

	inline GridSubPaving GridPavingCursor::operator*() {
		//IVAN S ZAPREEV:
		//NOTE: The first three parameters define the location of the _theStack[ _currentStackIndex ]
		//node with respect to the primary cell of the GridPaving.
		return GridSubPaving( _theCurrentGridCell._theGrid,
				_theCurrentGridCell._theHeight,
				_theCurrentGridCell._theWord,
				_theStack[ _currentStackIndex ] );
	}
	
	inline const GridSubPaving GridPavingCursor::operator*() const {
		return (const GridSubPaving & ) *(* this);
	}

/****************************************GridPavingConstIterator************************************/
	
	inline GridPavingConstIterator::GridPavingConstIterator( const GridPavingConstIterator& theGridPavingIter ) : 
			_is_in_end_state(theGridPavingIter._is_in_end_state), _pGridPavingCursor(theGridPavingIter._pGridPavingCursor) {
		//Copy everything
	}

	inline GridPavingConstIterator::GridPavingConstIterator( const GridSubPaving * pSubPaving, const tribool firstLastNone ):
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
	
	void GridPavingConstIterator::increment() {
		//If we are not done iterating
		if( ! _is_in_end_state){
			//We are at some enabled-leaf node and we want to find the next one.
			//The next node is somewhere on the right in the tree.
			find_next_enabled_leaf();
		}
	}
	
	inline GridPavingConstIterator::~GridPavingConstIterator() {
		//IVAN S ZAPREEV:
		//WARNING: There are no memory allocations in this class that have to be cleaned
	}
	
	inline bool GridPavingConstIterator::equal( GridPavingConstIterator const & theOtherIterator) const {
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
	
	inline GridPavingCell const& GridPavingConstIterator::dereference() const {
		return _pGridPavingCursor.cell();
	}
	
	inline GridPavingCursor const& GridPavingConstIterator::cursor() const {
		return _pGridPavingCursor;
	}
 
/*********************************************GridPavingCell***********************************************/
	
	inline GridPavingCell::GridPavingCell(const GridPavingCell& theGridPavingCell):
									_theGrid(theGridPavingCell._theGrid),
									_theHeight(theGridPavingCell._theHeight),
									_theWord(theGridPavingCell._theWord),
									_theBox(theGridPavingCell._theBox){
	}
	
	inline GridPavingCell::GridPavingCell(const Grid& theGrid, const uint theHeight, const BinaryWord& theWord) :
								_theGrid(theGrid), _theHeight(theHeight), _theWord(theWord),
								_theBox(compute_box(theGrid, theHeight, theWord)) {
		
	}

	inline GridPavingCell::GridPavingCell(const Grid& theGrid, const uint theHeight,
								const BinaryWord& theWord, const Box& theBox):
								_theGrid(theGrid), _theHeight(theHeight),
								_theWord(theWord), _theBox(theBox) {
	}
	
	inline void GridPavingCell::primary_cell_at_height( const uint theHeight, int & leftBottomCorner, int & rightTopCorner ) {
		if ( theHeight % 2 == 1 ) {
			leftBottomCorner = 2*leftBottomCorner - rightTopCorner;
		} else {
			rightTopCorner   = 2*rightTopCorner - leftBottomCorner;
		}
	}
	
	inline const Grid& GridPavingCell::grid() const {
		return _theGrid;
	}
		
	inline const uint& GridPavingCell::height() const {
		return _theHeight;
	}
	
	inline const BinaryWord& GridPavingCell::word() const {
		return _theWord;
	}
	
	inline const Box& GridPavingCell::box() const {
		return _theBox;
	}

	inline dimension_type GridPavingCell::dimension() const {
		return _theGrid.dimension();
	}

	inline Interval GridPavingCell::operator[](dimension_type i) const {
                return _theBox[i];
	}
	
/********************************************GridSubPaving******************************************/
	
	inline GridSubPaving::GridSubPaving( const Grid& theGrid, const uint theHeight,
		const BinaryWord& theWord, BinaryTreeNode * pRootTreeNode ): _theGridPavingCell(theGrid, theHeight, theWord) {
		_pRootTreeNode = pRootTreeNode;
	}

	inline GridSubPaving::~GridSubPaving() {
		//IVAN S ZAPREEV:
		//WARNING: This method should have no implementation what so ever
		//All the synamically allocatged data should be destroyed from the
		//corresponding Paving object
	}

	inline uint GridSubPaving::compute_number_subdiv( Float theWidth, const Float theMaxWidth) const{
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

	inline GridPavingCell GridSubPaving::cell() const {
		return _theGridPavingCell;
	}

	inline void GridSubPaving::mince( const uint theNewDepth ) {
		_pRootTreeNode->mince( theNewDepth );
	}

	inline void GridSubPaving::recombine() {
		_pRootTreeNode->recombine();
	}

	inline uint GridSubPaving::depth() const {
		return _pRootTreeNode->depth();
	}

	inline bool GridSubPaving::operator==(const GridSubPaving& anotherGridSubPaving) const{
		return this->_pRootTreeNode == anotherGridSubPaving._pRootTreeNode;
	}
	
	inline GridSubPaving::const_iterator GridSubPaving::begin() const {
		return GridSubPaving::const_iterator(this, true);
	}

	inline GridSubPaving::const_iterator GridSubPaving::end() const {
		return GridSubPaving::const_iterator(this, indeterminate);
	}

	/*! \brief Stream insertion operator. */
	inline std::ostream& operator<<(std::ostream& os, const GridPavingCell& gridPavingCell){
		return os << gridPavingCell.to_string();
	}

	inline const BinaryTreeNode * GridSubPaving::binary_tree() const {
		return _pRootTreeNode;
	}


/*********************************************GridPaving*********************************************/
	
	inline GridPaving::GridPaving( const Grid& theGrid, const bool enable  ) :
								GridSubPaving( theGrid, 0, BinaryWord(), new BinaryTreeNode( enable ) ){
	}

	inline GridPaving::GridPaving( const Grid& theGrid, const uint theHeight, BinaryTreeNode * pRootTreeNode ) : 
								GridSubPaving( theGrid, theHeight, BinaryWord(), pRootTreeNode ){
	}

	inline GridPaving::GridPaving( const GridPaving & theGridPaving ) :
							GridSubPaving( theGridPaving._theGridPavingCell.grid(), theGridPaving._theGridPavingCell.height(),
									theGridPaving._theGridPavingCell.word(), new BinaryTreeNode( *theGridPaving._pRootTreeNode )) {
		//Call the super constructor: Create an exact copy of the tree, copy the bounding box
	}

	inline GridPaving::GridPaving( const uint theDimension, const bool enable ) :
							GridSubPaving( Grid( theDimension, Float(1.0) ), 0, BinaryWord(), new BinaryTreeNode( enable )) {
		//We want a [0,1]x...[0,1] cell in N dimensional space with no sxaling or shift of coordinates:
		//1. Create a new non scaling grid with no shift of the coordinates
		//2. The higth of the primary cell is zero, since is is [0,1]x...[0,1] itself
		//3. The binary word that describes the path from the primary cell to the root
		//   of the tree is empty, because any paving always has a primary cell as a root
		//4. A new disabled binary tree node, gives us the root for the paving tree
	}

	inline GridPaving::GridPaving(const Grid& theGrid, const Box & theBoundingBox ) :
								GridSubPaving( theGrid, GridPavingCell::smallest_primary_cell( theGrid.dimension(), theBoundingBox ),
											BinaryWord(), new BinaryTreeNode( false ) ) {
		//1. The main point here is that we have to compute the smallest primary cell that contains theBoundingBox
		//2. This cell is defined by it's height and becomes the root of the GridPaving
		//3. Point 2. implies that the word to the root of GridSubPaving should be set to
		//   empty and we have only one disabled node in the binary tree
	}
	
	inline GridPaving::GridPaving( const Grid& theGrid, uint theHeight, const BooleanArray& theTree, const BooleanArray& theEnabledCells ) :
								GridSubPaving( theGrid, theHeight, BinaryWord(), new BinaryTreeNode( theTree, theEnabledCells ) ) {
		//Use the super class constructor and the binary tree constructed from the arrays: theTree and theEnabledCells
	}

	inline GridPaving* GridPaving::clone() const {
		return new GridPaving( *this );
	}
	
	inline GridPaving::~GridPaving() {
		if( GridSubPaving::_pRootTreeNode != NULL){
			delete GridSubPaving::_pRootTreeNode;
			GridSubPaving::_pRootTreeNode = NULL;
		}
	}
	
	inline void GridPaving::adjoin( const GridPavingCell& theCell ) {
		bool has_stopped = false;
		//Align the paving and the cell
		BinaryTreeNode* pBinaryTreeNode = align_with_cell( theCell, true, has_stopped );
		
		//If we are not trying to adjoin something into an enabled sub cell of the paving
		if( ! has_stopped ){
			//Add the enabled cell to the binary tree.
			pBinaryTreeNode->add_enabled( theCell.word() );
		}
	}
	
	inline void GridPaving::adjoin( const GridSubPaving& theOtherSubPaving ) {
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
	
	template<class Set> inline void GridPaving::adjoin_inner_approximation( const Set& theSet, const uint depth ) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}
	
	template<class Set> inline void GridPaving::adjoin_lower_approximation( const Set& theSet, const uint depth ) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}

/*************************************FRIENDS OF BinaryTreeNode*************************************/

	inline std::ostream& operator<<(std::ostream & theOstream, const BinaryTreeNode & theBinaryTreeRoot ){
		return theOstream << theBinaryTreeRoot.tree_to_string();
	}

/*************************************FRIENDS OF GridPaving*****************************************/


} // namespace Ariadne

