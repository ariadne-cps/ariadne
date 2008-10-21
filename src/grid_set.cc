/***************************************************************************
 *            grid_set.cc
 *
 *  Copyright  2008  Ivan S. Zapreev, Pieter Collins
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

#include <iostream>
#include <iomanip>

#include "macros.h"
#include "exceptions.h"
#include "stlio.h"
#include "list_set.h"
#include "grid_set.h"

#include "set_interface.h"

namespace Ariadne {

typedef size_t size_type;


/****************************************Grid**********************************************/

 
Grid::~Grid()
{
}



 
Grid::Grid()
  : _data(new Data())
{
}


 
Grid::Grid(const Grid& gr)
  : _data(gr._data)
{
}


 
Grid::Grid(uint d, Float l)
  : _data(new Data())
{
  array<Float> origin(d,Float(0));
  array<Float> lengths(d,l);
  this->_create(origin,lengths);
}

 

Grid::Grid(const Vector<Float>& lengths)
  : _data(new Data())
{
  array<Float> o(lengths.size(),0);
  const array<Float>& l=reinterpret_cast<const array<Float>&>(lengths.data());
  this->_create(o,l);
}


 
Grid::Grid(const Vector<Float>& origin, const Vector<Float>& lengths)
  : _data(new Data())
{
  if(origin.size() != lengths.size()) {
    throw IncompatibleSizes(__PRETTY_FUNCTION__);
  }
  const array<Float>& o=reinterpret_cast<const array<Float>&>(origin.data());
  const array<Float>& l=reinterpret_cast<const array<Float>&>(lengths.data());
  this->_create(o,l);
}


 



void
Grid::_create(const array<Float>& origin, const array<Float>& lengths) 
{
  this->_data->_origin=origin;
  this->_data->_lengths=lengths;
}




uint
Grid::dimension() const
{
  return this->_data->_lengths.size();
}



const Vector<Float>&
Grid::origin() const
{
  return reinterpret_cast<const Vector<Float>&>(this->_data->_origin);
}



const Vector<Float>&
Grid::lengths() const
{
  return reinterpret_cast<const Vector<Float>&>(this->_data->_lengths);
}





Float
Grid::coordinate(uint d, dyadic_type x) const 
{
  return add_approx(this->_data->_origin[d],mul_approx(this->_data->_lengths[d],x));
}


Float
Grid::subdivision_coordinate(uint d, dyadic_type x) const 
{
  return add_approx(this->_data->_origin[d],mul_approx(this->_data->_lengths[d],x));
}


Float
Grid::subdivision_coordinate(uint d, integer_type n) const 
{
  return add_approx(this->_data->_origin[d],mul_approx(this->_data->_lengths[d],n));
}


 
int
Grid::subdivision_index(uint d, const real_type& x) const 
{
  
  
  Float half=0.5;
  int n=floor(add_approx(div_approx(sub_approx(x,this->_data->_origin[d]),this->_data->_lengths[d]),half));
  Float sc=add_approx(this->_data->_origin[d],mul_approx(this->_data->_lengths[d],n));
    if(sc == x) { 
      return n; 
    } else {
      std::cerr << std::setprecision(20) << std::boolalpha
                << "sc=" << sc << " x=" << x << " sc-x=" << Interval(sc-x) << "\n"
                << "sc==x=" << (sc==x) << " sc!=x=" << (sc!=x)
                << " sc<x=" << (sc<x) << " sc>x=" << (sc>x) << " sc<=x=" << (sc<=x) << " sc>=x=" << (sc>=x) << std::endl; 
      ARIADNE_THROW(InvalidGridPosition,std::setprecision(20)<<"Grid::subdivision_index(uint d,real_type x)","d="<<d<<", x="<<x<<", this->origin[d]="<<this->_data->_origin[d]<<", this->lengths[d]="<<this->_data->_lengths[d]<<" (closest value is "<<sc<<")");
    }
}


 
int
Grid::subdivision_lower_index(uint d, const real_type& x) const 
{
  
  
  int n=floor(div_down(sub_down(x,this->_data->_origin[d]),this->_data->_lengths[d]));
  if(x>=add_approx(this->_data->_origin[d],mul_approx(this->_data->_lengths[d],(n+1)))) {
    return n+1;
  } else {
    return n;
  }
}


 
int
Grid::subdivision_upper_index(uint d, const real_type& x) const 
{
  
  
  int n=ceil(div_up(sub_up(x,this->_data->_origin[d]),this->_data->_lengths[d]));
  if(x<=add_approx(this->_data->_origin[d],mul_approx(this->_data->_lengths[d],(n-1)))) {
    return n-1;
  } else {
    return n;
  }
}



 
bool 
Grid::operator==(const Grid& g) const
{
  if(this->_data==g._data) { 
    return true; 
  } else {
    return this->_data->_origin==g._data->_origin && this->_data->_lengths==g._data->_lengths;
  }
}


 
bool 
Grid::operator!=(const Grid& g) const
{
  return !(*this==g);
}



array<double>
Grid::index(const Vector<Float>& pt) const
{
        array<double> res(pt.dimension());
        for(size_t i=0; i!=res.size(); ++i) {
                res[i]=subdivision_index(i,pt[i]);
        }
        return res;
}



array<double>
Grid::lower_index(const Vector<Interval>& bx) const {
        array<double> res(bx.dimension());
        for(size_t i=0; i!=res.size(); ++i) {
                res[i]=subdivision_lower_index(i,bx[i].lower());
        }
        return res;
}



array<double>
Grid::upper_index(const Vector<Interval>& bx) const {
        array<double> res(bx.size());
        for(size_type i=0; i!=res.size(); ++i) {
                res[i]=subdivision_upper_index(i,bx[i].upper());
        }
        return res;
}





Vector<Float>
Grid::point(const array<int>& a) const
{
  Vector<float> res(a.size());
  for(size_type i=0; i!=res.dimension(); ++i) {
          res[i]=this->_data->_origin[i]+this->_data->_lengths[i]*a[i];
  }
  return res;
}


Vector<Float>
Grid::point(const array<double>& a) const
{
  Vector<float> res(a.size());
  for(size_type i=0; i!=res.dimension(); ++i) {
          res[i]=this->_data->_origin[i]+this->_data->_lengths[i]*a[i];
  }
  return res;
}



Vector<Interval>
Grid::box(const array<double>& lower, const array<double>& upper) const
{
        Vector<Interval> res(lower.size());
        for(size_type i=0; i!=res.size(); ++i) {
                res[i]=Interval(this->subdivision_coordinate(i,lower[i]),
                                this->subdivision_coordinate(i,upper[i]));
        }
        return res;
}



std::ostream& 
operator<<(std::ostream& os, const Grid& gr) 
{
        os << "Grid( ";
        os << "origin=" << gr.origin() << ", ";
        os << "lengths=" << gr.lengths() << " )";
        return os;
}





/****************************************BinaryTreeNode**********************************************/

	void BinaryTreeNode::restore_node( BinaryTreeNode * theCurrentNode, uint & arr_index, uint & leaf_counter,
					const BooleanArray& theTree, const BooleanArray& theEnabledCells) {
		//If we are not done with the the tree yet
		if( arr_index < theTree.size() ) {
			//If we are in a non-leaf node then go further
			if( theTree[arr_index] ) {
					theCurrentNode->split();
					//IVAN S. ZAPREEV:
					//NOTE: We assume a correct input, i.e. both children are present
					//NOTE: We increase the arr_index before calling recursion for each
					//      of the subnodes of the tree
					restore_node( theCurrentNode->left_node(), ++arr_index, leaf_counter, theTree, theEnabledCells );
					restore_node( theCurrentNode->right_node(), ++arr_index, leaf_counter, theTree, theEnabledCells );
			} else {
				//If we are in a leaf node then chek if it needs to be enabled/disabled
				if( theEnabledCells[leaf_counter] ) {
					theCurrentNode->set_enabled();
				} else {
					theCurrentNode->set_disabled();
				}
				leaf_counter++;
			}
		}
	}

	void BinaryTreeNode::mince_node(BinaryTreeNode * pCurrentNode, const uint depth) {
		//If we need to mince further.
		if( depth > 0 ){
			//If the node is present, i.e. the parent is not a disabled node.
			if( pCurrentNode != NULL ){
				//If the current node is not disabled: enabled (leaf) or
				//indeterminate (non-leaf) then there is a work to do.
				if( ! pCurrentNode->is_disabled() ){
					pCurrentNode->split();
				
					const int remaining_depth = depth - 1;
				
					mince_node( pCurrentNode->_pLeftNode, remaining_depth );
					mince_node( pCurrentNode->_pRightNode, remaining_depth );
				}
			} else {
				throw std::runtime_error(__PRETTY_FUNCTION__);
			}
		}
	}

	void BinaryTreeNode::recombine_node(BinaryTreeNode * pCurrentNode) {
		//Just a safety check
		if( pCurrentNode != NULL ){
			//If it is not a leaf node then it should have both of it's subnodes != NULL
			if( ! pCurrentNode->is_leaf() ){
				BinaryTreeNode * pLeftNode = pCurrentNode->_pLeftNode;
				BinaryTreeNode * pRightNode = pCurrentNode->_pRightNode;
			
				//This recursive calls ensure that we do recombination from the bottom up
				recombine_node( pLeftNode );
				recombine_node( pRightNode );

				//Do the recombination for the leaf nodes rooted to pCurrentNode
				if( pLeftNode->is_leaf() && pRightNode->is_leaf() ){
					if( pLeftNode->_isEnabled == pRightNode->_isEnabled ){
						//Make it the leaf node with the derived _isEnabled value
						pCurrentNode->make_leaf( pLeftNode->_isEnabled );
					}
				}
			}
		} else {
			throw std::runtime_error(__PRETTY_FUNCTION__);
		}
	}
	
	uint BinaryTreeNode::depth() const {
		//The depth of the sub-tree rooted to this node
		uint result;
		
		if( ! this->is_leaf() ){
			//If the node is not a leaf, compute the depth of the sub-trees and take the maximum + 1
			//Note that, both left and right sub-nodes must exist, by to the way we construct the tree
			result = max(this->left_node()->depth(), this->right_node()->depth() ) + 1;
		} else {
			//If the node is a leaf then the depth of the sub-tree is zero
			result = 0;
		}
		
		return result;
	}

	string BinaryTreeNode::tree_to_string() const {
		string result = "";
		char trib_value;

		//Get the node's _isEnabled value
		if( is_enabled() ){
			trib_value = '+';
		} else {
			if( is_disabled() ){
				trib_value = '-';
			} else {
				trib_value = '?';
			}
		}
		result += trib_value;
		
		//Indicate whether it is a leaf or not
		if( is_leaf() ){
			result += "0";
		} else {
			result += "1";
			result += _pLeftNode->tree_to_string();
			result += "1";
			result += _pRightNode->tree_to_string();
		}
		return result;
	}
	
	void BinaryTreeNode::add_enabled( const BinaryTreeNode * pOtherSubTree, const BinaryWord & path ){
		//1. Locate the node, follow the path until it's end or until we meet an enabled node
		BinaryTreeNode * pCurrentSubTree = this;
		uint position = 0;
		while( ( position < path.size() ) && ! pCurrentSubTree->is_enabled() ){
			//Split the node, if it is not a leaf it will not be changed
			pCurrentSubTree->split();
			//Follow the path step
			pCurrentSubTree = ( path[position] ? pCurrentSubTree->_pRightNode : pCurrentSubTree->_pLeftNode );
			//Go to the next path element
			position ++;
		}
		//2. Now we are in the right node of this tree or we have met an enabled node on the path,
		//   Thus, if this node is not enabled then we go on with adding subTree.
		if( ! pCurrentSubTree->is_enabled() ){
			add_enabled( pCurrentSubTree, pOtherSubTree );
		}
	}
	
	void BinaryTreeNode::add_enabled( BinaryTreeNode* pToTreeRoot, const BinaryTreeNode* pFromTreeRoot ){
		if( pToTreeRoot->is_leaf() ){
			//If we are adding something to a leaf node
			if( pToTreeRoot->is_enabled() ){
				//Do nothing, adding to an enabled leaf node (nothing new can be added)
			} else {
				if( pFromTreeRoot->is_leaf() ){
					//If we are adding something to a disabled leaf node
					if( pFromTreeRoot->is_enabled() ){
						//Adding an enabled node: Enable the node of pToTreeRoot
						pToTreeRoot->set_enabled();
					} else {
						//Do nothing, adding a disabled leaf node to a disabled leaf node
					}
				} else {
					//Adding a subtree pFromTreeRoot to a disabled leaf node pToTreeRoot
					//Using copy constructors here, to avoid memory collisions
					pToTreeRoot->_pLeftNode = new BinaryTreeNode( *pFromTreeRoot->_pLeftNode );
					pToTreeRoot->_pRightNode = new BinaryTreeNode( *pFromTreeRoot->_pRightNode );
					//Set the leaf node as unknown, since we do not know what is below
					pToTreeRoot->set_unknown();
				}
			}
		} else {
			//If we are adding something to a non-leaf node
			if( pFromTreeRoot->is_leaf() ){
				//Adding a leaf to a non-leaf node
				if( pFromTreeRoot->is_enabled() ){
					//Make the enabled leaf node
					pToTreeRoot->make_leaf(true);
				} else {
					//Do nothing, adding a disabled node to a sub tree (nothing new can be added)
				}
			} else {
				//Adding a non-leaf node to a non-leaf node, do recursion
				add_enabled( pToTreeRoot->_pLeftNode, pFromTreeRoot->_pLeftNode );
				add_enabled( pToTreeRoot->_pRightNode, pFromTreeRoot->_pRightNode );
			}
		}
	}

	void BinaryTreeNode::add_enabled( BinaryTreeNode* pRootTreeNode, const BinaryWord& path, const uint position ) {
		if( position < path.size() ) {
			//There is still something to do
			if( pRootTreeNode->is_leaf() ){
				if( pRootTreeNode->is_enabled() ) {
					//This leaf is enabled so adding path will not change anything
					return;
				} else {
					//Split the disabled node
					pRootTreeNode->split();
				}
			}
			//Go left-right depenting on the specified path
			add_enabled( ( path[position] ? pRootTreeNode->_pRightNode : pRootTreeNode->_pLeftNode ), path, position + 1 );
		} else {
			//We are at the destination node
			if( pRootTreeNode->is_leaf() ){
				//Mark the node as enabled
				pRootTreeNode->set_enabled();
			} else {
				//If this is not a leaf node, then make it leaf
				//The leafs below are not interesting any more
				pRootTreeNode->make_leaf(true);
			}
		}
	}
	
	BinaryTreeNode * BinaryTreeNode::prepend_tree( const BinaryWord & rootNodePath, BinaryTreeNode * oldRootNode){
		//Create the new binary tree node
		BinaryTreeNode * pRootBinaryTreeNode = new BinaryTreeNode(), * pCurrentBinaryTreeNode = pRootBinaryTreeNode;
		uint i = 0;
		//Loop until the last path element, because it has to be treated in a different manner
		for( ; i < ( rootNodePath.size() - 1 ) ; i++ ){
			//Split the node
			pCurrentBinaryTreeNode->split();
			//Move to the appropriate subnode
			pCurrentBinaryTreeNode = (rootNodePath[i]) ? pCurrentBinaryTreeNode->_pRightNode : pCurrentBinaryTreeNode->_pLeftNode;
		}
		//Split the node for the last time
		pCurrentBinaryTreeNode->split();
		//Substitute the new primary cell with the one we had before
		if( rootNodePath[i] ){
			delete pCurrentBinaryTreeNode->_pRightNode;
			pCurrentBinaryTreeNode->_pRightNode = oldRootNode;
		} else {
			delete pCurrentBinaryTreeNode->_pLeftNode;
			pCurrentBinaryTreeNode->_pLeftNode = oldRootNode;
		}
		return pRootBinaryTreeNode;
	}

/********************************************GridTreeCursor***************************************/
	
	string GridTreeCursor::to_string() const{
		stringstream tmp_stream;
		string result = "\n The underlying subpaving:"+_pSubPaving->to_string();

		tmp_stream << _currentStackIndex;
		result += "\n The current stack index: " + tmp_stream.str();
		tmp_stream.str("");

		result += "\n The current stack data: \n";
		for(int i = 0; i <= _currentStackIndex; i++){
			tmp_stream << i;
			result += " Element " + tmp_stream.str() + ": "+ _theStack[i]->node_to_string()+"\n";
			tmp_stream.str("");
		}
		result += " The current grid cell: " + _theCurrentGridCell.to_string();
		return result;
	}


/****************************************GridTreeConstIterator************************************/

	void GridTreeConstIterator::find_next_enabled_leaf() {
		if( _pGridTreeCursor.is_left_child() ) {
			//Move to the parent node
			_pGridTreeCursor.move_up();
			//The right node must exist, due to the way we allocate the tree
			//Also, this node is the root of the branch which we did not investigate.
			_pGridTreeCursor.move_right();
			//Find the first enabled leaf on this sub-tree
			if( ! navigate_to(true) ){
				//If there are no enabled leafs in the subtree, then
				//move up and check the remaining branches recursively.
				find_next_enabled_leaf();
			}
		} else {
			if( _pGridTreeCursor.is_right_child() ) {
				//Move to the parent node
				_pGridTreeCursor.move_up();
				//Move up and check the remaining branches recursively.
				find_next_enabled_leaf();
			} else {
				//Is the root node already and we've seen all the leaf nodes
				_is_in_end_state = true;
			}
		}
	}

	bool GridTreeConstIterator::navigate_to(bool firstLast){
		bool isEnabledLeafFound = false;
		if( _pGridTreeCursor.is_leaf() ){
			if ( _pGridTreeCursor.is_enabled() ) {
				//If the leaf is enabled, then the search is over
				isEnabledLeafFound = true;
			}
		} else {
			//If it is not a leaf then we need to keep searching
			if( firstLast ) {
				//If we are looking for the first enabled
				//node, then go to the left sub node
				_pGridTreeCursor.move_left();
			} else {
				//Otherwise we go to the right
				_pGridTreeCursor.move_right();
			}
			//Do recursive check of the newly visited node
			isEnabledLeafFound = navigate_to( firstLast );
			//If the leaf node is not found yet
			if ( ! isEnabledLeafFound ) {
				if( firstLast ) {
					_pGridTreeCursor.move_right();
				} else {
					_pGridTreeCursor.move_left();
				}
				isEnabledLeafFound = navigate_to( firstLast );
			}
		}
		
		//If the enabled leaf is not found and we are not in the root then we go back
		if( ! isEnabledLeafFound && ! _pGridTreeCursor.is_root() ){
			_pGridTreeCursor.move_up();
		}
		
		return isEnabledLeafFound;
	}


/*********************************************GridCell***********************************************/

	Box GridCell::primary_cell( const uint theHeight, const dimension_type dimensions ) {
		int leftBottomCorner = 0, rightTopCorner = 1;
		//The zero level coordinates are known, so we need to iterate only for higher level primary cells
		for(uint i = 1; i <= theHeight; i++){
			primary_cell_at_height(i, leftBottomCorner, rightTopCorner);
		}
		// 1.2 Constructing and return the box defining the primary cell (relative to the grid).
		
		return Box( Vector< Interval >( dimensions, Interval( leftBottomCorner, rightTopCorner ) ) );
	}

	uint GridCell::smallest_primary_cell( const dimension_type dimensions, const Box& theBox ) {
		int leftBottomCorner = 0, rightTopCorner = 1;
		uint height = 0;
		//The zero level coordinates are known, so we need to iterate only for higher level primary cells
		do{
			//Check if the given box is a subset of a primary cell
			Box primaryCellBox( Vector< Interval >( dimensions, Interval( leftBottomCorner, rightTopCorner ) ) );
			if( definitely( subset(theBox, primaryCellBox ) ) ) {
				//If yes then we are done
				break;
			}
			//Otherwise increase the height and recompute the new borders
			primary_cell_at_height( ++height, leftBottomCorner, rightTopCorner);
		}while(true);
		
		return height;
	}

	/*! \brief Apply grid data \a theGrid to \a theGridBox in order to compute the box dimensions in the original space*/
	Box GridCell::grid_box_to_space(const Box & theGridBox, const Grid& theGrid ){
		const uint dimensions = theGrid.dimension();
		Box theTmpBox( dimensions );
		
		Point theGridOrigin( theGrid.origin() );
		Vector<Float> theGridLengths( theGrid.lengths() );
		
		for(uint current_dimension = 0; current_dimension < dimensions; current_dimension ++){
			const Float theDimLength = theGridLengths[current_dimension];
			const Float theDimOrigin = theGridOrigin[current_dimension];
			//Recompute the new dimension coordinates, detaching them from the grid 
			theTmpBox[current_dimension].set_lower( add_approx( theDimOrigin, mul_approx( theDimLength, theGridBox[current_dimension].lower() ) ) );
                        theTmpBox[current_dimension].set_upper( add_approx( theDimOrigin, mul_approx( theDimLength, theGridBox[current_dimension].upper() ) ) );
		}
		
		return theTmpBox;
	}

	//The box is not related to the Grid, whereas the binary tree is related to the Lattice of the grid:
	// 1. Compute the primary cell located the the height \a theHeight above the zero level,
	// 2. Compute the cell defined by the path \a theWord (from the primary cell).
	// 3. Use Grid data to compute the box coordinates in the original space.
	Box GridCell::compute_box(const Grid& theGrid, const uint theHeight, const BinaryWord& theWord){
		//1. Obtain the primary-cell box, related to some grid
		const uint dimensions = theGrid.dimension();
		Box  theTmpBox( primary_cell( theHeight , dimensions ) );

		//2. Compute the cell on some grid, corresponding to the binary path from the primary cell.
		uint current_dimension = 0;
		for(uint i = 0; i < theWord.size(); i++){
			//We move through the dimensions in a linear fasshion
			current_dimension = i % dimensions;
			//Compute the middle point of the box's projection onto
			//the dimension \a current_dimension (relative to the grid)
			Float middlePointInCurrDim = theTmpBox[current_dimension].midpoint();
			if( theWord[i] ){
				//Choose the right half
				theTmpBox[current_dimension].set_lower( middlePointInCurrDim );
			} else {
				//Choose the left half
				theTmpBox[current_dimension].set_upper( middlePointInCurrDim );
			}
		}
		
		// 3. Use Grid data to compute the box coordinates in the original space.
		return grid_box_to_space( theTmpBox, theGrid );
	}
	
	 BinaryWord GridCell::primary_cell_path( const uint dimensions, const uint topPCellHeight, const uint bottomPCellHeight) {
		BinaryWord theBinaryPath;
		
		//The path from one primary cell to another consists of alternating subsequences
		//of length \a dimensions. These subsequences consist either of ones or zeroes.
		//Odd primary cell height means that the first subsequence will consist of all
		//ones. Even primary cell height indicates the the first subsequence will consis
		//of all zeroes. This is due to the way we do the space subdivisions.
		if( topPCellHeight > bottomPCellHeight ){
			for( uint i = topPCellHeight; i > bottomPCellHeight; i-- ){
				bool odd_height = (i % 2) != 0;
				for( uint j = 0; j < dimensions; j++ ){
					theBinaryPath.push_back( odd_height );
				}
			}
		}
		
		return theBinaryPath;
	}
		
	/*! \brief Serialize the current object state, this is mostly needed for testing */
	string GridCell::to_string() const{
		//Write the grid data to the string stream
		stringstream tmp_stream;
		tmp_stream << _theGrid;
		
		//Get the Grid data
		string result = "\n The grid: " + tmp_stream.str();
		tmp_stream.str("");
		
		//Get the Primary cell higth
		tmp_stream << _theHeight;
		result += "\n Primary root cell's higth: " + tmp_stream.str();
		tmp_stream.str("");
		
		//Get the binary word path to the subpaving root cell
		tmp_stream << _theWord;
		result += "\n The path to the cell from the primary root cell: " + tmp_stream.str();
		tmp_stream.str("");
		
		//Print the Box corresponding to this cell
		tmp_stream << _theBox;
		result += "\n The cell's box: " + tmp_stream.str();
		tmp_stream.str("");
		
		return result;
	}




/********************************************GridTreeSubset*****************************************/


GridCellListSet::~GridCellListSet()
{
}


GridCellListSet::GridCellListSet(const Grid& g)
  : _grid(g), _list()
{
}




GridCellListSet::GridCellListSet(const GridTreeSet& gts)
  : _grid(gts.grid()), _list(gts.begin(),gts.end())
{
}






GridCellListSet::operator ListSet<Box>() const
{
  ListSet<Box> result(this->dimension());
  for(size_type i=0; i!=size(); ++i) {
    result.push_back((*this)[i].box());
  }
  return result;
}


const Grid& 
GridCellListSet::grid() const 
{
  return this->_grid; 
}


uint
GridCellListSet::dimension() const 
{
  return this->_grid.dimension(); 
}


uint
GridCellListSet::size() const 
{
  return this->_list.size(); 
}


GridCell
GridCellListSet::operator[](uint i) const 
{
  return this->_list[i]; 
}


GridCell
GridCellListSet::pop() 
{
  GridCell back=this->_list.back(); 
  this->_list.pop_back();
  return back;
}


void
GridCellListSet::clear()
{
  this->_list.clear();
}


std::ostream&     
operator<<(std::ostream& os, const GridCellListSet& gcls)
{
  os << "GridCellListSet("
     << " grid=" << gcls.grid() << ","
     << " size=" << gcls.size() << ","
     << " cells=" << gcls._list
     << " )";
  return os;
}





/********************************************GridTreeSubset*****************************************/

	void GridTreeSubset::subdivide( Float theMaxCellWidth ) {
		//1. Take the Box of this GridTreeSubset's GridCell
		//   I.e. the box that corresponds to the root cell of
		//   the GridTreeSubset in the original space.
		const Box& theRootCellBox = _theGridCell.box();
		
		//2. Compute the widths of the box in each dimension and the maximum number
		//   among the number of subdivisions that we need to do in each dimension
		//   in order to make the width in this dimension <= theMaxCellWidth.
		const uint dimensions = _theGridCell.dimension();
		uint max_num_subdiv_dim = 0, num_subdiv = 0, max_subdiv_dim = 0;
		
		for(uint i = 0; i < dimensions; i++){
			//Get the number of required subdivisions in this dimension
			//IVAN S ZAPREEV:
			//NOTE: We compute sub_up because we do not want to have insufficient number of subdivisions
			num_subdiv = compute_number_subdiv( sub_up( theRootCellBox[i].upper(), theRootCellBox[i].lower() ) , theMaxCellWidth );

			//Compute the max number of subdivisions and the dimension where to do them
			if( num_subdiv >= max_num_subdiv_dim ){
				max_num_subdiv_dim = num_subdiv;
				max_subdiv_dim = i;
			}
		}
		
		//3. Let the maximum number of subdivisions M has to be done in dimension K  with the total number of
		//   dimensions N: 1 <= K <= N. This means that from this cell down we have to do M splits for dimension K.
		uint needed_num_tree_subdiv = 0;
		//If we need to subdivide in one of the dimensions then
		if( max_num_subdiv_dim != 0 ){
			//3.1 Compute the dimension C for which we had the last split, we should start with the primary cell which is the root of
			//the GridTreeSet because from this cell we begin subdividing in dimension one by one: 1,2,...,N, then again 1,2,...,N.
			//The path to the root of the sub-paving is given by the binary word, its length gives the number of tree subdivisions:
			const uint pathLength = _theGridCell.word().size();
			//If pathLength == 0 then there were no subdivisions in the tree, so we assign last_subdiv_dim == -1
			const int last_subdiv_dim = ( pathLength == 0 ) ? -1 : ( pathLength - 1 )  % dimensions;
			
			//3.2 Compute the needed number of tree subdivisions by first computing how many subdivisions in the tree
			//we need to do to reach and split the dimension K one first time and then we should add the remaining
			//( M - 1 )*N tree subdevisions which will make shure that the K'th dimension is subdivided M times.
			uint first_subdiv_steps;
			if( last_subdiv_dim == static_cast<int>(max_subdiv_dim) ) {
				//If last_subdiv_dim == -1 then we will never get here
				first_subdiv_steps = dimensions; // C == K
			} else {
				//If last_subdiv_dim == -1 then we will add a needed extra subdivision
				first_subdiv_steps = max_subdiv_dim - last_subdiv_dim; // C < K
				if( last_subdiv_dim > static_cast<int>(max_subdiv_dim) ) {
					//If last_subdiv_dim == -1 then we will never get here
					first_subdiv_steps = dimensions - first_subdiv_steps; // C > K
				}
			}
			needed_num_tree_subdiv = first_subdiv_steps + ( max_num_subdiv_dim - 1 ) * dimensions;
		}
		
		//Mince to the computed number of tree levels
		mince(needed_num_tree_subdiv);
	}
	
        GridTreeSubset::operator ListSet<Box>() const {
                ListSet<Box> result(this->cell().dimension());
		
		//IVAN S ZAPREEV:
		//NOTE: Push back the boxes, note that BoxListSet uses a vector, that
		//in its turn uses std:vector to store boxes in it (via push_back method),
		//the latter stores the copies of the boxes, not the references to them.
		for (GridTreeSubset::const_iterator it = this->begin(), end = this->end(); it != end; it++ ) {
			result.push_back((*it).box());
		}
		
		return result;
	}

	GridTreeSubset::operator GridCellListSet() const {
		//GridCellListSet result( this->cell()->grid() );
		
		//IVAN S ZAPREEV:
		//WARNING: Implementation of this conversion operator does not make much sence.
		//The reason is that GridCellListSet contains GridCell elements, which are cells
		//related to the lattice. In our case, enabled cells in the (sub) paving are not
		//on the Lattice, moreover, we can have as small Cells on the grid as we like.
		
		throw NotImplemented(__PRETTY_FUNCTION__);
	}
			
	/*! \brief Serialize the current object state, this is mostly needed for testing */
	string GridTreeSubset::to_string() const{
		//Call the super class's method to get the Cell's information
		string result = _theGridCell.to_string();
		
		//Print the SubPaving's binary tree
		result += "\n The subpaving's tree: " + _pRootTreeNode->tree_to_string();
		
		return result;
	}

/*********************************************GridTreeSet*********************************************/
	
	BinaryTreeNode* GridTreeSet::align_with_cell( const GridCell& theCell, const bool stop_on_enabled, bool & has_stopped ) {
		const uint thisPavingPCellHight = this->cell().height();
		const uint otherPavingPCellHight = theCell.height();
		
		ARIADNE_ASSERT( this->cell().grid()==theCell.grid() );
		
		//The tree node on which we will call add_enabled( theCell.word() ) method to 
		//add the cell. This variable's value might change in the following if clauses
		BinaryTreeNode * pBinaryTreeNode = this->_pRootTreeNode;
		
		if( thisPavingPCellHight > otherPavingPCellHight ){

			//The primary cell of this paving is higher then the one of the other paving.
			//1. We locate the path to the primary cell node common with the other paving
			BinaryWord primaryCellPath = GridCell::primary_cell_path( this->cell().grid().dimension(), thisPavingPCellHight, otherPavingPCellHight );
			
			//2. Locate the binary tree node corresponding to this primnary cell
			uint position = 0;
			while( position < primaryCellPath.size() && !( has_stopped = ( pBinaryTreeNode->is_enabled() && stop_on_enabled ) ) ){
				//Split the node, if it is not a leaf it will not be changed
				pBinaryTreeNode->split();
				//Follow the next path step
				pBinaryTreeNode = (primaryCellPath[position]) ? pBinaryTreeNode->right_node() : pBinaryTreeNode->left_node();
				//Move to the next path element
				position++;
			}
		} else {
			if( thisPavingPCellHight < otherPavingPCellHight ) {
				//The primary cell of this paving is lower then the one in the other paving so this
				//paving's has to be rerooted to another primary cell and then we merge the pavings.
				//1. Compute the path 
				BinaryWord primaryCellPath = GridCell::primary_cell_path( this->cell().grid().dimension(), otherPavingPCellHight, thisPavingPCellHight );
				//2. Substitute the root node of the paiving with the extended tree
				pBinaryTreeNode = ( this->_pRootTreeNode = BinaryTreeNode::prepend_tree( primaryCellPath, this->_pRootTreeNode ) );
				//3. Update the GridCell that corresponds to the root of this GridTreeSubset
				this->_theGridCell = GridCell( this->_theGridCell.grid(), otherPavingPCellHight, BinaryWord() );
			} else {
				//If we are rooted to the same primary cell, then there
				//is nothing to be done, except adding the enabled cell
			}
		}
		return pBinaryTreeNode;
	}
	
	template<class Set> void GridTreeSet::adjoin_outer_approximation( BinaryTreeNode * pBinaryTreeNode, const uint primary_cell_height,
												const uint max_mince_depth,  const Set& theSet, BinaryWord * pPath ){
		//Compute the cell correspomding to the current node
		GridCell theCurrentCell( GridTreeSubset::_theGridCell.grid(), primary_cell_height, *pPath );

		if( bool( theSet.disjoint( theCurrentCell.box() ) ) ){
			//DO NOTHING: We are in the node's representation in the original space is disjoint
			//            from theSet and thus there will be nothing added to this cell.
		} else {
			//This node's cell is not disjoint from theSer, thus it is possible to adjoin elements
			//of its outer approximation, unless this node is already and enabled leaf node.
			if( pBinaryTreeNode->is_enabled() ){ //NOTE: A non-leaf node can not be enabled so this check suffices
				//DO NOTHING: If it is enabled, then we can not add anything new to it.
			} else {
				//If the node is no enabled, so may be we can add something from the outer approximation of theSet.
				if( pPath->size() < max_mince_depth ){
					//Since we still do not have the finest cells for the outer approximation of theSet, we split 
					pBinaryTreeNode->split(); //NOTE: splitting a non-leaf node does not do any harm
					//Check the left branch
					pPath->push_back(false);
					adjoin_outer_approximation( pBinaryTreeNode->left_node(), primary_cell_height, max_mince_depth, theSet, pPath );
					//Check the right branch
					pPath->push_back(true);
					adjoin_outer_approximation( pBinaryTreeNode->right_node(), primary_cell_height, max_mince_depth, theSet, pPath );
				} else {
					//We should not mince any further, so since the node is a leaf and
					//it's cell is not disjoint from theSet, we mark the node as enabled.
					if( pBinaryTreeNode->is_leaf() ){
						//If the node is not leaf, then we make it an enabled one
						pBinaryTreeNode->make_leaf(true);
					} else {
						//Just make the node enabled
						pBinaryTreeNode->set_enabled();
					}
				}
			}
		}
		//Return to the previous level, since the initial call is made
		//with the empty word, we check that it is not yet empty.
		if( pPath->size() > 0 ) {
			pPath->pop_back();
		}
	}

	void GridTreeSet::adjoin_outer_approximation( const CompactSetInterface& theSet, const uint depth ) {
		Grid theGrid( this->cell().grid() );
		ARIADNE_ASSERT( theSet.dimension() == this->cell().dimension());
		
                throw NotImplemented(__PRETTY_FUNCTION__);
                // PIETER COLLINS: TODO: GridBlock not in current version; preferably reimplement without
		//1. Computes the outer approximation (on the grid) of theSet
		//const GridBlock theSetGridBlock = outer_approximation( theSet.bounding_box(), theGrid );
		const GridCell theSetGridBlock = outer_approximation( theSet.bounding_box(), theGrid, depth );
		//2. Allocates the paving for this outer approximation
		//IVAN S. ZAPREEV
		//WARNING: We need to have the bounding box on the grid, in order to compute proper primary cell for the GridTreeSet
		// A) We need to get the box representing the block theSetGridBlock in the grid theGrid
		// B) We create a trivial grid with with no scaling
		// C) We use this grid and the Lattice block of the theSetGridBlock to create a new GridBlock
		// D) The resulting GridBlock will be in the grid theGrid i.e. it's bounding_box() method will
		//     return us it's box in the grid theGrid but not in the original space
		
                // FIXME: No LatticeSet defined
                //GridTreeSet theGridPavingOfSetOuterBox( theGrid, GridBlock( Grid( theGrid.dimension(), R(1.0) ), theSetGridBlock.lattice_set() ).bounding_box() );
		GridTreeSet theGridPavingOfSetOuterBox(theGrid);

		GridCell theSetBoxCell( theGridPavingOfSetOuterBox.cell() );

		//3. Align this paving and paving enclosing the provided set
		bool has_stopped = false;
		BinaryTreeNode* pBinaryTreeNode = align_with_cell( theSetBoxCell, true, has_stopped );
		
		//If the outer aproximation of the bounding box of the provided set is enclosed
		//in an enabled cell of this paving, then there is nothing to be done. The latter
		//is because adjoining the outer approx of the set will not change this paving.
		if( ! has_stopped ){
			//Compute the depth to which we must mince the outer approximation of the adjoining set.
			//This depth is relative to the root of the constructed paving, which has been alligned
			//with the binary tree node pBinaryTreeNode.
			const uint max_mince_depth = theSetBoxCell.height() * theSetBoxCell.dimension() + depth;
			
			//Adjoin the outer approximation, comoputing it on the fly.
			BinaryWord * pEmptyPath = new BinaryWord(); 
			adjoin_outer_approximation( pBinaryTreeNode, theSetBoxCell.height(), max_mince_depth, theSet, pEmptyPath );
			delete pEmptyPath;
		}
	}

	void GridTreeSet::restrict( const GridTreeSubset& theOtherSubPaving ) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}
	
	void GridTreeSet::remove( const GridCell& theCell ) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}
	
	void GridTreeSet::remove( const GridTreeSubset& theOtherSubPaving ) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}

/*************************************FRIENDS OF BinaryTreeNode*************************************/

/*************************************FRIENDS OF GridTreeSet*****************************************/

	bool subset( const GridCell& theCell, const GridTreeSubset& theSet ) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}
	
	bool overlap( const GridCell& theCell, const GridTreeSubset& theSet ) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}

	bool subset( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 ) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}
	
	bool overlap( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 ) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}

	bool subset( const Box& theBox, const GridTreeSubset& theSet ) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}
	
	bool disjoint( const Box& theBox, const GridTreeSubset& theSet ) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}
	
	bool intersects( const Box& theBox, const GridTreeSubset& theSet ) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}

	GridTreeSet join( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 ) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}
	
	GridTreeSet intersection( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 ) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}
	
	GridTreeSet difference( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 ) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}

        std::ostream& operator<<(std::ostream& os, const GridTreeSet& theSet) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}
 


} // namespace Ariadne

