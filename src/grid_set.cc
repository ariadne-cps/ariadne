/***************************************************************************
 *            grid_set.cc
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

#include <iostream>
#include <iomanip>

#include "macros.h"
#include "exceptions.h"
#include "stlio.h"
#include "function_set.h"
#include "list_set.h"
#include "grid_set.h"

#include "set_interface.h"

namespace Ariadne {

typedef size_t size_type;


/****************************************Grid**********************************************/
 
struct Grid::Data 
{ 
  Vector<Float> _origin; 
  Vector<Float> _lengths; 
};

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
 
Grid::Grid(uint d)
  : _data(new Data())
{
  Vector<Float> origin(d,Float(0));
  Vector<Float> lengths(d,Float(1));
  this->_create(origin,lengths);
}
 
Grid::Grid(uint d, Float l)
  : _data(new Data())
{
  Vector<Float> origin(d,Float(0));
  Vector<Float> lengths(d,l);
  this->_create(origin,lengths);
}

Grid::Grid(const Vector<Float>& lengths)
  : _data(new Data())
{
  Vector<Float> origin(lengths.size(),0);
  this->_create(origin,lengths);
}
 
Grid::Grid(const Vector<Float>& origin, const Vector<Float>& lengths)
  : _data(new Data())
{
  if(origin.size() != lengths.size()) {
    throw IncompatibleSizes(__PRETTY_FUNCTION__);
  }
  this->_create(origin,lengths);
}

void Grid::_create(const Vector<Float>& origin, const Vector<Float>& lengths) 
{
  this->_data->_origin=origin;
  this->_data->_lengths=lengths;
}

uint Grid::dimension() const
{
  return this->_data->_lengths.size();
}

const Vector<Float>& Grid::origin() const
{
  return this->_data->_origin;
}

const Vector<Float>& Grid::lengths() const
{
  return this->_data->_lengths;
}

Float Grid::coordinate(uint d, dyadic_type x) const 
{
  return add_approx(this->_data->_origin[d],mul_approx(this->_data->_lengths[d],x));
}

Float Grid::subdivision_coordinate(uint d, dyadic_type x) const 
{
  return add_approx(this->_data->_origin[d],mul_approx(this->_data->_lengths[d],x));
}

Float Grid::subdivision_coordinate(uint d, integer_type n) const 
{
  return add_approx(this->_data->_origin[d],mul_approx(this->_data->_lengths[d],n));
}

int Grid::subdivision_index(uint d, const real_type& x) const 
{
  Float half=0.5;
  int n=int(floor(add_approx(div_approx(sub_approx(x,this->_data->_origin[d]),this->_data->_lengths[d]),half)));
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
 
int Grid::subdivision_lower_index(uint d, const real_type& x) const 
{
  int n=int(floor(div_down(sub_down(x,this->_data->_origin[d]),this->_data->_lengths[d])));
  if(x>=add_approx(this->_data->_origin[d],mul_approx(this->_data->_lengths[d],(n+1)))) {
    return n+1;
  } else {
    return n;
  }
}
 
int Grid::subdivision_upper_index(uint d, const real_type& x) const 
{
	int n=int(ceil(div_up(sub_up(x,this->_data->_origin[d]),this->_data->_lengths[d])));
	if(x<=add_approx(this->_data->_origin[d],mul_approx(this->_data->_lengths[d],(n-1)))) {
		return n-1;
	} else {
		return n;
	}
}
 
bool Grid::operator==(const Grid& g) const
{
  if(this->_data==g._data) { 
    return true; 
  } else {
    return this->_data->_origin==g._data->_origin && this->_data->_lengths==g._data->_lengths;
  }
}
 
bool Grid::operator!=(const Grid& g) const
{
  return !(*this==g);
}

array<double> Grid::index(const Vector<Float>& pt) const
{
        array<double> res(pt.size());
        for(size_t i=0; i!=res.size(); ++i) {
                res[i]=subdivision_index(i,pt[i]);
        }
        return res;
}

array<double> Grid::lower_index(const Vector<Interval>& bx) const {
        array<double> res(bx.size());
        for(size_t i=0; i!=res.size(); ++i) {
                res[i]=subdivision_lower_index(i,bx[i].lower());
        }
        return res;
}

array<double> Grid::upper_index(const Vector<Interval>& bx) const {
        array<double> res(bx.size());
        for(size_type i=0; i!=res.size(); ++i) {
                res[i]=subdivision_upper_index(i,bx[i].upper());
        }
        return res;
}

Vector<Float> Grid::point(const array<int>& a) const
{
  Vector<Float> res(a.size());
  for(size_type i=0; i!=res.size(); ++i) {
    res[i]=this->_data->_origin[i]+this->_data->_lengths[i]*a[i];
  }
  return res;
}

Vector<Float> Grid::point(const array<double>& a) const
{
  Vector<float> res(a.size());
  for(size_type i=0; i!=res.size(); ++i) {
          res[i]=this->_data->_origin[i]+this->_data->_lengths[i]*a[i];
  }
  return res;
}

Vector<Interval> Grid::box(const array<double>& lower, const array<double>& upper) const
{
        Vector<Interval> res(lower.size());
        for(size_type i=0; i!=res.size(); ++i) {
                res[i]=Interval(this->subdivision_coordinate(i,lower[i]),
                                this->subdivision_coordinate(i,upper[i]));
        }
        return res;
}

std::ostream& operator<<(std::ostream& os, const Grid& gr) 
{
        os << "Grid( ";
        os << "origin=" << gr.origin() << ", ";
        os << "lengths=" << gr.lengths() << " )";
        return os;
}

/****************************************BinaryTreeNode**********************************************/
	
	bool BinaryTreeNode::has_enabled() const {
		if( is_leaf() ) {
			return is_enabled();
		} else {
			return ( left_node()->has_enabled() || right_node()->has_enabled() );
		}
	}
	
	bool BinaryTreeNode::all_enabled() const {
		if( is_leaf() ) {
			return is_enabled();
		} else {
			return ( left_node()->all_enabled() && right_node()->all_enabled() );
		}
	}

	bool BinaryTreeNode::is_enabled( const BinaryWord & path, const uint position) const {
		bool result = false ;
		
		if( this->is_leaf() ){
			//If we are in an enabled leaf node then the answer is true.
			//Since path.size() >= 0, the node defined by the path is a
			//subnode of an enabled node, thus we return true.
			result = is_enabled();
		} else {
			if( position < path.size() ) {
				//The path is not complete yet and we are in a
				//non-leaf node, so we follow the path in the tree
				if( path[ position ] ) {
					result = right_node()->is_enabled( path, position + 1 );
				} else {
					result = left_node()->is_enabled( path, position + 1 );
				}
			} else {
				//We are somewhere in the tree in a non-leaf node,
				//this node corresponds to the node given by the path
				//If both left and right sub trees are fully enabled,
				//then the cell defined by the binary path is "enabled"
				//in this tree otherwise it is not.
				result = all_enabled();
			}
		}
		return result;
	}

	bool BinaryTreeNode::is_equal_nodes( const BinaryTreeNode * pFirstNode, const BinaryTreeNode * pSecondNode ) {
		bool result = true;
		
		if( pFirstNode != pSecondNode){
			//If the pointers do not reference the same objects
			if( ( pFirstNode != NULL ) && ( pSecondNode != NULL ) ){
				//And both nodes are not null
				if( ! ( (* pFirstNode) == ( * pSecondNode ) ) ){
					//The objects referenced by the pointers do not contain equal data
					result = false;
				}
			} else {
				//One of the nodes is null and the other is not
				result = false;
			}
		}
		return result;
	}
	
	bool BinaryTreeNode::operator==(const BinaryTreeNode & otherNode ) const {
		return ( ( this->_isEnabled == otherNode._isEnabled ) ||
			 ( indeterminate( this->_isEnabled     ) &&
			   indeterminate( otherNode._isEnabled ) ) )		    &&
			is_equal_nodes( this->_pLeftNode , otherNode._pLeftNode )   &&
			is_equal_nodes( this->_pRightNode , otherNode._pRightNode );
	}

	void BinaryTreeNode::restrict( BinaryTreeNode * pThisNode, const BinaryTreeNode * pOtherNode ){
		if( ( pThisNode != NULL ) && ( pOtherNode != NULL ) ){
			if( pThisNode->is_leaf() && pOtherNode->is_leaf() ){
				//Both nodes are leaf nodes: Make a regular AND
				pThisNode->_isEnabled = (pThisNode->_isEnabled && pOtherNode->_isEnabled);
			} else {
				if( !pThisNode->is_leaf() && pOtherNode->is_leaf() ){
					if( pOtherNode->is_enabled() ){
						//DO NOTHING: The restriction will not affect pThisNode
					} else {
						//Turn the node a disabled leaf, since we do AND with false
						pThisNode->make_leaf(false);
					}
				} else {
					if( pThisNode->is_leaf() && !pOtherNode->is_leaf() ){
						if( pThisNode->is_enabled() ){
							//If this node is enabled then copy in the other node
							//Since it will be their intersection any ways
							pThisNode->copy_from( pOtherNode );
						} else {
							//DO NOTHING: The restriction is empty in this case
							//because this node is a disabled leaf
						}
					} else {
						//Both nodes are non-leaf nodes: Go recursively left and right
						restrict( pThisNode->_pLeftNode, pOtherNode->_pLeftNode );
						restrict( pThisNode->_pRightNode, pOtherNode->_pRightNode );
					}
				}
			}
		}
	}
	
	void BinaryTreeNode::remove( BinaryTreeNode * pThisNode, const BinaryTreeNode * pOtherNode ){
		if( ( pThisNode != NULL ) && ( pOtherNode != NULL ) ){
			if( pThisNode->is_leaf() && pOtherNode->is_leaf() ){
				if( pThisNode->is_enabled() && pOtherNode->is_enabled() ){
					//Both nodes are enabled leaf nodes: Make a regular subtraction, i.e. set the false
					pThisNode->_isEnabled = false;
				} else {
					//DO NOTHING: In all other cases there is nothing to be done
				}
			} else {
				if( !pThisNode->is_leaf() && pOtherNode->is_leaf() ){
					if( pOtherNode->is_enabled() ){
						//Turn the node into a disabled leaf, since we subtract all below
						pThisNode->make_leaf(false);
					} else {
						//DO NOTHING: We are trying to remove a disabled node
					}
				} else {
					if( pThisNode->is_leaf() ) {
						if( pThisNode->is_enabled() ){
							//This is an enabled leaf node and so we might subtract smth from it
							//The pOtherNode is not a leaf, due to previous checks so we split
							//pThisNode and then do the recursion as in case on two non-leaf nodes
							pThisNode->split();
						} else {
							//DO NOTHING: We are trying to remove from a disabled leaf node. Whatever
							//we are trying to remove, will not have any effect on the result. We return
							//because in all remaining cases we need to do recursion for sub trees.
							return;
						}
					} else {
						//We will have to do the recursion to remove the leaf nodes
					}
					//Both nodes are non-leaf nodes now: Go recursively left and right
					remove( pThisNode->_pLeftNode, pOtherNode->_pLeftNode );
					remove( pThisNode->_pRightNode, pOtherNode->_pRightNode );
				}
			}
		}
	}
	
	void BinaryTreeNode::restore_node( BinaryTreeNode * pCurrentNode, uint & arr_index, uint & leaf_counter,
					const BooleanArray& theTree, const BooleanArray& theEnabledCells) {
		//If we are not done with the the tree yet
		if( arr_index < theTree.size() ) {
			//If we are in a non-leaf node then go further
			if( theTree[arr_index] ) {
					pCurrentNode->split();
					//IVAN S. ZAPREEV:
					//NOTE: We assume a correct input, i.e. both children are present
					//NOTE: We increase the arr_index before calling recursion for each
					//      of the subnodes of the tree
					restore_node( pCurrentNode->_pLeftNode, ++arr_index, leaf_counter, theTree, theEnabledCells );
					restore_node( pCurrentNode->_pRightNode, ++arr_index, leaf_counter, theTree, theEnabledCells );
			} else {
				//If we are in a leaf node then chek if it needs to be enabled/disabled
				pCurrentNode->_isEnabled = theEnabledCells[leaf_counter];
				
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
				throw std::runtime_error( __PRETTY_FUNCTION__ );
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
			throw std::runtime_error( __PRETTY_FUNCTION__ );
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

        size_t BinaryTreeNode::count_enabled_leaf_nodes( const BinaryTreeNode* pNode ) {
                if(pNode->is_leaf()) { 
                        return pNode->is_enabled() ? 1u : 0u; 
                } else {
                        return count_enabled_leaf_nodes(pNode->left_node())
                            + count_enabled_leaf_nodes(pNode->right_node());
                }
        }
	
	void BinaryTreeNode::tree_to_binary_words( BinaryWord & tree, BinaryWord & leaves ) const {
		if( is_leaf() ) {
			tree.push_back( false );
			leaves.push_back( definitely( _isEnabled ) );
		} else {
			tree.push_back( true );
			_pLeftNode->tree_to_binary_words( tree, leaves );
			_pRightNode->tree_to_binary_words( tree, leaves );
		}
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

/****************************************GridTreeConstIterator************************************/

        GridTreeConstIterator::GridTreeConstIterator(  ) :
                _pGridTreeCursor()  {  
            
        }

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

        bool GridCell::operator==( const GridCell& other ) const {
                return this->_theGrid == other._theGrid && 
                        this->_theHeight == other._theHeight &&
			this->_theWord == other._theWord;
        }      

        // PIETER: TODO: Ivan, please take a look at operator< to make sure it's
        // independant of the presentation of the cell.
        bool GridCell::operator<(const GridCell& other) const {
                ARIADNE_ASSERT( this->_theGrid == other._theGrid );
                return this->_theHeight < other._theHeight || 
                        ( this->_theHeight == other._theHeight && this->_theWord < other._theWord );
        }      

	Box GridCell::primary_cell( const uint theHeight, const dimension_type dimensions ) {
		int leftBottomCorner = 0, rightTopCorner = 1;
		//The zero level coordinates are known, so we need to iterate only for higher level primary cells
		for(uint i = 1; i <= theHeight; i++){
			primary_cell_at_height(i, leftBottomCorner, rightTopCorner);
		}
		// 1.2 Constructing and return the box defining the primary cell (relative to the grid).
		
		return Box( Vector< Interval >( dimensions, Interval( leftBottomCorner, rightTopCorner ) ) );
	}

	uint GridCell::smallest_primary_cell_height( const Box& theBox ) {
                const dimension_type dimensions = theBox.dimension();
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
		}while( true );
		
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

        bool GridTreeSubset::superset( const Box& theBox ) const {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}
	
        bool GridTreeSubset::subset( const Box& theBox ) const {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}
	
        bool GridTreeSubset::disjoint( const Box& theBox ) const {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}
	
        bool GridTreeSubset::intersects( const Box& theBox ) const {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}
	
/*********************************************GridTreeSet*********************************************/

	void GridTreeSet::up_to_primary_cell( const uint toPCellHeight ){
		const uint fromPCellHeight = this->cell().height();
		
		//The primary cell of this paving is lower then the one in the other paving so this
		//paving's has to be rerooted to another primary cell and then we merge the pavings.
		//1. Compute the path 
		BinaryWord primaryCellPath = GridCell::primary_cell_path( this->cell().grid().dimension(), toPCellHeight, fromPCellHeight );
		//2. Substitute the root node of the paiving with the extended tree
		this->_pRootTreeNode = BinaryTreeNode::prepend_tree( primaryCellPath, this->_pRootTreeNode );
		//3. Update the GridCell that corresponds to the root of this GridTreeSubset
		this->_theGridCell = GridCell( this->_theGridCell.grid(), toPCellHeight, BinaryWord() );
	}

	BinaryTreeNode* GridTreeSet::align_with_cell( const uint otherPavingPCellHeight, const bool stop_on_enabled, const bool stop_on_disabled, bool & has_stopped ) {
		const uint thisPavingPCellHeight = this->cell().height();
		
		//The current root node of the GridTreeSet
		BinaryTreeNode * pBinaryTreeNode = this->_pRootTreeNode;
		
		if( thisPavingPCellHeight > otherPavingPCellHeight ){

			//The primary cell of this paving is higher then the one of the other paving.
			//1. We locate the path to the primary cell node common with the other paving
			BinaryWord primaryCellPath = GridCell::primary_cell_path( this->cell().grid().dimension(), thisPavingPCellHeight, otherPavingPCellHeight );
			
			//2. Locate the binary tree node corresponding to this primnary cell
			uint position = 0;
			while( position < primaryCellPath.size() && !( has_stopped = (
					( pBinaryTreeNode->is_enabled() && stop_on_enabled ) ||
					( pBinaryTreeNode->is_disabled() && stop_on_disabled ) ) ) ){
				//Split the node, if it is not a leaf it will not be changed
				pBinaryTreeNode->split();
				//Follow the next path step
				pBinaryTreeNode = (primaryCellPath[position]) ? pBinaryTreeNode->right_node() : pBinaryTreeNode->left_node();
				//Move to the next path element
				position++;
			}
		} else {
			if( thisPavingPCellHeight < otherPavingPCellHeight ) {
				up_to_primary_cell( otherPavingPCellHeight );
				pBinaryTreeNode = this->_pRootTreeNode;
			} else {
				//If we are rooted to the same primary cell, then there
				//is nothing to be done, except adding the enabled cell
			}
		}
		return pBinaryTreeNode;
	}
	
        void GridTreeSet::adjoin_outer_approximation( BinaryTreeNode * pBinaryTreeNode, const uint primary_cell_height,
                                                      const uint max_mince_depth,  const CompactSetInterface& theSet, BinaryWord * pPath ){
                //Compute the cell correspomding to the current node
                GridCell theCurrentCell( GridTreeSubset::_theGridCell.grid(), primary_cell_height, *pPath );

                const OpenSetInterface* pOpenSet=dynamic_cast<const OpenSetInterface*>(static_cast<const SetInterfaceBase*>(&theSet));

                if( bool( theSet.disjoint( theCurrentCell.box() ) ) ) {
                        //DO NOTHING: We are in the node's representation in the original space is disjoint
                        //            from pSet and thus there will be nothing added to this cell.
                } else if( pOpenSet && bool( pOpenSet->superset( theCurrentCell.box() ) ) ) {
                        pBinaryTreeNode->make_leaf(true);
                } else {
                        //This node's cell is not disjoint from theSet, thus it is possible to adjoin elements
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
                                        if( !pBinaryTreeNode->is_leaf() ){
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

        // FIXME: This method can fail if we cannot determine which of a node's children intersects
        // the set. In principle this can be solved by checking if one of the children intersects
	// the set, before doing recursion, and if none intersects then we mark the present node as
	// enabled and stop. Generally speaking, the present algorithm is not wrong, it also gives us
	// a lower approximation, but it is simply less accurate than it could be.
	// TODO:Think of another representation in terms of covers but not pavings, then this problem
	// can be cured in a different fashion.
        void GridTreeSet::adjoin_lower_approximation( BinaryTreeNode * pBinaryTreeNode, const uint primary_cell_height,
							const uint max_mince_depth,  const OvertSetInterface& theSet, BinaryWord * pPath ){
		//Compute the cell correspomding to the current node
		GridCell theCurrentCell( GridTreeSubset::_theGridCell.grid(), primary_cell_height, *pPath );

		if( bool( theSet.intersects( theCurrentCell.box() ) ) ) {
                        if( pPath->size() >= max_mince_depth ) {
				//We should not mince any further.
				//If the cell is not a leaf, then some subset is enabled,
				//so the lower approximation does not add any information.
				//If the cell is a leaf, we mark it as enabled.
				if( pBinaryTreeNode->is_leaf() ) {
					pBinaryTreeNode->set_enabled(); 
				}
			} else {
				//If the node is no enabled, so may be we can add something from the lower approximation of theSet.
				//Since we still do not have the finest cells for the lower approximation of theSet, we split 
				pBinaryTreeNode->split(); //NOTE: splitting a non-leaf node does not do any harm
				//Check the left branch
				pPath->push_back(false);
				adjoin_lower_approximation( pBinaryTreeNode->left_node(), primary_cell_height, max_mince_depth, theSet, pPath );
				//Check the right branch
				pPath->push_back(true);
				adjoin_lower_approximation( pBinaryTreeNode->right_node(), primary_cell_height, max_mince_depth, theSet, pPath );
			}
		}
		//Return to the previous level, since the initial call is made
		//with the empty word, we check that it is not yet empty.
		if( pPath->size() > 0 ) {
			pPath->pop_back();
		}
        }
	
        void GridTreeSet::adjoin_lower_approximation( BinaryTreeNode * pBinaryTreeNode, const uint primary_cell_height,
                                                      const uint max_mince_depth,  const OpenSetInterface& theSet, BinaryWord * pPath ){
		//Compute the cell correspomding to the current node
		GridCell theCurrentCell( GridTreeSubset::_theGridCell.grid(), primary_cell_height, *pPath );

		if( bool( theSet.superset( theCurrentCell.box() ) ) ) {
                        this->mince( max_mince_depth - pPath->size() );
                } else if ( bool( theSet.intersects( theCurrentCell.box() ) ) ) {
                        if( pPath->size() >= max_mince_depth ) {
                                 //We should not mince any further.
                                 //If the cell is not a leaf, then some subset is enabled,
                                 //so the lower approximation does not add any information.
                                 //If the cell is a leaf, we mark it as enabled.
                                if( pBinaryTreeNode->is_leaf() ) {
                                        pBinaryTreeNode->set_enabled(); 
                                }
                        } else {
                         //Since we still do not have the finest cells for the outer approximation of theSet, we split 
                                if( pBinaryTreeNode->is_leaf() ) { 
                                        pBinaryTreeNode->split();
                                }
                                //Check the left branch
                                pPath->push_back(false);
                                adjoin_lower_approximation( pBinaryTreeNode->left_node(), primary_cell_height, max_mince_depth, theSet, pPath );
                                //Check the right branch
                                pPath->push_back(false);
                                adjoin_lower_approximation( pBinaryTreeNode->right_node(), primary_cell_height, max_mince_depth, theSet, pPath );
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
		ARIADNE_ASSERT( theSet.dimension() == this->cell().dimension() );
		
		//1. Compute the smallest GridCell (corresponding to the primary cell)
		//   that encloses the theSet (after it is mapped onto theGrid).
		const GridCell theOverApproxCell = over_approximation( theSet.bounding_box(), theGrid );
		//Compute the height of the primary cell for the outer approximation stepping up by the number of dimensions
		const uint outer_approx_primary_cell_height = theOverApproxCell.height() + theGrid.dimension();

		//2. Align this paving and paving enclosing the provided set
		bool has_stopped = false;
		BinaryTreeNode* pBinaryTreeNode = align_with_cell( outer_approx_primary_cell_height, true, false, has_stopped );
		
		//If the outer aproximation of the bounding box of the provided set is enclosed
		//in an enabled cell of this paving, then there is nothing to be done. The latter
		//is because adjoining the outer approx of the set will not change this paving.
		if( ! has_stopped ){
			//Compute the depth to which we must mince the outer approximation of the adjoining set.
			//This depth is relative to the root of the constructed paving, which has been alligned
			//with the binary tree node pBinaryTreeNode.
			const uint max_mince_depth = outer_approx_primary_cell_height * theGrid.dimension() + depth;
			
			//Adjoin the outer approximation, computing it on the fly.
			BinaryWord * pEmptyPath = new BinaryWord(); 
			adjoin_outer_approximation( pBinaryTreeNode, outer_approx_primary_cell_height, max_mince_depth, theSet, pEmptyPath );

			delete pEmptyPath;
		}
	}
	
        // FIXME: The following three methods should adjoin an approximation at the given height and depth.
        // The use of GridCell(...) is a hack. In essence it would be better to align the binary tree with
	// the heightdepth and then do things like in:
	//	adjoin_lower_approximation( const OvertSetInterface& , const Box& , const uint  )
	// TODO:Think of another representation in terms of covers but not pavings, then the implementation
	// will be different, this is why, for now we do not fix these things.
	void GridTreeSet::adjoin_lower_approximation( const LocatedSetInterface& theSet, const uint heightdepth ) {
                this->adjoin_lower_approximation( theSet, GridCell( this->grid(), heightdepth, BinaryWord() ).box(), heightdepth );
        }
	
	void GridTreeSet::adjoin_lower_approximation( const OvertSetInterface& theSet, const uint heightdepth ) {
                this->adjoin_lower_approximation( theSet, GridCell( this->grid(), heightdepth, BinaryWord()).box(),heightdepth );
        }
	
	void GridTreeSet::adjoin_inner_approximation( const OpenSetInterface& theSet, const uint heightdepth ) {
                throw NotImplemented(__PRETTY_FUNCTION__);
        }
	
        void GridTreeSet::adjoin_lower_approximation( const OvertSetInterface& theSet, const Box& theBoundingBox, const uint depth ) {

		Grid theGrid( this->cell().grid() );
		ARIADNE_ASSERT( theSet.dimension() == this->cell().dimension() );
                  
                //1. Compute the smallest GridCell (corresponding to the primary cell)
		//   that encloses the theSet (after it is mapped onto theGrid).
		const GridCell theOverApproxCell = over_approximation( theBoundingBox, theGrid );
		//Compute the height of the primary cell for the outer approximation stepping up by the number of dimensions
		const uint outer_approx_primary_cell_height = theOverApproxCell.height() + theGrid.dimension();

		//2. Align this paving and paving enclosing the provided set
		bool has_stopped = false;
		BinaryTreeNode* pBinaryTreeNode = align_with_cell( outer_approx_primary_cell_height, true, false, has_stopped );
		
		//If the lower aproximation of the bounding box of the provided set is enclosed
		//in an enabled cell of this paving, then there is nothing to be done. The latter
		//is because adjoining the outer approx of the set will not change this paving.
		if( ! has_stopped ){
			//Compute the depth to which we must mince the outer approximation of the adjoining set.
			//This depth is relative to the root of the constructed paving, which has been alligned
			//with the binary tree node pBinaryTreeNode.
			const uint max_mince_depth = outer_approx_primary_cell_height * theGrid.dimension() + depth;
			
			//Adjoin the outer approximation, computing it on the fly.
			BinaryWord * pEmptyPath = new BinaryWord(); 
                        //const RegularSetInterface* theRegularVersionOfSet = dynamic_cast<const RegularSetInterface*>(&theSet);
                        const OpenSetInterface* theOpenVersionOfSet = dynamic_cast<const OpenSetInterface*>(&theSet);
                        //const LocatedSetInterface* theLocatedVersionOfSet = dynamic_cast<const LocatedSetInterface*>(&theSet);
                        const OvertSetInterface* theOvertVersionOfSet = dynamic_cast<const OvertSetInterface*>(&theSet);
                        if( theOpenVersionOfSet ) {
                                adjoin_lower_approximation( pBinaryTreeNode, outer_approx_primary_cell_height, max_mince_depth, *theOpenVersionOfSet, pEmptyPath );
                        } else {
                                adjoin_lower_approximation( pBinaryTreeNode, outer_approx_primary_cell_height, max_mince_depth, *theOvertVersionOfSet, pEmptyPath );
                        }
			delete pEmptyPath;
		}
	}

	void GridTreeSet::restrict_to_lower( const GridTreeSubset& theOtherSubPaving ){
		//The root of the binary tree of the current Paving
		BinaryTreeNode * pBinaryTreeNode = this->_pRootTreeNode;
		
		//The primary cell of this paving is higher then the one of the other paving.
		//1. We locate the path to the primary cell node common with the other paving
		BinaryWord rootNodePath = GridCell::primary_cell_path( this->cell().grid().dimension(), this->cell().height(), theOtherSubPaving.cell().height() );
		
		//2. Add the suffix path from the primary cell to the root node of
		//theOtherSubPaving. This is needed in order to be able to reach this root.
		rootNodePath.append( theOtherSubPaving.cell().word() );
		
		//3. Restrict this binary tree to the other one assuming the path prefix rootNodePath
		uint position = 0;
		//This will point to the children nodes that do not have chance for being in the restriction
		BinaryTreeNode * pBranchToDisable = NULL;
		//Iterate the path, get to the root cell of theOtherSubPaving
		while( position < rootNodePath.size() ){
			if( pBinaryTreeNode->is_leaf() ){
				//If we are in the leaf node then
				if( pBinaryTreeNode->is_disabled() ){
					//if it is disabled, then the intersection with
					//the other set is empty so we just return
					return;
				} else {
					//If it is an enabled leaf node then, because we still need to go further to reach the root cell of
					//theOtherSubPaving, we split this node and disable the leaf that does not intersect with the other set
					pBinaryTreeNode->split();
				}
			} else {
				//If this is not a leaf node then we need to follow the path and disable the child
				//node that is not on the path, because it will not make it into the intersection
			}
			//Follow the path and diable the other branch
			if( rootNodePath[position] ){
				pBranchToDisable = pBinaryTreeNode->left_node();
				pBinaryTreeNode = pBinaryTreeNode->right_node();
			} else {
				pBranchToDisable = pBinaryTreeNode->right_node();
				pBinaryTreeNode = pBinaryTreeNode->left_node();
			}
			//Move the unused branch into a disabled leaf node
			pBranchToDisable->make_leaf(false);
			
			//Move to the next path element
			position++;
		}
		if( pBinaryTreeNode->is_enabled() ){
			//If we ended up in a leaf node that is disabled, this meand that it is the only enabled node
			//in this GridTreeSet. At this point it corresponds to the root node of the theOtherSubPaving
			//and since we need to do the restriction to that set, we just need to copy it to this node.
			pBinaryTreeNode->copy_from( theOtherSubPaving.binary_tree() );
			//Return since there is nothing else we need to do
			return;
		} else {
			if( pBinaryTreeNode->is_disabled() ){
				//if we are in a disabled leaf node the the result of the
				//restriction is an empty set and we just return
				return;
			} else {
				//If it is two binary tree we have, then we need to restrict
				//pBinaryTreeNode to theOtherSubPaving._pRootTreeNode
				BinaryTreeNode::restrict( pBinaryTreeNode, theOtherSubPaving.binary_tree() );
			}
		}
	}

	void GridTreeSet::remove_from_lower( const GridTreeSubset& theOtherSubPaving ){
		//The root of the binary tree of the current Paving
		BinaryTreeNode * pBinaryTreeNode = this->_pRootTreeNode;
		
		//The primary cell of this paving is higher then the one of the other paving.
		//1. We locate the path to the primary cell node common with the other paving
		BinaryWord rootNodePath = GridCell::primary_cell_path( this->cell().grid().dimension(), this->cell().height(), theOtherSubPaving.cell().height() );
		
		//2. Add the suffix path from the primary cell to the root node of
		//theOtherSubPaving. This is needed in order to be able to reach this root.
		rootNodePath.append( theOtherSubPaving.cell().word() );
		
		//3. Remove theOtherSubPaving from this binary tree assuming the path prefix rootNodePath
		uint position = 0;
		//Iterate the path, get to the root cell of theOtherSubPaving
		while( position < rootNodePath.size() ){
			if( pBinaryTreeNode->is_leaf() ){
				//If we are in the leaf node then
				if( pBinaryTreeNode->is_disabled() ){
					//if it is disabled, then we are removing from not-enabled cells
					//so there is nothing to be done and thus we simply terminate
					return;
				} else {
					//If it is an enabled leaf node then, because we still need to go further
					//to reach the root cell of theOtherSubPaving, we split this node
					pBinaryTreeNode->split();
				}
			} else {
				//If this is not a leaf node then we need to follow the path and then do removal
			}
			//Follow the path and do the removal, the other branch stays intact
			pBinaryTreeNode = (rootNodePath[position]) ? pBinaryTreeNode->right_node() : pBinaryTreeNode->left_node();
			
			//Move to the next path element
			position++;
		}
		if( pBinaryTreeNode->is_disabled() ){
			//If we are in a disabled leaf node the the result of the
			//removal does not change this set and so we just return
			return;
		} else {
			//We have two aligned binary trees we have to subtract enabled
			//nodes of theOtherSubPaving._pRootTreeNode from pBinaryTreeNode
			BinaryTreeNode::remove( pBinaryTreeNode, theOtherSubPaving.binary_tree() );
		}
	}
	
	void GridTreeSet::restrict( const GridTreeSubset& theOtherSubPaving ) {
		const uint thisPavingPCellHeight = this->cell().height();
		const uint otherPavingPCellHeight = theOtherSubPaving.cell().height();
		
		ARIADNE_ASSERT( this->grid() == theOtherSubPaving.grid() );

		//In case theOtherSubPaving has the primary cell that is higher then this one
		//we extend it, i.e. reroot it to the same height primary cell.		
		if( thisPavingPCellHeight < otherPavingPCellHeight ){
			up_to_primary_cell( otherPavingPCellHeight );
		}
		
		//Now it is simple to restrict this set to another, since this set's
		//primary cell is not lower then for the other one
		restrict_to_lower( theOtherSubPaving );
	}
	
	void GridTreeSet::remove( const GridCell& theCell ) {
		ARIADNE_ASSERT( this->grid() == theCell.grid() );
		
		//If needed, extend the tree of this paving and then find it's the
		//primary cell common with the primary cell of the provided GridCell
		//If we encounter a disabled node then we do not move on, since then
		//there is nothing to remove, otherwise we split the node and go down
		bool has_stopped = false;
		BinaryTreeNode* pCurrentPrimaryCell = align_with_cell( theCell.height(), false, true, has_stopped );
		
		if( ! has_stopped ) {
			//Follow theCell.word() path in the tree rooted to pCommonPrimaryCell,
			//do that until we encounter a leaf node, then stop
			BinaryWord path = theCell.word();
			int position = 0;
			for(; position < path.size(); position++ ) {
				//If we are in a leaf node then
				if( pCurrentPrimaryCell->is_leaf() ) {
					break;
				}
				pCurrentPrimaryCell = path[position] ? pCurrentPrimaryCell->right_node() : pCurrentPrimaryCell->left_node();
			}
			
			//Check if we stopped because it was a leaf node
			if( pCurrentPrimaryCell->is_leaf() ) {
				if( pCurrentPrimaryCell->is_enabled() ) {
					//If the node is a leaf and enabled then continue following the path by splitting nodes
					for(; position < path.size(); position++ ) {
						pCurrentPrimaryCell->split();
						pCurrentPrimaryCell = path[position] ? pCurrentPrimaryCell->right_node() : pCurrentPrimaryCell->left_node();
					}
					//Disable the sub-tree rooted to the tree node we navigated to
					pCurrentPrimaryCell->set_disabled();
				} else {
					//DO NOTHING: The leaf node turns out to be already off,
					//so there is nothing to remove from
				}
			} else {
				//We followed the path to the cell, but we are still in some non-leaf node at this
				//point it does no matter what is below, and we just remove the entire sub-tree.
				pCurrentPrimaryCell->make_leaf( false );
			}
		} else {
			//DO NOTHING: If we stopped it means that we were in a
			//disabled node, so there is nothing to remove from
		}
	}
	
	void GridTreeSet::remove( const GridTreeSubset& theOtherSubPaving ) {
		const uint thisPavingPCellHeight = this->cell().height();
		const uint otherPavingPCellHeight = theOtherSubPaving.cell().height();
		
		ARIADNE_ASSERT( this->grid() == theOtherSubPaving.grid() );

		//In case theOtherSubPaving has the primary cell that is higher then this one
		//we extend it, i.e. reroot it to the same height primary cell.		
		if( thisPavingPCellHeight < otherPavingPCellHeight ){
			up_to_primary_cell( otherPavingPCellHeight );
		}
		
		//Now it is simple to remove theOtherSubPaving elements from this set,
		//since this set's primary cell is not lower then for the other one.
		remove_from_lower( theOtherSubPaving );
	}
	
	void GridTreeSet::restrict_to_height( const uint theHeight ) {
		const uint thisPavingPCellHeight = this->cell().height();
		
		if( thisPavingPCellHeight > theHeight){
			std::cerr << "Warning: restricting GridTreeSet of height " << this->cell().height() << " to height " << theHeight << ".\n";

			BinaryWord pathToPCell = GridCell::primary_cell_path( this->dimension(), thisPavingPCellHeight, theHeight );
			
			//Go throught the tree and disable all the leaves that
			//are not rooted to the primary cell defined by this path
			BinaryTreeNode * pCurrentNode = _pRootTreeNode;
			for( int i = 0; i < pathToPCell.size(); i++ ) {
				if( pCurrentNode->is_leaf() ){
					if( pCurrentNode->is_enabled() ){
						//If we are in an enabled leaf node then we split.
						//There are still cells to remove.
						pCurrentNode->split();
					} else {
						//If we are in a disabled leaf node then we stop.
						//There is nothing to be done, because there are
						//no more enabled cells to remove
						break;
					}
				}
				//If we are here, then we are in a non-leaf node
				if( pathToPCell[i] ){
					//Go to the right, and remove anything on the left
					pCurrentNode->left_node()->make_leaf(false);
					pCurrentNode = pCurrentNode->right_node();
				} else {
					//Go to the left, and remove anything on the right
					pCurrentNode->right_node()->make_leaf(false);
					pCurrentNode = pCurrentNode->left_node();
				}
			}
		}
	}
	
/*************************************FRIENDS OF BinaryTreeNode*************************************/

/*************************************FRIENDS OF GridCell*****************************************/
	
        GridCell over_approximation(const Box& theBox, const Grid& theGrid) {
		Box theDyadicBox( theBox.size() );
		//Convert the box to theGrid coordinates
		for( uint i = 0; i != theBox.size(); ++i ) {
			theDyadicBox[i] = ( theBox[i] - theGrid.origin()[i] ) / theGrid.lengths()[i];
		}
		//Compute the smallest primary cell, enclosing this grid
		uint height = GridCell::smallest_primary_cell_height( theDyadicBox );
		//Create the GridCell, corresponding to this cell
		return GridCell( theGrid, height, BinaryWord() );
        }
	
	bool subset(const GridCell& theCellOne, const GridCell& theCellTwo, BinaryWord * pPathPrefixOne, BinaryWord * pPathPrefixTwo, uint *pPrimaryCellHeight ) {
		//Test that the Grids are equal
		ARIADNE_ASSERT( theCellOne.grid() == theCellTwo.grid() );

		//Test that the binary words are empty, otherwise the results of computations are undefined
		bool isDefault = false;
		if( ( pPathPrefixOne == NULL ) && ( pPathPrefixTwo == NULL ) ) {
			//If both parameters are null, then the user does not
			//need this data, so we allocate tepmorary objects
			pPathPrefixOne = new BinaryWord();
			pPathPrefixTwo = new BinaryWord();
			isDefault = true;
		} else {
			if( ( pPathPrefixOne == NULL ) || ( pPathPrefixTwo == NULL ) ) {
				//TODO: Find some better way to notify the user that
				//only one of the required pointers is not NULL.
				ARIADNE_ASSERT( false );
			} else {
				//DO NOTHING: Both parameters are not NULL, so it
				//is user data and isDefault should stay false.
			}
		}
		
		//Check if theCellOne is a subset of theCellTwo:
		
		//01 Align the cell's primary cells by finding path prefixes
		uint primary_cell_height;
		if( theCellOne.height() < theCellTwo.height() ) {
			primary_cell_height = theCellTwo.height();
			pPathPrefixOne->append( GridCell::primary_cell_path( theCellOne.grid().dimension(), theCellTwo.height(), theCellOne.height() ) );
		} else {
			primary_cell_height = theCellOne.height();
			if( theCellOne.height() > theCellTwo.height() ){
				pPathPrefixTwo->append( GridCell::primary_cell_path( theCellOne.grid().dimension(), theCellOne.height(), theCellTwo.height() ) );
			} else {
				//DO NOTHING: The cells are rooted to the same primary cell
			}
		}
		if( pPrimaryCellHeight != NULL ){
			*pPrimaryCellHeight = primary_cell_height;
		}
		
		//02 Add the rest of the paths to the path prefixes to get the complete paths
		pPathPrefixOne->append( theCellOne.word() );
		pPathPrefixTwo->append( theCellTwo.word() );
		
		//03 theCellOne is a subset of theCellTwo if pathPrefixTwo is a prefix of pathPrefixOne
		bool result = pPathPrefixTwo->is_prefix( *pPathPrefixOne );
		
		//04 If working with default parameters, then release memory
		if( isDefault ) {
			delete pPathPrefixOne;
			delete pPathPrefixTwo;
		}

		return result;
	}

/*************************************FRIENDS OF GridTreeSubset*****************************************/
	
	bool subset( const GridCell& theCell, const GridTreeSubset& theSet ) {
		bool result = false;
		
		//Test that the Grids are equal
		ARIADNE_ASSERT( theCell.grid() == theSet.grid() );
		
		//Test if theCell is a subset of theSet, first check
		//if the cell can be a subset of the given tree.
		BinaryWord pathPrefixCell, pathPrefixSet;
		if( subset( theCell, theSet.cell(), &pathPrefixCell, &pathPrefixSet) ) {
			//It can and thus pathPrefixSet is a prefix of pathPrefixCell. Both
			//paths start in the same primary cell. Also note that pathPrefixSet
			//is a path from primary_cell_height to the root node of the binary
			//tree in theSet. Therefore, removing 0..pathPrefixSet.size() elements
			//from pathPrefixCell will give us a path to theCell in the tree of theSet.
			pathPrefixCell.erase( pathPrefixCell.begin(), pathPrefixCell.begin() + pathPrefixSet.size() );
			
			//Check that the cell given by pathPrefixCell is enabled in the tree of theSet
			result = theSet.binary_tree()->is_enabled( pathPrefixCell );
		} else {
			//DO NOTHING: the cell is a strict superset of the tree
		}
		return result;
	}
	
	bool overlap( const GridCell& theCell, const GridTreeSubset& theSet ) {
		bool result = false;
		
		//Test that the Grids are equal
		ARIADNE_ASSERT( theCell.grid() == theSet.grid() );
		
		//If the primary cell of the theCell is lower that that of theSet
		//Then we re-root theCell to the primary cell theSet.cell().height()
		const GridCell * pWorkGridCell;
		const int theSetsPCellHeigh = theSet.cell().height();
		if( theSetsPCellHeigh > theCell.height() ) {
			//Compute the path from the primary cell of theSet to the primary cell of theCell
			BinaryWord pathFromSetsPCellToCell = GridCell::primary_cell_path( theCell.dimension(), theSetsPCellHeigh, theCell.height() );
			pathFromSetsPCellToCell.append( theCell.word() );
			pWorkGridCell = new GridCell( theCell.grid(), theSetsPCellHeigh, pathFromSetsPCellToCell );
		} else {
			pWorkGridCell = &theCell;
		}

		//Compute the path for the primary cell of theCell to the primary cell of theSet
		BinaryWord pathFromPCellCellToSetsRootNode = GridCell::primary_cell_path( theCell.dimension(), pWorkGridCell->height(), theSetsPCellHeigh );
		//Append the path from the primary cell node to the root binary tree node of theSet
		pathFromPCellCellToSetsRootNode.append( theSet.cell().word() );
		
		const BinaryWord & workCellWord = pWorkGridCell->word();
		if( workCellWord.is_prefix( pathFromPCellCellToSetsRootNode ) ) {
			//If the path (from some primary cell) to the cell is a prefix
			//of the path (from the same primary cell) to the root of the
			//sub-paving, then theCell contains theSet
			
			result = theSet.binary_tree()->has_enabled();
		} else {
			if( pathFromPCellCellToSetsRootNode.is_prefix( workCellWord ) ) {
				//If the path (from some primary cell) to the root of the
				//binary tree node (defining theSet) is the prefix of the
				//path (from the same primary cell) to the cell, then the
				//cell might be somwhere within the tree and we should
				//check if it intersects with the tree
				
				const BinaryTreeNode *pCurrentNode = theSet.binary_tree();
				//Note that, pathFromPCellCellToSetsRootNode.size() < workCellWord.size()
				//Because we already checked for workCellWord.is_prefix( pathFromPCellCellToSetsRootNode )
				//Here we try to find the node corresponding to theCell in the binary tree of theSet
				//in case we encounter a leaf node then we just stop, because it is enough information for us
				for( int i = pathFromPCellCellToSetsRootNode.size(); i < workCellWord.size(); i++ ) {
					if( pCurrentNode->is_leaf() ) {
						//We reached the leaf node and theCell is it's subset, so we stop now
						break;
					} else {
						//Follow the path to the node corresponding to theCell within the binary tree of theSet
						pCurrentNode = ( workCellWord[i] ? pCurrentNode->right_node() : pCurrentNode->left_node() );
					}
				}
				
				//At this point we have the following cases:
				// 1. pCurrentNode - is a leaf node, contains theCell as a subset, is an enabled node
				// RESULT: theSet and theCell overlap
				// 2. pCurrentNode - is a leaf node, contains theCell as a subset, is a disabled node
				// RESULT: theSet and theCell do not overlap
				// 3. pCurrentNode - is a non-leaf node, corresponds to theCell, contains enabled sub-nodes
				// RESULT: theSet and theCell overlap
				// 4. pCurrentNode - is a non-leaf node, corresponds to theCell, contains no enabled sub-nodes
				// RESULT: theSet and theCell do not overlap
				result = pCurrentNode->has_enabled();
			} else {
				//DO NOTHING: The pathes to the cell and to the root of
				//the binary tree (from the same primary cell) diverge
				//This means that theCell and theSet do not overlap
			}
		}
		
		if( theSetsPCellHeigh > theCell.height() ) {
			delete pWorkGridCell;
		}
		
		return result;
	}

	bool subset( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 ) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}
	
	bool overlap( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 ) {
		throw NotImplemented(__PRETTY_FUNCTION__);
	}

/*************************************FRIENDS OF GridTreeSet*****************************************/

        GridTreeSet outer_approximation(const Box& theBox, const Grid& theGrid, const uint depth) {
          return outer_approximation(ImageSet(theBox),theGrid,depth);
        }
	
	GridTreeSet join( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 ) {
		//Test that the Grids are equal
		ARIADNE_ASSERT( theSet1.grid() == theSet2.grid() );
		
		//Compute the highest primary cell 
		const uint heightSet1 = theSet1.cell().height();
		const uint heightSet2 = theSet2.cell().height();
		const uint maxPCHeight = ( heightSet1 <  heightSet2 ) ? heightSet2 : heightSet1;
		
		//Create the resulting GridTreeSet
		GridTreeSet resultSet( theSet1.grid(), maxPCHeight, new BinaryTreeNode() );
		
		//Adjoin the sets
		resultSet.adjoin( theSet1 );
		resultSet.adjoin( theSet2 );
		
		//Return the result
		return resultSet;
	}
	
	GridTreeSet intersection( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 ) {
		//Test that the Grids are equal
		ARIADNE_ASSERT( theSet1.grid() == theSet2.grid() );
		
		//Compute the highest primary cell 
		const uint heightSet1 = theSet1.cell().height();
		const uint heightSet2 = theSet2.cell().height();
		const uint maxPCHeight = ( heightSet1 <  heightSet2 ) ? heightSet2 : heightSet1;
		
		//Create the resulting GridTreeSet
		GridTreeSet resultSet( theSet1.grid(), maxPCHeight, new BinaryTreeNode() );
		
		//Adjoin the first set
		resultSet.adjoin( theSet1 );
		//Intersect the result with the second set
		resultSet.restrict( theSet2 );
		
		//Return the result
		return resultSet;
	}
	
	GridTreeSet difference( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 ) {
		//Test that the Grids are equal
		ARIADNE_ASSERT( theSet1.grid() == theSet2.grid() );
		
		//Compute the highest primary cell 
		const uint heightSet1 = theSet1.cell().height();
		const uint heightSet2 = theSet2.cell().height();
		const uint maxPCHeight = ( heightSet1 <  heightSet2 ) ? heightSet2 : heightSet1;
		
		//Create the resulting GridTreeSet
		GridTreeSet resultSet( theSet1.grid(), maxPCHeight, new BinaryTreeNode() );
		
		//Adjoin the first set
		resultSet.adjoin( theSet1 );
		//Remove the second set from the result set
		resultSet.remove( theSet2 );
		
		//Return the result
		return resultSet;
	}

} // namespace Ariadne

