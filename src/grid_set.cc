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
        throw IncompatibleSizes(ARIADNE_PRETTY_FUNCTION);
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
               indeterminate( otherNode._isEnabled ) ) )            &&
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
    
void BinaryTreeNode::remove( BinaryTreeNode * pThisNode, const BinaryTreeNode * pOtherNode ) {
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
            throw std::runtime_error( ARIADNE_PRETTY_FUNCTION );
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
        throw std::runtime_error( ARIADNE_PRETTY_FUNCTION );
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

bool BinaryTreeNode::overlap( const BinaryTreeNode * pRootNodeOne, const BinaryTreeNode * pRootNodeTwo ) {
    bool result = false;

    bool isNodeOneALeaf = pRootNodeOne->is_leaf();
    bool isNodeTwoALeaf = pRootNodeTwo->is_leaf();
        
    if( isNodeOneALeaf && isNodeTwoALeaf ) {
        //If both nodes are leaves, then the trees overlap if only both of the nodes are enabled
        result = pRootNodeOne->is_enabled() && pRootNodeTwo->is_enabled(); 
    } else {
        if( ! isNodeOneALeaf && isNodeTwoALeaf ){
            //If the second node is a leaf then the trees overlap is it is enabled
            //and the first node has an enabled sub-node
            result = pRootNodeTwo->is_enabled() && pRootNodeOne->has_enabled();
        } else {
            if( isNodeOneALeaf && ! isNodeTwoALeaf ){
                //If the first node is a leaf then the trees overlap is it is enabled
                //and the second node has an enabled sub-node
                result = pRootNodeOne->is_enabled() && pRootNodeTwo->has_enabled();
            } else {
                //Both nodes are non-lead nodes, then the trees overlap if
                //either their left or right branches overlap
                result = overlap( pRootNodeOne->left_node(), pRootNodeTwo->left_node() ) || 
                    overlap( pRootNodeOne->right_node(), pRootNodeTwo->right_node() );
            }
        }
    }
        
    return result;
}

bool BinaryTreeNode::subset( const BinaryTreeNode * pRootNodeOne, const BinaryTreeNode * pRootNodeTwo ) {
    bool result = false;

    bool isNodeOneALeaf = pRootNodeOne->is_leaf();
    bool isNodeTwoALeaf = pRootNodeTwo->is_leaf();
        
    if( isNodeOneALeaf && isNodeTwoALeaf ) {
        //If both nodes are leaves, then pRootNodeOne is a subset of pRootNodeTwo if:
        // 1. both of the nodes are enabled or
        // 2. pRootNodeOne is disabled (represents an empty set)
        result = ( ! pRootNodeOne->is_enabled() ) || pRootNodeTwo->is_enabled(); 
    } else {
        if( ! isNodeOneALeaf && isNodeTwoALeaf ){
            //If pRootNodeTwo is a leaf then pRootNodeOne is a subset of pRootNodeTwo if:
            // 1. pRootNodeTwo is enabled or
            // 2. pRootNodeOne has no enabled sub-nodes
            result = pRootNodeTwo->is_enabled() || ( ! pRootNodeOne->has_enabled() );
        } else {
            if( isNodeOneALeaf && ! isNodeTwoALeaf ){
                //If pRootNodeOne is a leaf then pRootNodeOne is a subset of pRootNodeTwo if:
                // 1. pRootNodeOne is disabled or
                // 2. all of the pRootNodeTwo's leaf nodes are enabled
                result = ( ! pRootNodeOne->is_enabled() ) || pRootNodeTwo->all_enabled();
            } else {
                //Both nodes are non-lead nodes, then pRootNodeOne is a subset of pRootNodeTwo
                //if sub-trees of pRootNodeOne are subsets of the subtrees of pRootNodeTwo
                result = subset( pRootNodeOne->left_node(), pRootNodeTwo->left_node() )
                         &&
                         subset( pRootNodeOne->right_node(), pRootNodeTwo->right_node() );
            }
        }
    }
        
    return result;
}

/********************************************GridTreeCursor***************************************/

/****************************************GridTreeConstIterator************************************/

GridTreeConstIterator::GridTreeConstIterator(  ) : _pGridTreeCursor()  {  
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
    return    this->_theGrid == other._theGrid && 
        this->_theHeight == other._theHeight &&
        this->_theWord == other._theWord;
}      

// PIETER: TODO: Ivan, please take a look at operator< to make sure it's
// independant of the presentation of the cell.
bool GridCell::operator<(const GridCell& other) const {
    ARIADNE_ASSERT( this->_theGrid == other._theGrid );
    return    this->_theHeight < other._theHeight || 
        ( this->_theHeight == other._theHeight && this->_theWord < other._theWord );
}      

Vector<Interval> GridCell::primary_cell_lattice_box( const uint theHeight, const dimension_type dimensions ) {
    int leftBottomCorner = 0, rightTopCorner = 1;
    //The zero level coordinates are known, so we need to iterate only for higher level primary cells
    for(uint i = 1; i <= theHeight; i++){
        primary_cell_at_height(i, leftBottomCorner, rightTopCorner);
    }
    // 1.2 Constructing and return the box defining the primary cell (relative to the grid).
        
    return Vector< Interval >( dimensions, Interval( leftBottomCorner, rightTopCorner ) );
}

uint GridCell::smallest_enclosing_primary_cell_height( const Vector<Interval>& theLatticeBox ) {
    const dimension_type dimensions = theLatticeBox.size();
    int leftBottomCorner = 0, rightTopCorner = 1;
    uint height = 0;
    //The zero level coordinates are known, so we need to iterate only for higher level primary cells
    do{
        //Check if the given box is a subset of a primary cell
        Vector< Interval > primaryCellBox( dimensions, Interval( leftBottomCorner, rightTopCorner ) );
        if(  inside(theLatticeBox, primaryCellBox ) ) {
            //If yes then we are done
            break;
        }
        //Otherwise increase the height and recompute the new borders
        primary_cell_at_height( ++height, leftBottomCorner, rightTopCorner);
    }while( true );
        
    return height;
}

/*! \brief Apply grid data \a theGrid to \a theLatticeBox in order to compute the box dimensions in the original space*/
Box GridCell::lattice_box_to_space(const Vector<Interval> & theLatticeBox, const Grid& theGrid ){
    const uint dimensions = theGrid.dimension();
    Box theTmpBox( dimensions );
        
    Vector<Float> theGridOrigin( theGrid.origin() );
    Vector<Float> theGridLengths( theGrid.lengths() );
        
    for(uint current_dimension = 0; current_dimension < dimensions; current_dimension ++){
        const Float theDimLength = theGridLengths[current_dimension];
        const Float theDimOrigin = theGridOrigin[current_dimension];
        //Recompute the new dimension coordinates, detaching them from the grid 
        theTmpBox[current_dimension].set_lower( add_approx( theDimOrigin, mul_approx( theDimLength, theLatticeBox[current_dimension].lower() ) ) );
        theTmpBox[current_dimension].set_upper( add_approx( theDimOrigin, mul_approx( theDimLength, theLatticeBox[current_dimension].upper() ) ) );
    }
        
    return theTmpBox;
}

//The box is not related to the Grid, whereas the binary tree is related to the Lattice of the grid:
// 1. Compute the primary cell located the the height \a theHeight above the zero level,
// 2. Compute the cell defined by the path \a theWord (from the primary cell).
// 3. Use Grid data to compute the box coordinates in the original space.
Box GridCell::compute_box(const Grid& theGrid, const uint theHeight, const BinaryWord& theWord) {
    //1. Obtain the primary-cell box, related to some grid
    const uint dimensions = theGrid.dimension();
    Vector<Interval>  theTmpLatticeBox( primary_cell_lattice_box( theHeight , dimensions ) );

    //2. Compute the cell on some grid, corresponding to the binary path from the primary cell.
    uint current_dimension = 0;
    for(uint i = 0; i < theWord.size(); i++){
        //We move through the dimensions in a linear fasshion
        current_dimension = i % dimensions;
        //Compute the middle point of the box's projection onto
        //the dimension \a current_dimension (relative to the grid)
        Float middlePointInCurrDim = theTmpLatticeBox[current_dimension].midpoint();
        if( theWord[i] ){
            //Choose the right half
            theTmpLatticeBox[current_dimension].set_lower( middlePointInCurrDim );
        } else {
            //Choose the left half
            theTmpLatticeBox[current_dimension].set_upper( middlePointInCurrDim );
        }
    }
        
    // 3. Use Grid data to compute the box coordinates in the original space.
    return lattice_box_to_space( theTmpLatticeBox, theGrid );
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

tribool GridTreeSubset::covers( const BinaryTreeNode* pCurrentNode, const Grid& theGrid,
                                const uint theHeight, BinaryWord &theWord, const Box& theBox ) {
    tribool result;
    
    //Check if the current node's cell intersects with theBox
    Box theCellsBox = GridCell::compute_box( theGrid, theHeight, theWord );
    tribool doIntersect = theCellsBox.overlaps( theBox );
    
    if( ! doIntersect ) {
        //If theBox does not intersect with the cell then for the covering relation
        //is is not important if we add or remove this cell, so we return true
        result = true;
    } else {
        //If the cell possibly intersects with theBox then
        if( pCurrentNode->is_leaf() ) {
            if( pCurrentNode->is_enabled() ){
                //If this is an enabled node that possibly or definitely intersects
                //with theBox then the covering property is fine, so we return true
                result = true;
            } else {
                //If the node is disabled then if it definitely intersects with theBox
                //we have to report false, as there is not covering, in case it is only
                //possible intersecting with theBox then we return possibly. The latter
                //is because we are not completely sure, if the intersection does not
                //have place then the covering property is not broken.
                result = ! doIntersect;
            }
        } else {
            //The node is not a leaf so we need to go down and see if the cell
            //falls into sub cells for which we can sort things out
            theWord.push_back(false);
            const tribool result_left = covers( pCurrentNode->left_node(), theGrid, theHeight, theWord, theBox );
            theWord.pop_back();
            
            if( ! result_left) {
                //If there is definitely no covering property then this is all
                //we need to know so there is not need to check the other branch.
                result = false;
            } else {
                //If the covering property holds or is possible, then we still
                //need to check the second branch because it can change the outcome.
                theWord.push_back(true);
                const tribool result_right = covers( pCurrentNode->right_node(), theGrid, theHeight, theWord, theBox );
                theWord.pop_back();
                
                if( !result_right ) {
                    //IF: The right sub-node reports false, then the result is false
                    //NOTE: We already sorted out the case of ( (!result_left) == true) before
                    result = false;
                } else {
                    if( result_left && result_right ) {
                        //ELSE: If both sub-nodes report true then it is true
                        result = true;
                    } else {
                        if( indeterminate(result_left) || indeterminate(result_right) ) {
                            //ELSE: If one of the sub-nodes reports indeterminate then indeterminate,
                            result = indeterminate;
                        } else {
                            //ELSE: An impossible situation
                            ARIADNE_ASSERT( false );
                        }
                    }
                }
            }
        }
    }
    
    return result;
}

tribool GridTreeSubset::subset( const BinaryTreeNode* pCurrentNode, const Grid& theGrid,
                                const uint theHeight, BinaryWord &theWord, const Box& theBox ) {
    tribool result;
    
    //Check if the current node overlaps with theBox
    Box theCellsBox = GridCell::compute_box( theGrid, theHeight, theWord );
    tribool isASubset = theCellsBox.subset( theBox );
    
    if( isASubset ){
        //It does not matter if pCurrentNode has enableds leaves or not we already know that the cell
        //corresponding to this node is geometrically a subset of theBox.
        result = true;
    } else {
        //If the cell corresponding to pCurrentNode is not a subset, then
        if( pCurrentNode->is_leaf() && ! isASubset ) {
            //If pCurrentNode is a leaf node and geometrically the cell (corresponding to the node theCellsBox) is not a
            //subset of theBox, then: if it is enabled then pCurrentNode is not a subset of theBox but otherwise it is.
            result = ! pCurrentNode->is_enabled();
        } else {
            if( pCurrentNode->is_leaf() && indeterminate( isASubset ) ) {
                //If we are in a leaf node but we do not know for sure if the given cell
                //is a subset of theBox then we can only check if it is enabled or not
                if( pCurrentNode->is_enabled() ){
                    //For an enabled-leaf node (a filled cell) we do not know if it is a subset of theBox
                    result = indeterminate;
                } else {
                    //The node is disabled, so it represents an empty set, which is a subset of any set
                    result = true;
                }
            } else {
                //The node is not a leaf, and we either know that the cell of pCurrentNode is not a geometrical subset
                //of theBox or we are not sure that it is, This means that we can do recursion to sort things out.
                theWord.push_back(false);
                const tribool result_left = subset( pCurrentNode->left_node(), theGrid, theHeight, theWord, theBox );
                theWord.pop_back();
                
                if( !result_left ) {
                    //If the left branch is not a subset, then there is no need to check the right one
                    result = false;
                } else {
                    //if we still do not know the answer, then we check the right branch
                    theWord.push_back(true);
                    const tribool result_right = subset( pCurrentNode->right_node(), theGrid, theHeight, theWord, theBox );
                    theWord.pop_back();
                    
                    if( !result_right ) {
                        //IF: The right sub-node reports false, then the result is false
                        //NOTE: We already sorted out the case of ( (!result_left) == true) before
                        result = false;
                    } else {
                        if( result_left && result_right ) {
                            //ELSE: If both sub-nodes report true then it is true
                            result = true;
                        } else {
                            if( indeterminate(result_left) || indeterminate(result_right) ) {
                                //ELSE: If one of the sub-nodes reports indeterminate then indeterminate,
                                result = indeterminate;
                            } else {
                                //ELSE: An impossible situation
                                ARIADNE_ASSERT( false );
                            }
                        }
                    }
                }
            }
        }
    }
    
    return result;
}

tribool GridTreeSubset::disjoint( const BinaryTreeNode* pCurrentNode, const Grid& theGrid,
                                  const uint theHeight, BinaryWord &theWord, const Box& theBox ) {
    tribool intersect;
    
    //Check if the current node overlaps with theBox
    Box theCellsBox = GridCell::compute_box( theGrid, theHeight, theWord );
    tribool doPossiblyIntersect = !theCellsBox.disjoint( theBox );
    
    if( doPossiblyIntersect || indeterminate( doPossiblyIntersect ) ) {
        //If there is a possible intersection then we do the checking
        if( pCurrentNode->is_leaf() ) {
            //If this is a leaf node then
            if( pCurrentNode->is_enabled() ){
                //If the node is enabled, then we have a possible intersection
                intersect = doPossiblyIntersect;
            } else {
                //Since the node is disabled, there can be no intersection
                intersect = false;
            }
        } else {
            //The node is not a leaf and the intersection is possible so check the left sub-node
            theWord.push_back(false);
            const tribool intersect_left = overlaps( pCurrentNode->left_node(), theGrid, theHeight, theWord, theBox );
            theWord.pop_back();
            
            //
            //WARNING: I know how to write a shorter code, like:
            //          if( ! ( intersect = definitely( intersect_left ) ) ) {
            //              ...
            //          }
            // I DO NOT DO THIS ON PERPOSE, KEEP THE CODE EASILY UNDERSTANDABLE!
            //
            if( intersect_left ) {
                //If we definitely have intersection for the left branch then answer is true
                intersect = true;
            } else {
                //If we still not sure/ or do not know then try to search further, i.e. check the right node
                theWord.push_back(true);
                const tribool intersect_right = overlaps( pCurrentNode->right_node(), theGrid, theHeight, theWord, theBox );
                theWord.pop_back();
                if( intersect_right ) {
                    //If we definitely have intersection for the right branch then answer is true
                    intersect = true;
                } else {
                    //Now either we have indeterminate answers or we do not know, if one
                    //if the answers for one of the branches was indeterminate then we
                    //report indeterminate, otherwise it is definitely false.
                    if( indeterminate( intersect_left ) || indeterminate( intersect_right )  ) {
                        intersect = indeterminate;
                    } else {
                        intersect = false;
                    }
                    //ERROR: Substituting the above if statement with the following conditional assignement DOES NOT WORK:
                    //    intersect = ( ( indeterminate( intersect_left ) || indeterminate( intersect_right ) ) ? indeterminate : false );
                    //In this case, if we need to assign false, the proper branch of the conditional statement is executed,
                    //but somehow the value of the "intersect" variable becomes indeterminate!
                }
            }
        }
    } else {
        //If there is no intersection then we just stop with a negative intersect
        intersect = false;
    }
    
    return !intersect;
}

tribool GridTreeSubset::overlaps( const BinaryTreeNode* pCurrentNode, const Grid& theGrid,
                                  const uint theHeight, BinaryWord &theWord, const Box& theBox ) {
    tribool result;
    
    //Check if the current node overlaps with theBox
    Box theCellsBox = GridCell::compute_box( theGrid, theHeight, theWord );
    tribool doPossiblyIntersect = theCellsBox.overlaps( theBox );
    
    if( doPossiblyIntersect || indeterminate( doPossiblyIntersect ) ) {
        //If there is a possible intersection then we do the checking
        if( pCurrentNode->is_leaf() ) {
            //If this is a leaf node then
            if( pCurrentNode->is_enabled() ){
                //If the node is enabled, then we have a possible intersection
                result = doPossiblyIntersect;
            } else {
                //Since the node is disabled, there can be no intersection
                result = false;
            }
        } else {
            //The node is not a leaf and the intersection is possible so check the left sub-node
            theWord.push_back(false);
            const tribool result_left = overlaps( pCurrentNode->left_node(), theGrid, theHeight, theWord, theBox );
            theWord.pop_back();
            
            //
            //WARNING: I know how to write a shorter code, like:
            //          if( ! ( result = definitely( result_left ) ) ) {
            //              ...
            //          }
            // I DO NOT DO THIS ON PERPOSE, KEEP THE CODE EASILY UNDERSTANDABLE!
            //
            if( result_left ) {
                //If we definitely have intersection for the left branch then answer is true
                result = true;
            } else {
                //If we still not sure/ or do not know then try to search further, i.e. check the right node
                theWord.push_back(true);
                const tribool result_right = overlaps( pCurrentNode->right_node(), theGrid, theHeight, theWord, theBox );
                theWord.pop_back();
                if( result_right ) {
                    //If we definitely have intersection for the right branch then answer is true
                    result = true;
                } else {
                    //Now either we have indeterminate answers or we do not know, if one
                    //if the answers for one of the branches was indeterminate then we
                    //report indeterminate, otherwise it is definitely false.
                    if( indeterminate( result_left ) || indeterminate( result_right )  ) {
                        result = indeterminate;
                    } else {
                        result = false;
                    }
                    //ERROR: Substituting the above if statement with the following conditional assignement DOES NOT WORK:
                    //    result = ( ( indeterminate( result_left ) || indeterminate( result_right ) ) ? indeterminate : false );
                    //In this case, if we need to assign false, the proper branch of the conditional statement is executed,
                    //but somehow the value of the "result" variable becomes indeterminate!
                }
            }
        }
    } else {
        //If there is no intersection then we just stop with a negative result
        result = false;
    }
    
    return result;
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

BinaryTreeNode* GridTreeSet::align_with_cell( const uint otherPavingPCellHeight, const bool stop_on_enabled,
                                              const bool stop_on_disabled, bool & has_stopped ) {
    const uint thisPavingPCellHeight = this->cell().height();
        
    //The current root node of the GridTreeSet
    BinaryTreeNode * pBinaryTreeNode = this->_pRootTreeNode;
        
    if( thisPavingPCellHeight > otherPavingPCellHeight ){

        //The primary cell of this paving is higher then the one of the other paving.
        //1. We locate the path to the primary cell node common with the other paving
        BinaryWord primaryCellPath = GridCell::primary_cell_path( this->cell().grid().dimension(), thisPavingPCellHeight, otherPavingPCellHeight );
            
        //2. Locate the binary tree node corresponding to this primnary cell
        uint position = 0;
        while( position < primaryCellPath.size() &&
               !( has_stopped = ( ( pBinaryTreeNode->is_enabled() && stop_on_enabled ) ||
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
    
void GridTreeSet::_adjoin_outer_approximation( const Grid & theGrid, BinaryTreeNode * pBinaryTreeNode, const uint primary_cell_height,
                                               const uint max_mince_depth,  const CompactSetInterface& theSet, BinaryWord * pPath ){
    //Compute the cell correspomding to the current node
    GridCell theCurrentCell( theGrid, primary_cell_height, *pPath );

    const OpenSetInterface* pOpenSet=dynamic_cast<const OpenSetInterface*>(static_cast<const SetInterfaceBase*>(&theSet));
    
    if( bool( theSet.disjoint( theCurrentCell.box() ) ) ) {
        //DO NOTHING: We are in the node whoes representation in the original space is
        //disjoint from pSet and thus there will be nothing added to this cell.
    } else if( pOpenSet && bool( pOpenSet->covers( theCurrentCell.box() ) ) ) {
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
                _adjoin_outer_approximation( theGrid, pBinaryTreeNode->left_node(), primary_cell_height, max_mince_depth, theSet, pPath );
                //Check the right branch
                pPath->push_back(true);
                _adjoin_outer_approximation( theGrid, pBinaryTreeNode->right_node(), primary_cell_height, max_mince_depth, theSet, pPath );
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

// FIXME: This method can fail if we cannot determine which of a node's children overlaps
// the set. In principle this can be solved by checking if one of the children overlaps
// the set, before doing recursion, and if none overlaps then we mark the present node as
// enabled and stop. Generally speaking, the present algorithm is not wrong, it also gives us
// a lower approximation, but it is simply less accurate than it could be.
// TODO:Think of another representation in terms of covers but not pavings, then this problem
// can be cured in a different fashion.
void GridTreeSet::_adjoin_lower_approximation( const Grid & theGrid, BinaryTreeNode * pBinaryTreeNode, const uint primary_cell_height,
                                               const uint max_mince_depth,  const OvertSetInterface& theSet, BinaryWord * pPath ){
    //Compute the cell correspomding to the current node
    GridCell theCurrentCell( theGrid, primary_cell_height, *pPath );

    if( bool( theSet.overlaps( theCurrentCell.box() ) ) ) {
        if( pPath->size() >= max_mince_depth ) {
            //We should not mince any further. If the cell is not a leaf, then some subset is enabled, so the
            //lower approximation does not add any information. If the cell is a leaf, we mark it as enabled.
            if( ! pBinaryTreeNode->has_enabled() ) {
                pBinaryTreeNode->make_leaf(true);
            }
        } else {
            //If the node is no enabled, so may be we can add something from the lower approximation of theSet.
            //Since we still do not have the finest cells for the lower approximation of theSet, we split 
            pBinaryTreeNode->split(); //NOTE: splitting a non-leaf node does not do any harm
            //Check the left branch
            pPath->push_back(false);
            _adjoin_lower_approximation( theGrid, pBinaryTreeNode->left_node(), primary_cell_height, max_mince_depth, theSet, pPath );
            //Check the right branch
            pPath->push_back(true);
            _adjoin_lower_approximation( theGrid, pBinaryTreeNode->right_node(), primary_cell_height, max_mince_depth, theSet, pPath );
        }
    }
    //Return to the previous level, since the initial call is made
    //with the empty word, we check that it is not yet empty.
    if( pPath->size() > 0 ) {
        pPath->pop_back();
    }
}
    
void GridTreeSet::_adjoin_lower_approximation( const Grid & theGrid, BinaryTreeNode * pBinaryTreeNode, const uint primary_cell_height,
                                               const uint max_mince_depth,  const OpenSetInterface& theSet, BinaryWord * pPath ){
    //Compute the cell corresponding to the current node
    GridCell theCurrentCell( theGrid, primary_cell_height, *pPath );
    
    if( bool( theSet.covers( theCurrentCell.box() ) ) ) {
        pBinaryTreeNode->make_leaf(true);
        pBinaryTreeNode->mince( max_mince_depth - pPath->size() );
    } else if ( bool( theSet.overlaps( theCurrentCell.box() ) ) ) {
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
            pBinaryTreeNode->split();  //NOTE: splitting a non-leaf node does not do any harm
            
            //Check the left branch
            pPath->push_back(false);
            _adjoin_lower_approximation( theGrid, pBinaryTreeNode->left_node(), primary_cell_height, max_mince_depth, theSet, pPath );
            //Check the right branch
            pPath->push_back(true);
            _adjoin_lower_approximation( theGrid, pBinaryTreeNode->right_node(), primary_cell_height, max_mince_depth, theSet, pPath );
        }        
    }
    //Return to the previous level, since the initial call is made
    //with the empty word, we check that it is not yet empty.
    if( pPath->size() > 0 ) {
        pPath->pop_back();
    }
}
    
void GridTreeSet::adjoin_over_approximation( const Box& theBox, const uint depth ) {
    // FIXME: This adjoins an outer approximation; change to ensure only overlapping cells are adjoined
    for(size_t i=0; i!=theBox.dimension(); ++i) {
        if(theBox[i].lower()>=theBox[i].upper()) {
            ARIADNE_THROW(std::runtime_error,"GridTreeSet::adjoin_over_approximation(Box,uint)","Box "<<theBox<<" has empty interior.");
        }
    }
    this->adjoin_outer_approximation(theBox,depth);
}

void GridTreeSet::adjoin_outer_approximation( const CompactSetInterface& theSet, const uint depth ) {
    Grid theGrid( this->cell().grid() );
    ARIADNE_ASSERT( theSet.dimension() == this->cell().dimension() );
        
    //1. Compute the smallest GridCell (corresponding to the primary cell)
    //   that encloses the theSet (after it is mapped onto theGrid).
    const GridCell theOverApproxCell = GridCell::smallest_enclosing_primary_cell( theSet.bounding_box(), theGrid );
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
        _adjoin_outer_approximation( GridTreeSubset::_theGridCell.grid(), pBinaryTreeNode, outer_approx_primary_cell_height,
                                     max_mince_depth, theSet, pEmptyPath );

        delete pEmptyPath;
    }
}

// TODO:Think of another representation in terms of covers but not pavings, then the implementation
// will be different, this is why, for now we do not fix these things.
void GridTreeSet::adjoin_lower_approximation( const LocatedSetInterface& theSet, const uint depth ) {
    this->adjoin_lower_approximation( theSet, theSet.bounding_box(), depth );
}

void GridTreeSet::adjoin_lower_approximation( const OvertSetInterface& theSet, const uint height, const uint depth ) {
    Grid theGrid( this->cell().grid() );
    ARIADNE_ASSERT( theSet.dimension() == this->cell().dimension() );
    
    //Align this paving and paving at the given height
    bool has_stopped = false;
    BinaryTreeNode* pBinaryTreeNode = align_with_cell( height, true, false, has_stopped );
        
    //If the lower aproximation of the bounding box of the provided set is enclosed
    //in an enabled cell of this paving, then there is nothing to be done. The latter
    //is because adjoining the outer approx of the set will not change this paving.
    if( ! has_stopped ){
        //Compute the depth to which we must mince the outer approximation of the adjoining set.
        //This depth is relative to the root of the constructed paving, which has been alligned
        //with the binary tree node pBinaryTreeNode.
        const uint max_mince_depth = height * theGrid.dimension() + depth;
            
        //Adjoin the outer approximation, computing it on the fly.
        BinaryWord * pEmptyPath = new BinaryWord(); 
        //const RegularSetInterface* theRegularVersionOfSet = dynamic_cast<const RegularSetInterface*>(&theSet);
        const OpenSetInterface* theOpenVersionOfSet = dynamic_cast<const OpenSetInterface*>(&theSet);
        //const LocatedSetInterface* theLocatedVersionOfSet = dynamic_cast<const LocatedSetInterface*>(&theSet);
        const OvertSetInterface* theOvertVersionOfSet = dynamic_cast<const OvertSetInterface*>(&theSet);
        if( theOpenVersionOfSet ) {
            _adjoin_lower_approximation( GridTreeSubset::_theGridCell.grid(), pBinaryTreeNode, height, max_mince_depth, *theOpenVersionOfSet, pEmptyPath );
        } else {
            _adjoin_lower_approximation( GridTreeSubset::_theGridCell.grid(), pBinaryTreeNode, height, max_mince_depth, *theOvertVersionOfSet, pEmptyPath );
        }
        delete pEmptyPath;
    }
}

void GridTreeSet::adjoin_lower_approximation( const OvertSetInterface& theSet, const Box& theBoundingBox, const uint depth ) {
    Grid theGrid( this->cell().grid() );
    ARIADNE_ASSERT( theSet.dimension() == this->cell().dimension() );
    ARIADNE_ASSERT( theBoundingBox.dimension() == this->cell().dimension() );
    
    //Compute the smallest the primary cell that encloses the theSet (after it is mapped onto theGrid).
    const GridCell theOverApproxCell = GridCell::smallest_enclosing_primary_cell( theBoundingBox, theGrid );
    
    //Adjoin the lower approximation with the bounding cell being the primary cell at the given height.
    adjoin_lower_approximation( theSet, theOverApproxCell.height(), depth );
}

void GridTreeSet::_adjoin_inner_approximation( const Grid & theGrid, BinaryTreeNode * pBinaryTreeNode, const uint primary_cell_height,
                                               const uint max_mince_depth, const OpenSetInterface& theSet, BinaryWord * pPath ) {
    //Compute the cell corresponding to the current node
    GridCell theCurrentCell( theGrid, primary_cell_height, *pPath );

    if( ! pBinaryTreeNode->is_enabled() ) {
        //If this it is not an enabled leaf node then we can add something to it.
        if( definitely( theSet.covers( theCurrentCell.box() ) ) ) {
            //If this node's box is a subset of theSet then it belongs to the inner approximation
            //Thus we need to make it an enabled leaf and then do the mincing to the maximum depth
            pBinaryTreeNode->make_leaf( true );
        } else if ( possibly( theSet.overlaps( theCurrentCell.box() ) ) ) {
            //If theSet overlaps with the box corresponding to the given node in the original
            //space, then there might be something to add from the inner approximation of theSet.
            if( pPath->size() >= max_mince_depth ) {
                //DO NOTHING: Since it is the maximum depth to which we were allowed to mince
                //and we know that the node's box only overlaps with theSet but is not it's
                //subset we have to exclude it from the inner approximation of theSet.
            } else {
                //Since we still do not have the finest cells for the outer approximation of theSet, we split 
                pBinaryTreeNode->split();  //NOTE: splitting a non-leaf node does not do any harm
                
                //Check the left branch
                pPath->push_back(false);
                _adjoin_inner_approximation( theGrid, pBinaryTreeNode->left_node(), primary_cell_height, max_mince_depth, theSet, pPath );
                //Check the right branch
                pPath->push_back(true);
                _adjoin_inner_approximation( theGrid, pBinaryTreeNode->right_node(), primary_cell_height, max_mince_depth, theSet, pPath );
            }
        } else {
            //DO NOTHING: the node's box is disjoint from theSet and thus it or it's
            //sub nodes can not be the part of the inner approximation of theSet.
        }
    } else {
        //DO NOTHING: If this is an enabled leaf node we are in, then there is nothing to add to it
        //TODO: We could mince the node until the depth (max_mince_depth - pPath->size()) but I do
        //not know if it makes any sence to do this, since this does not change the result.
    }
    
    //Return to the previous level, since the initial call is made
    //with the empty word, we check that it is not yet empty.
    if( pPath->size() > 0 ) {
        pPath->pop_back();
    }
}

void GridTreeSet::adjoin_inner_approximation( const OpenSetInterface& theSet, const uint height, const uint depth ) {
    Grid theGrid( this->cell().grid() );
    ARIADNE_ASSERT( theSet.dimension() == this->cell().dimension() );

    //Align this paving and paving at the given height
    bool has_stopped = false;
    BinaryTreeNode* pBinaryTreeNode = align_with_cell( height, true, false, has_stopped );
        
    //If the inner aproximation of the bounding box of the provided set is enclosed
    //in an enabled cell of this paving, then there is nothing to be done. The latter
    //is because adjoining the inner approx of the set will not change this paving.
    if( ! has_stopped ){
        //Compute the depth to which we must mince the inner approximation of the adjoining set.
        //This depth is relative to the root of the constructed paving, which has been alligned
        //with the binary tree node pBinaryTreeNode.
        const uint max_mince_depth = height * theGrid.dimension() + depth;
        
        //Adjoin the inner approximation, computing it on the fly.
        BinaryWord * pEmptyPath = new BinaryWord(); 
        _adjoin_inner_approximation( GridTreeSubset::_theGridCell.grid(), pBinaryTreeNode, height, max_mince_depth, theSet, pEmptyPath );
        delete pEmptyPath;
    }
}

void GridTreeSet::adjoin_inner_approximation( const OpenSetInterface& theSet, const Box& theBoundingBox, const uint depth ) {
    Grid theGrid( this->cell().grid() );
    ARIADNE_ASSERT( theSet.dimension() == this->cell().dimension() );
    ARIADNE_ASSERT( theBoundingBox.dimension() == this->cell().dimension() );
    
    //Compute the smallest the primary cell that encloses the theSet (after it is mapped onto theGrid).
    const GridCell theOverApproxCell = GridCell::smallest_enclosing_primary_cell( theBoundingBox, theGrid );
    
    //Compute the height of the smallest primary cell that encloses the box
    //Note that, since we will need to adjoin inner approximation bounded by
    //the theBoundingBox, it is enough to take this cell's height and not to
    //go higher. Remember that the inner approximations consists of the cells
    //that are subsets of the set theSet Adjoin the inner approximation with
    //the bounding cell being the primary cell at the given height.
    adjoin_inner_approximation( theSet, theOverApproxCell.height(), depth );
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
    
void GridTreeSet::clear(  ) {
    // TODO: Pieter: This implementation may not be optimal. Please could you
    // check this, Ivan.
    *this=GridTreeSet(this->grid());
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
        uint position = 0;
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
        for( uint i = 0; i < pathToPCell.size(); i++ ) {
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
    
GridCell GridCell::smallest_enclosing_primary_cell( const Box& theBox, const Grid& theGrid) {
    Vector<Interval> theLatticeBox( theBox.size() );
    //Convert the box to theGrid coordinates
    for( uint i = 0; i != theBox.size(); ++i ) {
        theLatticeBox[i] = ( theBox[i] - theGrid.origin()[i] ) / theGrid.lengths()[i];
    }
    //Compute the smallest primary cell, enclosing this grid
    uint height = GridCell::smallest_enclosing_primary_cell_height( theLatticeBox );
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
    const uint theSetsPCellHeight = theSet.cell().height();
    if( theSetsPCellHeight > theCell.height() ) {
        //Compute the path from the primary cell of theSet to the primary cell of theCell
        BinaryWord pathFromSetsPCellToCell = GridCell::primary_cell_path( theCell.dimension(), theSetsPCellHeight, theCell.height() );
        pathFromSetsPCellToCell.append( theCell.word() );
        pWorkGridCell = new GridCell( theCell.grid(), theSetsPCellHeight, pathFromSetsPCellToCell );
    } else {
        pWorkGridCell = &theCell;
    }

    //Compute the path for the primary cell of theCell to the primary cell of theSet
    BinaryWord pathFromPCellCellToSetsRootNode = GridCell::primary_cell_path( theCell.dimension(), pWorkGridCell->height(), theSetsPCellHeight );
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
            //check if it overlaps with the tree
                
            const BinaryTreeNode *pCurrentNode = theSet.binary_tree();
            //Note that, pathFromPCellCellToSetsRootNode.size() < workCellWord.size()
            //Because we already checked for workCellWord.is_prefix( pathFromPCellCellToSetsRootNode )
            //Here we try to find the node corresponding to theCell in the binary tree of theSet
            //in case we encounter a leaf node then we just stop, because it is enough information for us
            for( uint i = pathFromPCellCellToSetsRootNode.size(); i < workCellWord.size(); i++ ) {
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
        
    if( theSetsPCellHeight > theCell.height() ) {
        delete pWorkGridCell;
    }
        
    return result;
}

/*! \brief This is a helper method, it receives two GridTreeSubset elements
 *  then it computes the primary cell that is common to them in a sence that
 *  these sets can be rooted to it. After that the method updates the paths
 *  with the information about the paths from the found primary cell to the
 *  root binary tree nodes of both sets.
 */
static void common_primary_cell_path(const GridTreeSubset& theSet1, const GridTreeSubset& theSet2, BinaryWord &pathCommonPCtoRC1, BinaryWord &pathCommonPCtoRC2 ) {
    //Get the root cells for the subsets
    GridCell rootCell1 = theSet1.cell();
    GridCell rootCell2 = theSet2.cell();
    //Get the heights of the primary cells for both subsets
    const int heightPC1 = rootCell1.height();
    const int heightPC2 = rootCell2.height();
    if( heightPC2 > heightPC1 ) {
        //Compute the path from the common primary cell to the root cell of theSet1
        pathCommonPCtoRC1 = GridCell::primary_cell_path( rootCell1.dimension(), heightPC2, heightPC1 );
        pathCommonPCtoRC1.append( rootCell1.word() );
        //Compute the path from the common primary cell to the root cell of theSet2
        pathCommonPCtoRC2 = rootCell2.word();
    } else {
        //Compute the path from the common primary cell to the root cell of theSet2
        pathCommonPCtoRC1 = rootCell1.word();
        //Compute the path from the common primary cell to the root cell of theSet1
        pathCommonPCtoRC2 = GridCell::primary_cell_path( rootCell1.dimension(), heightPC1, heightPC2 );
        pathCommonPCtoRC2.append( rootCell2.word() );
    }
}

/*! \brief this method locates the node in the tree rooted to pSuperTreeRootNode that corresponds to
 *  the path pathFromSuperToSub. In case we encounter a leaf node then we stop and return the node
 */
static const BinaryTreeNode * locate_node( const BinaryTreeNode * pSuperTreeRootNode, const BinaryWord& pathFromSuperToSub ){
    //Locate the node in the pSuperTreeRootNode such that it corresponds to pSubTreeRootNode
    const BinaryTreeNode * pCurrentSuperTreeNode = pSuperTreeRootNode;
    for( uint i = 0; i < pathFromSuperToSub.size(); i++ ) {
        if( pCurrentSuperTreeNode->is_leaf() ) {
            //We are in the leaf node and we have not yet reached the node
            //corresponding to pSubTreeRootNode, because i < pathFromSuperToSub.size()
            break;
        } else {
            pCurrentSuperTreeNode = ( pathFromSuperToSub[i] ? pCurrentSuperTreeNode->right_node() : pCurrentSuperTreeNode->left_node() );
        }
    }
    return pCurrentSuperTreeNode;
}

/*! \brief This is a helper functiuon for bool subset( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 )
 *  pSuperTreeRootNode is the root tree node, pathFromSuperToSub is the path from this root node to the root node
 *  pSubTreeRootNode. In general we see if the set represented by pSubTreeRootNode is a subset of pSuperTreeRootNode.
 *  If pSubTreeRootNode has no enabled leaf nodes then the result is always true, else if pSuperTreeRootNode has no
 *  enabled leaf nodes then the result is always false.
 */
static bool subset(const BinaryTreeNode * pSubTreeRootNode, const BinaryTreeNode * pSuperTreeRootNode, const BinaryWord & pathFromSuperToSub) {
    bool result = false;

    //Check if both sets are not empty
    if( pSubTreeRootNode->has_enabled() ) {
        if( pSuperTreeRootNode->has_enabled() ) {
            //Locate the node in pSuperTreeRootNode, by following the path pathFromSuperToSub
            const BinaryTreeNode * pCurrentSuperTreeNode = locate_node( pSuperTreeRootNode, pathFromSuperToSub );
            
            if( pCurrentSuperTreeNode->is_leaf() ) {
                //If we've reached the leaf node then pSubTreeRootNode is a subset of pSuperTreeRootNode if this
                //node is enabled, otherwise not. This is because we know that pSubTreeRootNode is not empty.
                result = pCurrentSuperTreeNode->is_enabled();
            } else {
                //At this point pCurrentSuperTreeNode corresponds to pSubTreeRootNode
                //So the trees are aligned now and we can continue checking further.
                result = BinaryTreeNode::subset( pSubTreeRootNode, pCurrentSuperTreeNode );
            }
        } else {
            //Nothing is a subset of an empty set except for an empty
            //set, but we know that pSubTreeRootNode is not empty.
            result = false;
        }
    } else {
        //An empty set is a subset of any set including the empty set itself
        result = true;
    }

    return result;
}

/*! \brief This is a helper functiuon for bool subset( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 )
 *  pSuperTreeRootNode is the root tree node, pathFromSuperToSub is the path from this root node to the root node
 *  pSubTreeRootNode. In general we see if the set represented by pSuperTreeRootNode is a subset of pSubTreeRootNode.
 *  If pSuperTreeRootNode has no enabled leaf nodes then the result is always true, else if pSubTreeRootNode has no
 *  enabled leaf nodes then the result is always false.
 */
static bool subset( const BinaryTreeNode * pSuperTreeRootNode, const BinaryWord & pathFromSuperToSub, const BinaryTreeNode * pSubTreeRootNode ) {
    bool result = false;
    
    //First we iterate throught the path pathFromSuperToSub trying to reach the common node with pSubTreeRootNode
    //Since we want to know if pSuperTreeRootNode is a subset of pSubTreeRootNode, the branches of pSuperTreeRootNode
    //that we ommit traveling the path pathFromSuperToSub, should contain no enabled leaf nodes. Because otherwise 
    //pSuperTreeRootNode is not a subset of pSubTreeRootNode. Apart from this, once we encounter a leaf node on the
    //path we stop because the we can already decide if pSuperTreeRootNode is a subset of pSubTreeRootNode.
    uint path_element = 0;
    bool areExtraLeavesDisabled = true;
    const BinaryTreeNode *pCurrentSuperTreeNode = pSuperTreeRootNode;
    while( ( path_element < pathFromSuperToSub.size() ) && areExtraLeavesDisabled ) {
        if( pCurrentSuperTreeNode->is_leaf() ){
            //We ended up in a leaf node so we have to stop iterating through the path
            break;
        } else {
            if( pathFromSuperToSub[ path_element ] ) {
                //The path goes right, check if the left branch has no enabled leaves
                areExtraLeavesDisabled = ! pCurrentSuperTreeNode->left_node()->has_enabled();
                pCurrentSuperTreeNode = pCurrentSuperTreeNode->right_node();
            } else {
                //The path goes left, check if the right branch has no enabled leaves
                areExtraLeavesDisabled = ! pCurrentSuperTreeNode->right_node()->has_enabled();
                pCurrentSuperTreeNode = pCurrentSuperTreeNode->left_node();
            }
        }
        path_element++;
    }
    
    if( areExtraLeavesDisabled ) {
        //If pSuperTreeRootNode does not have enabled leaves that are outside the bounding cell of
        //pSubTreeRootNode, this means that pSuperTreeRootNode can be a subset of pSubTreeRootNode.
        //Now we have to check that:
        if( pCurrentSuperTreeNode->is_leaf() ) {
            //A) We reached a leaf node, when following the path, then
            if( pCurrentSuperTreeNode->is_enabled() ) {
                //1. If it is enabled then it depends on whether we followed the path to the end
                if( path_element < pathFromSuperToSub.size() ) {
                    //1.1 If the path was not finished then pSuperTreeRootNode is a superset of pSubTreeRootNode
                    result = false;
                } else {
                    //1.2 If the path was ended, i.e. pCurrentSuperTreeNode and pSubTreeRootNode are
                    //aligned, then we simply need to check if one tree is a subset of another:
                    //     bool subset(const BinaryTreeNode *, const BinaryTreeNode * )
                    result = BinaryTreeNode::subset( pCurrentSuperTreeNode, pSubTreeRootNode );
                }
            } else {
                //2. if it is disabled then pSuperTreeRootNode is empty and thus is a subset of theSet2
                result = true;
            }
        } else {
            //B) We are in a non-leaf node, so we definitely reached the end of the path,
            //thus proceed like in case A.1.2 (the root nodes of both trees are aligned )
            result = BinaryTreeNode::subset( pCurrentSuperTreeNode, pSubTreeRootNode );
        }
    } else {
        //If there are enabled leaf nodes in theSet1 that are outside of the bounding cell of
        //pSubTreeRootNode, then clearly pSuperTreeRootNode is not a subset of pSubTreeRootNode.
        result = false;
    }
    return result;
}

bool subset( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 ) {
    bool result = false;
        
    //Test that the Grids are equal
    ARIADNE_ASSERT( theSet1.grid() == theSet2.grid() );
    
    //Define paths for the root cells of theSet1 and theSet2 from the common primary cell
    BinaryWord pathCommonPCtoRC1;
    BinaryWord pathCommonPCtoRC2;
    //Get the paths from the common primary cell to the root nodes of the set's binary trees
    common_primary_cell_path( theSet1, theSet2, pathCommonPCtoRC1, pathCommonPCtoRC2 );
    
    //At this point we know paths from the common primary cell
    //to the root nodes of both subsets. If one of these paths is a prefix of
    //the other one, then there is a chance that theSet1 is a subset of theSet2.
    //If not, then they definitely do not overlap.
    if( pathCommonPCtoRC1.is_prefix( pathCommonPCtoRC2 ) ) {
        //In this case theSet2 is a subset of the bounding cell of theSet1. Still 
        //it is possible that theSet1 is a subset of theSet2 if all cells of theSet1
        //outside the bounding box of theSet2 are disabled cells. This we check by
        //following the path from the foor of theSet1 to the root of theSet2.
        pathCommonPCtoRC2.erase_prefix( pathCommonPCtoRC1.size() );
        result = subset( theSet1.binary_tree(), pathCommonPCtoRC2, theSet2.binary_tree() );
    } else {
        if( pathCommonPCtoRC2.is_prefix( pathCommonPCtoRC1 ) ) {
            //Since pathCommonPCtoRC2 is a prefix of pathCommonPCtoRC1,
            //theSet1 can be a subset of theSet2. This is because theSet1
            //lies within the bounding cell of theSet2
            pathCommonPCtoRC1.erase_prefix( pathCommonPCtoRC2.size() );
            result = subset( theSet1.binary_tree(), theSet2.binary_tree(), pathCommonPCtoRC1 );
        } else {
            //theSet1 is a definitely not a subset of theSet2 Since their bounding boxes
            //are disjoint, due to paths (from the common primary cell to their root nodes)
            //that have different suffixes.
            result = false;
        }
    }

    return result;
}
    
/*! \brief This is a helper functiuon for bool overlap( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 )
 *  pSuperTreeRootNode is the root tree node, pathFromSuperToSub is the path from this root node to the root node
 *  pSubTreeRootNode. In general we see if the sets represented by pSuperTreeRootNode and pSubTreeRootNode overlap.
 *  If at least one of pSuperTreeRootNode and pSubTreeRootNode have no enabled leaf nodes then the result is always false.
 */
static bool overlap(const BinaryTreeNode * pSuperTreeRootNode, const BinaryWord & pathFromSuperToSub, const BinaryTreeNode * pSubTreeRootNode) {
    bool result = false;
        
    //Check if both sets are not empty
    if( pSuperTreeRootNode->has_enabled() && pSubTreeRootNode->has_enabled() ) {

        //Locate the node in pSuperTreeRootNode, by following the path pathFromSuperToSub
        const BinaryTreeNode * pCurrentSuperTreeNode = locate_node( pSuperTreeRootNode, pathFromSuperToSub );
            
        if( pCurrentSuperTreeNode->is_leaf() ) {
            //If we've reached the leaf node then the sets overlap if this node is enabled,
            //otherwise not. This is because we know that pSubTreeRootNode is not empty.
            result = pCurrentSuperTreeNode->is_enabled();
        } else {
            //At this point pCurrentSuperTreeNode corresponds to pSubTreeRootNode
            //So the trees are aligned now and we can continue checking further.
            result = BinaryTreeNode::overlap( pCurrentSuperTreeNode, pSubTreeRootNode );
        }
    }
        
    return result;
}
    
bool overlap( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 ) {
    bool result = false;
        
    //Test that the Grids are equal
    ARIADNE_ASSERT( theSet1.grid() == theSet2.grid() );
    
    //Define paths for the root cells of theSet1 and theSet2 from the common primary cell
    BinaryWord pathCommonPCtoRC1;
    BinaryWord pathCommonPCtoRC2;
    //Get the paths from the common primary cell to the root nodes of the set's binary trees
    common_primary_cell_path( theSet1, theSet2, pathCommonPCtoRC1, pathCommonPCtoRC2 );

    //At this point we know paths from the common primary cell
    //to the root nodes of both subsets. If one of these paths is a prefix of
    //the other one, then there is a chance that the subsets overlap.
    //If not, then they definitely do not overlap.
    if( pathCommonPCtoRC1.is_prefix( pathCommonPCtoRC2 ) ){
        //theSet2 is located somewhere within the bounding box of theSet1
        pathCommonPCtoRC2.erase_prefix( pathCommonPCtoRC1.size() );
        result = overlap( theSet1.binary_tree(), pathCommonPCtoRC2, theSet2.binary_tree() );
    } else {
        if( pathCommonPCtoRC2.is_prefix( pathCommonPCtoRC1 ) ){
            //theSet1 is located somewhere within the bounding box of theSet2
            pathCommonPCtoRC1.erase_prefix( pathCommonPCtoRC2.size() );
            result = overlap( theSet2.binary_tree(), pathCommonPCtoRC1, theSet1.binary_tree() );
        } else {
            //The sets do not overlap
            result = false;
        }
    }
        
    return result;
}

/*************************************FRIENDS OF GridTreeSet*****************************************/

GridTreeSet outer_approximation(const Box& theBox, const Grid& theGrid, const uint depth) {
    return outer_approximation(ImageSet(theBox),theGrid,depth);
}

GridTreeSet outer_approximation(const Box& theBox, const uint depth) {
    return outer_approximation(theBox, Grid(theBox.dimension()), depth);
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

