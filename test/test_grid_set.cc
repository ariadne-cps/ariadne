/***************************************************************************
 *            test_grid_set.cc
 *
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
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "stlio.h"

#include "macros.h"
#include "tribool.h"
#include "numeric.h"

#include "function_set.h"
#include "grid_set.h"
#include "set_interface.h"


#include "test.h"

using namespace Ariadne;
using namespace std;

void test_binary_tree() {
		
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Allocate an ebabled node and manipulate it");

	BinaryTreeNode theBinaryTreeRoot(true);
	ARIADNE_PRINT_TEST_COMMENT("Initial tree node");
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.is_leaf(), true );
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.is_enabled(), true );

	theBinaryTreeRoot.set_disabled();
	ARIADNE_PRINT_TEST_COMMENT("Disable the node");
	BinaryTreeNode expected_node1(false);
	ARIADNE_TEST_COMPARE( expected_node1, ==, theBinaryTreeRoot );
	
	ARIADNE_PRINT_TEST_COMMENT("Making the leaf node intermediate should cause the IsALeafNodeException exception");
	ARIADNE_TEST_THROW( theBinaryTreeRoot.set_unknown(), IsALeafNodeException);

	ARIADNE_PRINT_TEST_COMMENT("The node still have to be disabled");
	ARIADNE_TEST_COMPARE( expected_node1, ==, theBinaryTreeRoot );
	
	ARIADNE_PRINT_TEST_COMMENT("Enable the node");
	theBinaryTreeRoot.set_enabled();
	BinaryTreeNode expected_node2(true);
	ARIADNE_TEST_COMPARE( expected_node2, ==, theBinaryTreeRoot );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Split the enabled node");
	ARIADNE_PRINT_TEST_COMMENT("Should mark the node as intermediate and create two enabled subnodes");
	theBinaryTreeRoot.split();
	ARIADNE_TEST_COMPARE( theBinaryTreeRoot.left_node(), !=, NULL );
	ARIADNE_TEST_COMPARE( theBinaryTreeRoot.right_node(), !=, NULL );
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.is_enabled(), false );
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.is_disabled(), false );
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.left_node()->is_leaf(), true );
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.right_node()->is_leaf(), true );
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.left_node()->is_enabled(), true );
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.right_node()->is_enabled(), true );

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Split a splited node");
	ARIADNE_PRINT_TEST_COMMENT("We expect an unchanged tree, since it is already split");
	theBinaryTreeRoot.mince(1);
	ARIADNE_TEST_COMPARE( theBinaryTreeRoot.left_node(), !=, NULL );
	ARIADNE_TEST_COMPARE( theBinaryTreeRoot.right_node(), !=, NULL );
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.is_enabled(), false );
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.is_disabled(), false );
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.left_node()->is_leaf(), true );
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.right_node()->is_leaf(), true );
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.left_node()->is_enabled(), true );
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.right_node()->is_enabled(), true );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Copy the node by using the copy constructor");
	ARIADNE_PRINT_TEST_COMMENT("The entire subtree should be copied");
	BinaryTreeNode theBinaryTreeCopy( theBinaryTreeRoot );
	ARIADNE_TEST_COMPARE( theBinaryTreeCopy.left_node(), !=, NULL );
	ARIADNE_TEST_COMPARE( theBinaryTreeCopy.right_node(), !=, NULL );
	ARIADNE_TEST_EQUAL( theBinaryTreeCopy.is_enabled(), false );
	ARIADNE_TEST_EQUAL( theBinaryTreeCopy.is_disabled(), false );
	ARIADNE_TEST_EQUAL( theBinaryTreeCopy.left_node()->is_leaf(), true );
	ARIADNE_TEST_EQUAL( theBinaryTreeCopy.right_node()->is_leaf(), true );
	ARIADNE_TEST_EQUAL( theBinaryTreeCopy.left_node()->is_enabled(), true );
	ARIADNE_TEST_EQUAL( theBinaryTreeCopy.right_node()->is_enabled(), true );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Recombine the new tree and see what happens");
	ARIADNE_PRINT_TEST_COMMENT("The new tree should be reduced to a node");
	ARIADNE_PRINT_TEST_COMMENT("The old tree should remain unchanged");
	theBinaryTreeCopy.recombine();
	ARIADNE_PRINT_TEST_COMMENT( "Reduced tree" );
	ARIADNE_TEST_EQUAL( theBinaryTreeCopy.is_leaf(), true );
	ARIADNE_TEST_EQUAL( theBinaryTreeCopy.is_enabled(), true );
	ARIADNE_PRINT_TEST_COMMENT( "The original tree" );
	ARIADNE_TEST_COMPARE( theBinaryTreeRoot.left_node(), !=, NULL );
	ARIADNE_TEST_COMPARE( theBinaryTreeRoot.right_node(), !=, NULL );
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.is_enabled(), false );
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.is_disabled(), false );
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.left_node()->is_leaf(), true );
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.right_node()->is_leaf(), true );
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.left_node()->is_enabled(), true );
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.right_node()->is_enabled(), true );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Disable the left tree leaf, recombine the tree and split the disabled leaf to level 3");
	ARIADNE_PRINT_TEST_COMMENT("A tree with disabled left node");
	theBinaryTreeRoot.left_node()->set_disabled();
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.left_node()->is_disabled(), true );

	ARIADNE_PRINT_TEST_COMMENT("A tree after recombination" );
	theBinaryTreeRoot.recombine();
	ARIADNE_TEST_COMPARE( theBinaryTreeRoot.left_node(), !=, NULL );
	ARIADNE_TEST_COMPARE( theBinaryTreeRoot.right_node(), !=, NULL );

	ARIADNE_PRINT_TEST_COMMENT( "A tree after splitting the disabled left node" );
	theBinaryTreeRoot.left_node()->mince(3);
	ARIADNE_TEST_EQUAL( theBinaryTreeRoot.left_node()->is_leaf(), true );
	
	ARIADNE_PRINT_TEST_COMMENT("A tree after splitting the root node" );
	theBinaryTreeRoot.mince(3);
	BooleanArray treeArray(9), leafArray(5);
	treeArray[0] = true; 
	treeArray[1] = false; leafArray[0] = false;
	treeArray[2] = true;
	treeArray[3] = true;
	treeArray[4] = false; leafArray[1] = true;
	treeArray[5] = false; leafArray[2] = true;
	treeArray[6] = true;
	treeArray[7] = false; leafArray[3] = true;
	treeArray[8] = false; leafArray[4] = true;
	BinaryTreeNode expected_binary_tree( treeArray , leafArray );
	ARIADNE_TEST_COMPARE( expected_binary_tree, ==, theBinaryTreeRoot );

	ARIADNE_PRINT_TEST_COMMENT("Testing the depth of the binary tree");
	ARIADNE_TEST_EQUAL(theBinaryTreeRoot.depth(), 3);

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Disable the left most enabled leaf, recombine the tree");
	theBinaryTreeRoot.right_node()->left_node()->left_node()->set_disabled();
	ARIADNE_PRINT_TEST_COMMENT("A tree with enabled left-most leaf");
	leafArray[1] = false;
	BinaryTreeNode expected_binary_tree1( treeArray , leafArray );
	ARIADNE_TEST_COMPARE( expected_binary_tree1, ==, theBinaryTreeRoot );
	
	ARIADNE_PRINT_TEST_COMMENT("A tree after recombination");
	theBinaryTreeRoot.recombine();
	expected_binary_tree1.right_node()->right_node()->make_leaf(true);
	ARIADNE_TEST_COMPARE( expected_binary_tree1, ==, theBinaryTreeRoot );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Testing the add_enabled( BinaryWord& ) method");
	ARIADNE_PRINT_TEST_COMMENT("Adding enabled node to an enabled node, the former one is defined by an empty path");

	BinaryTreeNode theNewBinaryTreeRoot(true);
	BinaryWord binaryWord;
	theNewBinaryTreeRoot.add_enabled( binaryWord );
	ARIADNE_PRINT_TEST_COMMENT("The binary tree should stay intact, i.e. consist of one enabled node");
	BinaryTreeNode expected_binary_tree2( true );
	ARIADNE_TEST_COMPARE( expected_binary_tree2, ==, theNewBinaryTreeRoot );
	
	ARIADNE_PRINT_TEST_COMMENT("Adding enabled sub node, path: (true,false) to a disabled node");
	theNewBinaryTreeRoot.set_disabled();
	binaryWord.push_back(true);
	binaryWord.push_back(false);
	theNewBinaryTreeRoot.add_enabled( binaryWord );
	ARIADNE_PRINT_TEST_COMMENT("The binary tree have five nodes, four leafs and one (true, false) is enabled");
	expected_binary_tree2.set_disabled();
	expected_binary_tree2.split();
	expected_binary_tree2.right_node()->split();
	expected_binary_tree2.right_node()->left_node()->set_enabled();
	ARIADNE_TEST_COMPARE( expected_binary_tree2, ==, theNewBinaryTreeRoot );

	ARIADNE_PRINT_TEST_COMMENT("Adding enabled sub node, path: (true) to a disabled node");
	binaryWord = BinaryWord();
	binaryWord.push_back(true);
	theNewBinaryTreeRoot.add_enabled( binaryWord );
	ARIADNE_PRINT_TEST_COMMENT("The binary tree have three nodes, two leafs and one (true) is enabled");
	expected_binary_tree2.right_node()->make_leaf(true);
	ARIADNE_TEST_COMPARE( expected_binary_tree2, ==, theNewBinaryTreeRoot );

	ARIADNE_PRINT_TEST_COMMENT("Adding enabled sub node, path: (false, true) to a disabled node");
	binaryWord = BinaryWord();
	binaryWord.push_back(false);
	binaryWord.push_back(true);
	theNewBinaryTreeRoot.add_enabled( binaryWord );
	ARIADNE_PRINT_TEST_COMMENT("The binary tree have five nodes, four leafs and two (false, true) and (true) are enabled");
	expected_binary_tree2.left_node()->split();
	expected_binary_tree2.left_node()->right_node()->set_enabled();
	ARIADNE_TEST_COMPARE( expected_binary_tree2, ==, theNewBinaryTreeRoot );
}

void test_grid_paving_cursor(){
	
	//Allocate the Grid
	Grid theGrid( make_vector<Float>("[-0.25, 0.25, 1.5]"), make_vector<Float>("[0.25, 0.25, 0.25]") );

	//Define the higth of the primary root cell.
	const uint theHeight = 2;

	//Create the binary tree;
	BinaryTreeNode * pRootTreeNode = new BinaryTreeNode(true);

	//Mince the tree to the certain higth
	pRootTreeNode->mince(3);

	//Mark some nodes as disabled
	pRootTreeNode->left_node()->right_node()->left_node()->set_disabled();
	pRootTreeNode->left_node()->left_node()->right_node()->set_disabled();
	pRootTreeNode->right_node()->right_node()->right_node()->set_disabled();

	//Create the path to the root node of the subpaving, which
	//corresponds to pRootTreeNode->left_node()->right_node()
	BinaryWord thePathToSubPavingRoot;
	thePathToSubPavingRoot.push_back(false);
	thePathToSubPavingRoot.push_back(true);
	
	//Create the GridTreeSubset
	GridTreeSubset theGridSubPavingSmall( theGrid, theHeight, thePathToSubPavingRoot, pRootTreeNode->left_node()->right_node() );
	//Create the Cursor for the GridTreeSubset
	GridTreeCursor theGSPCursorSmall( &theGridSubPavingSmall );

	//Create the GridTreeSubset
	GridTreeSubset theGridSubPavingLarge( theGrid, theHeight, BinaryWord(), pRootTreeNode );
	//Create the Cursor for the GridTreeSubset
	GridTreeCursor theGSPCursorLarge( &theGridSubPavingLarge );

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Check the created Cursor data on a small subpaving");
	Box expected_box = make_box("[-0.5,0]x[0.5,1]x[1.25,2.25]");
	ARIADNE_TEST_EQUAL( expected_box, theGSPCursorSmall.cell().box() );

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Check the created Cursor data on a large subpaving");
	expected_box = make_box("[-0.5,0.5]x[0,1]x[1.25,2.25]");
	ARIADNE_TEST_EQUAL( expected_box, theGSPCursorLarge.cell().box() );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Moving up from the Cursor's root node should cause the NotAllowedMoveException exception");
	ARIADNE_TEST_THROW( theGSPCursorSmall.move_up(), NotAllowedMoveException);
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Let's move the cursor of theGSPCursorLarge: left, right, we should be in the root of the smaller subpaving then");
	theGSPCursorLarge.move_left();
	theGSPCursorLarge.move_right();
	GridTreeSubset tmpSubPaving = *theGSPCursorLarge;
	GridTreeSubset expected_tree_subset( theGrid, theHeight, make_binary_word("01"), new BinaryTreeNode( *theGridSubPavingSmall.binary_tree()) );
	ARIADNE_PRINT_TEST_COMMENT("The cursor has moved and we have the following subpaving: ");
	ARIADNE_TEST_EQUAL( expected_tree_subset, tmpSubPaving);

	expected_box = make_box("[-0.5,0]x[0.5,1]x[1.25,2.25]");
	ARIADNE_PRINT_TEST_COMMENT("The box pointed by the cursor is: ");
	ARIADNE_TEST_EQUAL( expected_box, theGSPCursorLarge.cell().box() );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Check the is_enabled/disabled/leaf/root functions, while moving the cursor in the large and small subpavings");
	ARIADNE_TEST_EQUAL(theGSPCursorSmall.is_root(), true);
	ARIADNE_TEST_EQUAL(theGSPCursorSmall.is_leaf(), false);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_root(), false);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_leaf(), false);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_enabled(), false);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_disabled(), false);
	ARIADNE_PRINT_TEST_COMMENT("Moving to the left, to a disabled node");
	theGSPCursorLarge.move_left();
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_leaf(), true);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_root(), false);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_enabled(), false);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_disabled(), true);
	ARIADNE_PRINT_TEST_COMMENT("Moving up and to the right, to an enabled node");
	theGSPCursorLarge.move_up();
	theGSPCursorLarge.move_right();
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_leaf(), true);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_root(), false);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_enabled(), true);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_disabled(), false);
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test that moving to the left and right with the leaf nodes causes the NotAllowedMoveException exception");
	ARIADNE_PRINT_TEST_COMMENT("Get to the leaf node on the left");
	theGSPCursorSmall.move_left();

	ARIADNE_PRINT_TEST_COMMENT("Test the present state of theGSPCursorSmall");
	expected_box = make_box("[-0.5,0]x[0.5,1]x[1.25,1.75]");
	ARIADNE_PRINT_TEST_COMMENT("The box pointed by the cursor is: ");
	ARIADNE_TEST_EQUAL( expected_box, theGSPCursorSmall.cell().box() );

	ARIADNE_PRINT_TEST_COMMENT("Try to get to the node on the left");
	ARIADNE_TEST_THROW(theGSPCursorSmall.move_left(), NotAllowedMoveException );

	ARIADNE_PRINT_TEST_COMMENT("Try to get to the node on the right");
	ARIADNE_TEST_THROW(theGSPCursorSmall.move_right(), NotAllowedMoveException );

	ARIADNE_PRINT_TEST_COMMENT("The state of theGSPCursorSmall is supposed to remain unchanged");
	ARIADNE_PRINT_TEST_COMMENT("The box pointed by the cursor is: ");
	ARIADNE_TEST_EQUAL( expected_box, theGSPCursorSmall.cell().box() );

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test the set_enabled/set_disabled methods of the cursor");
	delete pRootTreeNode;
	pRootTreeNode = new BinaryTreeNode(true);
	ARIADNE_PRINT_TEST_COMMENT("Mince the enabled binary tree node to level 2");
	pRootTreeNode->mince(2);
	GridTreeSet theGridTreeSet( theGrid, theHeight, pRootTreeNode );
	//Create the Cursor for the GridTreeSubset
	GridTreeCursor theGPCursor( &theGridTreeSet );
	//Move to the leaf node
	ARIADNE_PRINT_TEST_COMMENT("Move the cursor to the right most leaf and check it");
	theGPCursor.move_right();
	theGPCursor.move_right();
	ARIADNE_TEST_EQUAL(theGPCursor.is_leaf(), true);
	ARIADNE_TEST_EQUAL(theGPCursor.is_root(), false);
	ARIADNE_TEST_EQUAL(theGPCursor.is_enabled(), true);
	ARIADNE_TEST_EQUAL(theGPCursor.is_disabled(), false);
	ARIADNE_PRINT_TEST_COMMENT("Set the leaf as disabled");
	theGPCursor.set_disabled();
	ARIADNE_TEST_EQUAL(theGPCursor.is_disabled(), true);
	ARIADNE_PRINT_TEST_COMMENT("Set the leaf (back) as enabled");
	theGPCursor.set_enabled();
	ARIADNE_TEST_EQUAL(theGPCursor.is_enabled(), true);
	ARIADNE_PRINT_TEST_COMMENT("Move one node up, to a non-leaf node and try to enable it, the NotALeafNodeException should be thrown");
	theGPCursor.move_up();
	ARIADNE_TEST_THROW( theGPCursor.set_disabled(), NotALeafNodeException);
}

void test_grid_paving_const_iterator(){
	std::vector< GridCell *> expected_result( 8 );
	
	//Allocate the Grid
	Grid theGrid( make_vector<Float>("[-0.25, 0.25]"), make_vector<Float>("[0.25, 0.25]") );
	
	//Define the higth of the primary root cell.
	const uint theHeight = 2;

	//Create the binary tree;
	BinaryTreeNode * pRootTreeNode = new BinaryTreeNode(true);

	//Create the GridTreeSubset
	GridTreeSubset theGridSubPavingLarge( theGrid, theHeight, BinaryWord(), pRootTreeNode );

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test the sequence in which GridPavingIterator goes through the tree leafs ");
	ARIADNE_PRINT_TEST_COMMENT("The tree depth is 0 and all leaf nodes are enabled");

	expected_result[0] = new GridCell( theGrid, 2, make_binary_word("") );
	
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result, theGridSubPavingLarge, 1 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result );

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test the sequence in which GridPavingIterator goes through the tree leafs ");
	ARIADNE_PRINT_TEST_COMMENT("The tree depth is 1 and all leaf nodes are enabled");

	//Mince the tree to the certain higth
	pRootTreeNode->mince(1);

	expected_result[0] =  new GridCell( theGrid, 2, make_binary_word("0") );
	expected_result[1] =  new GridCell( theGrid, 2, make_binary_word("1") );

	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result, theGridSubPavingLarge, 2 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result );

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test the sequence in which GridPavingIterator goes through the tree leafs ");
	ARIADNE_PRINT_TEST_COMMENT("The tree depth is 2 and all leaf nodes are enabled");

	//Mince the tree to the certain higth
	pRootTreeNode->mince(2);

	expected_result[0] =  new GridCell( theGrid, 2, make_binary_word("00") );
	expected_result[1] =  new GridCell( theGrid, 2, make_binary_word("01") );
	expected_result[2] =  new GridCell( theGrid, 2, make_binary_word("10") );
	expected_result[3] =  new GridCell( theGrid, 2, make_binary_word("11") );

	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result, theGridSubPavingLarge, 4 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result );

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test the sequence in which GridPavingIterator goes through the tree leafs ");
	ARIADNE_PRINT_TEST_COMMENT("The tree depth is 3 and all leaf nodes are enabled");

	//Mince the tree to the certain higth
	pRootTreeNode->mince(3);

	expected_result[0] = new GridCell( theGrid, 2, make_binary_word("000") );
	expected_result[1] = new GridCell( theGrid, 2, make_binary_word("001") );
	expected_result[2] = new GridCell( theGrid, 2, make_binary_word("010") );
	expected_result[3] = new GridCell( theGrid, 2, make_binary_word("011") );
	expected_result[4] = new GridCell( theGrid, 2, make_binary_word("100") );
	expected_result[5] = new GridCell( theGrid, 2, make_binary_word("101") );
	expected_result[6] = new GridCell( theGrid, 2, make_binary_word("110") );
	expected_result[7] = new GridCell( theGrid, 2, make_binary_word("111") );
	
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result, theGridSubPavingLarge, 8 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Disable some of the leaf nodes and test GridPavingIterator");
	ARIADNE_PRINT_TEST_COMMENT("The tree depth is 3 and some tree nodes are disabled");

	//Mark some nodes as disabled
	pRootTreeNode->left_node()->left_node()->left_node()->set_disabled();
	pRootTreeNode->left_node()->right_node()->right_node()->set_disabled();
	pRootTreeNode->right_node()->left_node()->left_node()->set_disabled();
	pRootTreeNode->right_node()->left_node()->right_node()->set_disabled();

	//Reuse the previous results
	std::vector< GridCell *>  expected_result_tmp(4);
	expected_result_tmp[0] = new GridCell( theGrid, 2, make_binary_word("001") );
	expected_result_tmp[1] = new GridCell( theGrid, 2, make_binary_word("010") );
	expected_result_tmp[2] = new GridCell( theGrid, 2, make_binary_word("110") );
	expected_result_tmp[3] = new GridCell( theGrid, 2, make_binary_word("111") );

	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_tmp, theGridSubPavingLarge, 4 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_tmp );

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Recombine the tree and test GridPavingIterator");
	ARIADNE_PRINT_TEST_COMMENT("The tree depth is 3 but some nodes are at depth 2");
	
	//Recombine the paving
	theGridSubPavingLarge.recombine();
	
	//Reuse some of the previous results
	expected_result_tmp[0] = new GridCell( theGrid, 2, make_binary_word("001") );
	expected_result_tmp[1] = new GridCell( theGrid, 2, make_binary_word("010") );
	expected_result_tmp[2] = new GridCell( theGrid, 2, make_binary_word("11") );
	
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_tmp, theGridSubPavingLarge, 3 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_tmp );

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Mince the tree back to level 3, enable/disable some nodes, recombine and and test GridPavingIterator");
	ARIADNE_PRINT_TEST_COMMENT("The tree depth is 3 but some nodes are at depth 2");

	//Mince the tree to the certain higth, (only for enabled nodes)
	pRootTreeNode->mince(3);
	//Split the disabled node
	pRootTreeNode->right_node()->left_node()->split();

	//Mark some nodes as disabled/enabled
	pRootTreeNode->right_node()->right_node()->left_node()->set_disabled();
	pRootTreeNode->right_node()->left_node()->left_node()->set_enabled();
	pRootTreeNode->right_node()->left_node()->right_node()->set_enabled();

	//Recombine the paving
	theGridSubPavingLarge.recombine();
	
	//Reuse some of the previous result strings
	expected_result_tmp[0] = new GridCell( theGrid, 2, make_binary_word("001") );
	expected_result_tmp[1] = new GridCell( theGrid, 2, make_binary_word("010") );
	expected_result_tmp[2] = new GridCell( theGrid, 2, make_binary_word("10") );
	expected_result_tmp[3] = new GridCell( theGrid, 2, make_binary_word("111") );

	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_tmp, theGridSubPavingLarge, 4 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_tmp );

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test how the constant Cursor can be retrieved from the Constant iterator");
	GridTreeSubset::const_iterator it = theGridSubPavingLarge.begin();
	const GridTreeCursor theGPCursor = it.cursor();
	ARIADNE_TEST_EQUAL(theGPCursor.is_leaf(), true);
	ARIADNE_TEST_EQUAL(theGPCursor.is_root(), false);
	ARIADNE_TEST_EQUAL(theGPCursor.is_enabled(), true);
	ARIADNE_TEST_EQUAL(theGPCursor.is_disabled(), false);
	ARIADNE_PRINT_TEST_COMMENT("Set the leaf as disabled");
	theGPCursor.set_disabled();
	ARIADNE_TEST_EQUAL(theGPCursor.is_disabled(), true);
	ARIADNE_PRINT_TEST_COMMENT("Set the leaf (back) as enabled");
	theGPCursor.set_enabled();
	ARIADNE_TEST_EQUAL(theGPCursor.is_enabled(), true);
}

void test_grid_sub_paving(){
	
	//Allocate the Grid, one Dimension
	Vector<Float> originOne(1); originOne.set(0, 0.0);
	Vector<Float> lengthsOne(1); lengthsOne.set(0, 0.5);
	const Grid theOneDimGrid( originOne, lengthsOne );

	//Define the higth of the primary root cell.
	const uint theHeight = 2;

	//Create the binary tree;
	BinaryTreeNode * pRootTreeNode = new BinaryTreeNode(true);

	//Create the path to the root node of the subpaving, which
	//corresponds to pRootTreeNode->left_node()->right_node()
	BinaryWord thePathToSubPavingRoot;
	thePathToSubPavingRoot.push_back(false);
	thePathToSubPavingRoot.push_back(true);

	//Create the GridTreeSubset
	GridTreeSubset theGridSPOneDim( theOneDimGrid, theHeight, thePathToSubPavingRoot, pRootTreeNode );
	
	//
	ARIADNE_PRINT_TEST_CASE_TITLE("Test Mincing operations of GridTreeSubset on the one dimensional Grid");
	ARIADNE_PRINT_TEST_COMMENT("The initial paving: ");
	ARIADNE_TEST_EQUAL( theGridSPOneDim.grid() , theOneDimGrid );
	ARIADNE_TEST_EQUAL( theGridSPOneDim.cell().height() , theHeight );
	ARIADNE_TEST_EQUAL( theGridSPOneDim.cell().word() , thePathToSubPavingRoot );
	ARIADNE_TEST_EQUAL( (*theGridSPOneDim.binary_tree()) , (*pRootTreeNode) );
	ARIADNE_TEST_EQUAL( theGridSPOneDim.cell().box() , make_box("[0,0.5]") );
	
	//Mince to the level two from the root of the paving, this should give
	//us 4 leaf nodes with intervals of width 0.125 (in the original space)
	theGridSPOneDim.mince(2);

	ARIADNE_PRINT_TEST_COMMENT("Minced the sub-paving to depth 2");
	std::vector< GridCell* > expected_result_arr(4);

	expected_result_arr[0] = new GridCell( theOneDimGrid, theHeight, make_binary_word("0100") );
	expected_result_arr[1] = new GridCell( theOneDimGrid, theHeight, make_binary_word("0101") );
	expected_result_arr[2] = new GridCell( theOneDimGrid, theHeight, make_binary_word("0110") );
	expected_result_arr[3] = new GridCell( theOneDimGrid, theHeight, make_binary_word("0111") );

	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theGridSPOneDim, 4 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );
	
	ARIADNE_PRINT_TEST_COMMENT("Recombine and subdivide the sub-paving to cell width 1.1");
	theGridSPOneDim.recombine();
	theGridSPOneDim.subdivide(1.1);
	ARIADNE_PRINT_TEST_COMMENT( "We should have a single node as in the initial sub paving: ");
	ARIADNE_TEST_EQUAL( theGridSPOneDim.binary_tree()->is_leaf(), true );
	ARIADNE_TEST_EQUAL( theGridSPOneDim.binary_tree()->is_enabled(), true );
	
	ARIADNE_PRINT_TEST_COMMENT("Subdivide the sub-paving to cell width 0.4, this should give us two sub cells");
	theGridSPOneDim.subdivide(0.4);

	expected_result_arr[0] = new GridCell( theOneDimGrid, theHeight, make_binary_word("010") );
	expected_result_arr[1] = new GridCell( theOneDimGrid, theHeight, make_binary_word("011") );

	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theGridSPOneDim, 2 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );
	
	ARIADNE_PRINT_TEST_COMMENT("Subdivide the sub-paving to cell width 0.126, this should give us four sub cells");
	theGridSPOneDim.subdivide(0.126);
	
	expected_result_arr[0] = new GridCell( theOneDimGrid, theHeight, make_binary_word("0100") );
	expected_result_arr[1] = new GridCell( theOneDimGrid, theHeight, make_binary_word("0101") );
	expected_result_arr[2] = new GridCell( theOneDimGrid, theHeight, make_binary_word("0110") );
	expected_result_arr[3] = new GridCell( theOneDimGrid, theHeight, make_binary_word("0111") );

	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theGridSPOneDim, 4 );
	//NOTE: To not clear the vector, the expected result is reused in the next test;

	ARIADNE_PRINT_TEST_COMMENT("Recombine and Subdivide the sub-paving to cell width 0.126, this should give us four sub cells");
	theGridSPOneDim.recombine();
	theGridSPOneDim.subdivide(0.126);

	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theGridSPOneDim, 4 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test Mincing operations of GridTreeSubset on the two dimensional Grid");
	//Allocate the Grid, one Dimension
	const Grid theTwoDimGrid( make_vector<Float>("[-0.25, 0.5]"), make_vector<Float>("[0.25, 0.5]") );

	//Create the GridTreeSubset
	GridTreeSubset theGridSPTwoDim( theTwoDimGrid, theHeight, thePathToSubPavingRoot, pRootTreeNode );

	ARIADNE_PRINT_TEST_COMMENT("Recombine and Mince the sub-paving to depth 2, this should give us four sub cells");
	theGridSPTwoDim.recombine();
	theGridSPTwoDim.mince(2);

	expected_result_arr[0] =  new GridCell( theTwoDimGrid, theHeight, make_binary_word("0100") );
	expected_result_arr[1] =  new GridCell( theTwoDimGrid, theHeight, make_binary_word("0101") );
	expected_result_arr[2] =  new GridCell( theTwoDimGrid, theHeight, make_binary_word("0110") );
	expected_result_arr[3] =  new GridCell( theTwoDimGrid, theHeight, make_binary_word("0111") );

	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theGridSPTwoDim, 4 );
	//Do not clean the array yet, the result well be reused

	ARIADNE_PRINT_TEST_COMMENT("Recombine and Subdivide the sub-paving to cell width 0.5, this should give us four sub cells");
	//At this moment the coordinate cell widths are: for x -- 0.25 and for y -- 0.5 
	theGridSPTwoDim.recombine();
	theGridSPTwoDim.subdivide(0.51);

	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theGridSPTwoDim, 4 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );

	ARIADNE_PRINT_TEST_COMMENT("Subdivide the sub-paving to cell width 0.4, this should give us sixteen sub cells");
	//At this moment the coordinate cell widths are: for x -- 0.25 and for y -- 0.5 
	theGridSPTwoDim.subdivide(0.4);

	BinaryTreeNode binaryTree( make_binary_word("1111001001100100111001001100100"), make_binary_word("1111111111111111") );
	GridTreeSubset expected_result( theTwoDimGrid, theHeight, make_binary_word("01"), &binaryTree );
	ARIADNE_PRINT_TEST_COMMENT("The sub paving with 16 leaf nodes: ");
	ARIADNE_TEST_EQUAL( expected_result, theGridSPTwoDim );
	
	ARIADNE_PRINT_TEST_COMMENT("The depth of the sub-paving's tree should be 4");
	ARIADNE_TEST_EQUAL( theGridSPTwoDim.depth(), 4 );
}

void test_grid_paving(){
	//Create a trivial 4-dimensional Grid
	Grid trivialFourDimGrid( make_vector<Float>("[0,0,0,0]"), make_vector<Float>("[1,1,1,1]") );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test allocation of a trivial GridTreeSet");
	GridTreeSet * pTrivialPaving = new GridTreeSet(4, true);
	GridTreeSet expected_result1( trivialFourDimGrid , 0, make_binary_word("0"), make_binary_word("1") );
	ARIADNE_PRINT_TEST_COMMENT("A trivial paving for [0,1]x[0,1]x[0,1]x[0,1], enabled cell: ");
	ARIADNE_TEST_EQUAL( expected_result1, (*pTrivialPaving) );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test GridTreeSet copy constructor");
	GridTreeSet theTrivialPaving( ( *pTrivialPaving ) );
	pTrivialPaving->mince(1);
	
	ARIADNE_PRINT_TEST_COMMENT("Minced trivial paving for [0,1]x[0,1]x[0,1]x[0,1], enabled cell: ");
	GridTreeSet expected_result_minced(trivialFourDimGrid , 0, make_binary_word("100"), make_binary_word("11") );
	ARIADNE_TEST_EQUAL( expected_result_minced, (*pTrivialPaving) );

	ARIADNE_PRINT_TEST_COMMENT("A copy of the original paving, should stay unchanged: ");
	ARIADNE_TEST_EQUAL( expected_result1, theTrivialPaving );

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test GridTreeSet (Grid, Box) constructor");
	//Allocate the Grid, one Dimension
	const Grid theTwoDimGrid( make_vector<Float>("[-0.25,0.5]"), make_vector<Float>("[0.25,0.5]") );
	//Note: the box is related to the grid, but not to the original space
	GridTreeSet theTwoDimPaving( theTwoDimGrid, make_box("[0,1.5]x[-1.5,3.5]") );
	ARIADNE_PRINT_TEST_COMMENT("The resulting GridTreeSet: ");
	GridTreeSet expected_result2( theTwoDimGrid, 4, make_binary_word("0"), make_binary_word("0") );
	ARIADNE_TEST_EQUAL( expected_result2, theTwoDimPaving );

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test GridTreeSet (Grid, Height, BooleanArray, BooleanArray ) constructor");
	BooleanArray theTree(9), theEnabledLeafs(5);
	theTree[0] = true;
	theTree[1] = true;
	theTree[2] = false;  theEnabledLeafs[0] = true;
	theTree[3] = false;  theEnabledLeafs[1] = false;
	theTree[4] = true;
	theTree[5] = true;
	theTree[6] = false;  theEnabledLeafs[2] = true;
	theTree[7] = false;  theEnabledLeafs[3] = false;
	theTree[8] = false;  theEnabledLeafs[4] = true;
	GridTreeSet theTwoDimTreePaving( theTwoDimGrid, 2, theTree, theEnabledLeafs );
	ARIADNE_PRINT_TEST_COMMENT("The resulting GridTreeSet: ");
	BinaryTreeNode * pRootTreeNode = new BinaryTreeNode(false);
	pRootTreeNode->split();
	pRootTreeNode->left_node()->split();
	pRootTreeNode->left_node()->left_node()->set_enabled();
	pRootTreeNode->right_node()->split();
	pRootTreeNode->right_node()->right_node()->set_enabled();
	pRootTreeNode->right_node()->left_node()->split();
	pRootTreeNode->right_node()->left_node()->left_node()->set_enabled();
	GridTreeSet expected_result3( theTwoDimGrid, 2, pRootTreeNode );
	ARIADNE_TEST_EQUAL( expected_result3, theTwoDimTreePaving );

	// !!!
	//Create a trivial 2-dimensional Grid
	Grid trivialTwoDimGrid( make_vector<Float>("[0,0]"), make_vector<Float>("[1,1]") );
	ARIADNE_PRINT_TEST_CASE_TITLE("Test GridTreeSet::restrict_to_height( const uint theHeight )");
	GridTreeSet initialTreeSetOne( trivialTwoDimGrid, 2, make_binary_word("1111001000100"), make_binary_word("10010001") );
	GridTreeSet initialTreeSetOneCopyOne( initialTreeSetOne );
	GridTreeSet initialTreeSetOneCopyTwo( initialTreeSetOne );
	
	ARIADNE_PRINT_TEST_COMMENT("The initialTreeSetOne.cell().height() == 2, theHeight == 3, nothing should change");
	initialTreeSetOne.restrict_to_height( 3 );
	ARIADNE_TEST_EQUAL( initialTreeSetOneCopyOne, initialTreeSetOne );
	
	ARIADNE_PRINT_TEST_COMMENT("The initialTreeSetOneCopyOne.cell().height() == 2, theHeight == 1");
	initialTreeSetOneCopyOne.restrict_to_height( 1 );
	GridTreeSet expectedTreeSetOne( trivialTwoDimGrid, 2, make_binary_word("11110010000"), make_binary_word("100100") );
	ARIADNE_TEST_EQUAL( expectedTreeSetOne, initialTreeSetOneCopyOne );
	
	ARIADNE_PRINT_TEST_COMMENT("The initialTreeSetOneCopyTwo.cell().height() == 2, theHeight == 0");
	initialTreeSetOneCopyTwo.restrict_to_height( 0 );
	GridTreeSet expectedTreeSetTwo( trivialTwoDimGrid, 2, make_binary_word("111010000"), make_binary_word("00100") );
	ARIADNE_TEST_EQUAL( expectedTreeSetTwo, initialTreeSetOneCopyTwo );
	
	ARIADNE_PRINT_TEST_COMMENT("The initialTreeSetTwo.cell().height() == 1, theHeight == 0");
	GridTreeSet initialTreeSetTwo( trivialTwoDimGrid, 1, make_binary_word("11000"), make_binary_word("100") );
	initialTreeSetTwo.restrict_to_height( 0 );
	GridTreeSet expectedTreeSetThree( trivialTwoDimGrid, 1, make_binary_word("100"), make_binary_word("00") );
	ARIADNE_TEST_EQUAL( expectedTreeSetThree, initialTreeSetTwo );
}

void test_grid_paving_cell(){
	BinaryWord expected_result;

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test the static methods of the GridCell");
	BinaryWord theBinaryPath;
	
	theBinaryPath = GridCell::primary_cell_path( 1, 0, 0 );
	ARIADNE_PRINT_TEST_COMMENT( "Dimension: 1, topCellHeight: 0, bottomCellHeight: 0"  );
	ARIADNE_TEST_EQUAL( expected_result , theBinaryPath );
	
	expected_result = make_binary_word( "1" ) ;
	theBinaryPath = GridCell::primary_cell_path( 1, 1, 0 );
	ARIADNE_PRINT_TEST_COMMENT( "Dimension: 1, topCellHeight: 1, bottomCellHeight: 0" );
	ARIADNE_TEST_EQUAL( expected_result , theBinaryPath );
	
	expected_result = make_binary_word( "01" );
	theBinaryPath = GridCell::primary_cell_path( 1, 2, 0 );
	ARIADNE_PRINT_TEST_COMMENT( "Dimension: 1, topCellHeight: 2, bottomCellHeight: 0" );
	ARIADNE_TEST_EQUAL( expected_result , theBinaryPath );
	
	expected_result = make_binary_word( "0" );
	theBinaryPath = GridCell::primary_cell_path( 1, 2, 1 );
	ARIADNE_PRINT_TEST_COMMENT( "Dimension: 1, topCellHeight: 2, bottomCellHeight: 1" );
	ARIADNE_TEST_EQUAL( expected_result , theBinaryPath );
	
	expected_result.clear();
	theBinaryPath = GridCell::primary_cell_path( 1, 2, 2 );
	ARIADNE_PRINT_TEST_COMMENT( "Dimension: 1, topCellHeight: 2, bottomCellHeight: 2" );
	ARIADNE_TEST_EQUAL( expected_result , theBinaryPath );
}

void test_adjoin_operation_one(){
	//Allocate a trivial Grid
	Grid theGrid(2, 1.0);
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test adjoining a GridCell to the GridTreeSet");
	//Define the GridCell that is rooted to a high primary cell
	const int theHigherCellHeight = 2;
	BinaryWord theHigherCellPath;
	theHigherCellPath.push_back(true);
	theHigherCellPath.push_back(false);
	theHigherCellPath.push_back(true);
	theHigherCellPath.push_back(false);
	GridCell theHigherLevelCell( theGrid, theHigherCellHeight, theHigherCellPath );

	ARIADNE_PRINT_TEST_COMMENT("The GridCell with the primary root cell height = 2");
	Box expected_box = make_box("[2,3]x[-1,0]");
	ARIADNE_PRINT_TEST_COMMENT("The initial GridCell, as given by it's box: ");
	ARIADNE_TEST_EQUAL( expected_box, theHigherLevelCell.box() );
	
	//Define the higth of the primary root cell.
	//Create the GridTreeSet with the box related to the grid, but not to the original space
	GridTreeSet theOneCellPaving( theGrid, true );
	ARIADNE_PRINT_TEST_COMMENT("The GridTreeSet with the primary root cell height = 0");
	GridTreeSet expected_one_cell_paving( theGrid, 0, make_binary_word("0"), make_binary_word("1") );
	ARIADNE_TEST_EQUAL( expected_one_cell_paving, theOneCellPaving );
	
	ARIADNE_PRINT_TEST_COMMENT("The GridTreeSet after adding the cell: ");
	theOneCellPaving.adjoin( theHigherLevelCell );
	GridTreeSet expected_two_cell_paving( theGrid, 2, make_binary_word("111010001101000"), make_binary_word("00100100") );
	ARIADNE_TEST_EQUAL( expected_two_cell_paving, theOneCellPaving );
}

void test_adjoin_operation_two(){
	
	//Allocate a trivial Grid
	Grid theGrid(2, 1.0);
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test adjoining a GridCell to the GridTreeSet");
	//Define the GridCell that is rooted to the lower primary cell
	const int theLowerCellHeight = 1;
	BinaryWord theLowerCellPath;
	theLowerCellPath.push_back(true);
	theLowerCellPath.push_back(true);
	GridCell theLowerLevelCell( theGrid, theLowerCellHeight, theLowerCellPath );

	ARIADNE_PRINT_TEST_COMMENT("The GridCell with the primary root cell height = 1");
	Box expected_box = make_box("[0,1]x[0,1]");
	ARIADNE_PRINT_TEST_COMMENT("The box of the initial GridCell: ");
	ARIADNE_TEST_EQUAL( expected_box, theLowerLevelCell.box() );
	
	//Define the higth of the primary root cell.
	const uint theHeight = 2;
	//Create the binary tree;
	BinaryTreeNode * pRootTreeNode = new BinaryTreeNode(false);
	pRootTreeNode->split();
	pRootTreeNode->right_node()->split();
	pRootTreeNode->right_node()->left_node()->split();
	pRootTreeNode->right_node()->left_node()->right_node()->split();
	pRootTreeNode->right_node()->left_node()->right_node()->left_node()->set_enabled();
	
	//Create the GridTreeSet with the box is related to the grid, but not to the original space
	GridTreeSet theOneCellPaving( theGrid, theHeight, pRootTreeNode );
	ARIADNE_PRINT_TEST_COMMENT("The GridTreeSet with the primary root cell height = 2");
	GridTreeSet expected_tree_set( theGrid, theHeight, make_binary_word("101101000"), make_binary_word("00100") );
	ARIADNE_TEST_EQUAL( expected_tree_set, theOneCellPaving );
	
	ARIADNE_PRINT_TEST_COMMENT("The GridTreeSet after adding the cell: ");
	theOneCellPaving.adjoin( theLowerLevelCell );
	GridTreeSet expected_tree_set_result( theGrid, theHeight, make_binary_word("111010001101000"), make_binary_word("00100100") );
	ARIADNE_TEST_EQUAL( expected_tree_set_result, theOneCellPaving );
}

void test_adjoin_operation_three(){
	string expected_result;
	
	//Allocate a trivial Grid
	Grid theGrid(2, 1.0);
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test adjoining a GridCell to the GridTreeSet");
	//Define the GridCell that is rooted to the same primary cell
	const int theLowerCellHeight = 2;
	BinaryWord theLowerCellPath;
	theLowerCellPath.push_back(false);
	theLowerCellPath.push_back(false);
	theLowerCellPath.push_back(true);
	theLowerCellPath.push_back(true);
	GridCell theLowerLevelCell( theGrid, theLowerCellHeight, theLowerCellPath );
	
	ARIADNE_PRINT_TEST_COMMENT("The GridCell with the primary root cell height = 2");
	Box expected_box = make_box("[0,1]x[0,1]");
	ARIADNE_PRINT_TEST_COMMENT("The initial GridCell: ");
	ARIADNE_TEST_EQUAL( expected_box, theLowerLevelCell.box() );
		
	//Define the higth of the primary root cell.
	const uint theHeight = 2;
	//Create the binary tree;
	BinaryTreeNode * pRootTreeNode = new BinaryTreeNode(false);
	pRootTreeNode->split();
	pRootTreeNode->right_node()->split();
	pRootTreeNode->right_node()->left_node()->split();
	pRootTreeNode->right_node()->left_node()->right_node()->split();
	pRootTreeNode->right_node()->left_node()->right_node()->left_node()->set_enabled();
	
	//Create the GridTreeSet with the box is related to the grid, but not to the original space
	ARIADNE_PRINT_TEST_COMMENT("The GridTreeSet with the primary root cell height = 2");
	GridTreeSet theOneCellPaving( theGrid, theHeight, pRootTreeNode );
	GridTreeSet expected_tree_set( theGrid, theHeight, make_binary_word("101101000"), make_binary_word("00100") );
	ARIADNE_TEST_EQUAL( expected_tree_set, theOneCellPaving );
	
	ARIADNE_PRINT_TEST_COMMENT("The GridTreeSet after adding the cell: ");
	theOneCellPaving.adjoin( theLowerLevelCell );
	GridTreeSet expected_tree_set_result( theGrid, theHeight, make_binary_word("111010001101000"), make_binary_word("00100100") );
	ARIADNE_TEST_EQUAL( expected_tree_set_result, theOneCellPaving );
}

void test_adjoin_outer_approximation_operation(){
	//Allocate a trivial Grid
	Grid theTrivialGrid(2, 1.0);
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test adjoining_outer_approximation a SetInterface to the GridTreeSet");
        ImageSet initialRectangle( make_box("[-0.5,1.5]x[-0.3,1.0]") );
	
	//Define the higth of the primary root cell.
	const uint theHeight = 2;
	//Create the binary tree;
	BinaryTreeNode * pRootTreeNode = new BinaryTreeNode(false);
	pRootTreeNode->split();
	pRootTreeNode->right_node()->split();
	pRootTreeNode->right_node()->left_node()->split();
	pRootTreeNode->right_node()->left_node()->right_node()->split();
	pRootTreeNode->right_node()->left_node()->right_node()->left_node()->set_enabled();
	
	//Create the GridTreeSet with the box related to the grid, but not to the original space
	ARIADNE_PRINT_TEST_COMMENT("The initial GridTreeSet with the primary root cell height = 2");
	GridTreeSet theOneCellPaving( theTrivialGrid, theHeight, pRootTreeNode );
	BinaryWord tree = make_binary_word("101101000");
	BinaryWord leaves = make_binary_word("00100");
	GridTreeSet expected_grid_tree_set1( theTrivialGrid, theHeight, tree, leaves );
	ARIADNE_TEST_EQUAL( expected_grid_tree_set1, theOneCellPaving );
	
	ARIADNE_PRINT_TEST_COMMENT("The GridTreeSet after adding the cell: ");
	theOneCellPaving.adjoin_outer_approximation( static_cast<const LocatedSetInterface&>(initialRectangle), 2 );
	tree = make_binary_word("1110101111110010011001001110010011001001111001000111001000111110010011001001001111001000000");
	leaves = make_binary_word("0001011111010111111010010100010111111010100000");
	GridTreeSet expected_grid_tree_set2( theTrivialGrid, 4, tree, leaves );
	ARIADNE_TEST_EQUAL( expected_grid_tree_set2, theOneCellPaving );

	ARIADNE_PRINT_TEST_COMMENT("Recombined GridTreeSet after adding the cell: ");
	std::vector< GridCell* > expected_result_arr(16);
	expected_result_arr[0] = new GridCell( theTrivialGrid, 4, make_binary_word("[0,0,1,1,0,0,0,0,0,1]") );
	expected_result_arr[1] = new GridCell( theTrivialGrid, 4, make_binary_word("[0,0,1,1,0,0,0,0,1,1]") );
	expected_result_arr[2] = new GridCell( theTrivialGrid, 4, make_binary_word("[0,0,1,1,0,0,0,1]") );
	expected_result_arr[3] = new GridCell( theTrivialGrid, 4, make_binary_word("[0,0,1,1,0,0,1,0,0,1]") );
	expected_result_arr[4] = new GridCell( theTrivialGrid, 4, make_binary_word("[0,0,1,1,0,0,1,0,1,1]") );
	expected_result_arr[5] = new GridCell( theTrivialGrid, 4, make_binary_word("[0,0,1,1,0,0,1,1]") );
	expected_result_arr[6] = new GridCell( theTrivialGrid, 4, make_binary_word("[0,0,1,1,0,1,0,0,0,0]") );
	expected_result_arr[7] = new GridCell( theTrivialGrid, 4, make_binary_word("[0,0,1,1,0,1,0,0,1,0]") );
	expected_result_arr[8] = new GridCell( theTrivialGrid, 4, make_binary_word("[0,0,1,1,0,1,1,0,0,0]") );
	expected_result_arr[9] = new GridCell( theTrivialGrid, 4, make_binary_word("[0,0,1,1,0,1,1,0,1,0]") );
	expected_result_arr[10] = new GridCell( theTrivialGrid, 4, make_binary_word("[0,0,1,1,1,0,0,0,0,1]") );
	expected_result_arr[11] = new GridCell( theTrivialGrid, 4, make_binary_word("[0,0,1,1,1,0,0,0,1,1]") );
	expected_result_arr[12] = new GridCell( theTrivialGrid, 4, make_binary_word("[0,0,1,1,1,0,0,1]") );
	expected_result_arr[13] = new GridCell( theTrivialGrid, 4, make_binary_word("[0,0,1,1,1,0,1,0]") );
	expected_result_arr[14] = new GridCell( theTrivialGrid, 4, make_binary_word("[0,0,1,1,1,1,0,0,0,0]") );
	expected_result_arr[15] = new GridCell( theTrivialGrid, 4, make_binary_word("[0,0,1,1,1,1,0,0,1,0]") );
	theOneCellPaving.recombine();
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theOneCellPaving, 16 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Create an outer_approximation of the rectangle on the scaling grid and get the GridTreeSet");
	Grid theScalingGrid(2, 2.0);
	GridTreeSet theOuterApproxGridTreeSet = outer_approximation( static_cast<LocatedSetInterface&>(initialRectangle), theScalingGrid, 2 );
	//IVAN S. ZAPREEV
	//NOTE: The recombination is needed because in the scaling Grid doing
	//	outer_approximation( theScalingGrid, initialRectangle, 2 )
	//will result in subdivisions that are equal to unit cell in the original
	//space in this case we will get much more elements, e.g. the cells
	//  [-1,0]x[0,2], [0,2]x[0,2] will be subdivided as well
	theOuterApproxGridTreeSet.recombine();
	expected_result_arr[0] = new GridCell( theScalingGrid, 3, make_binary_word("[1,1,0,0,0,0,1,1]") );
	expected_result_arr[1] = new GridCell( theScalingGrid, 3, make_binary_word("[1,1,0,0,0,1,1]") );
	expected_result_arr[2] = new GridCell( theScalingGrid, 3, make_binary_word("[1,1,0,0,1,0,0,1]") );
	expected_result_arr[3] = new GridCell( theScalingGrid, 3, make_binary_word("[1,1,0,0,1,0,1,1]") );
	expected_result_arr[4] = new GridCell( theScalingGrid, 3, make_binary_word("[1,1,0,0,1,1]") );
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theOuterApproxGridTreeSet, 5 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );
}

void test_restrict() {
	std::vector< GridCell* > expected_result_arr(3);
	
	//Allocate a trivial Grid
	Grid theTrivialGrid(2, 1.0);
	
	//Define the higths of the primary root cell.
	const uint theHeightTwo = 2;
	const uint theHeightThree = 3;

	//Create the binary tree with two enabled nodes;
	BinaryTreeNode * pTwoEnabledNodeTreeH2 = new BinaryTreeNode(false);
	pTwoEnabledNodeTreeH2->split();
	pTwoEnabledNodeTreeH2->right_node()->split();
	pTwoEnabledNodeTreeH2->right_node()->left_node()->split();
	pTwoEnabledNodeTreeH2->right_node()->left_node()->left_node()->set_enabled();
	pTwoEnabledNodeTreeH2->right_node()->left_node()->right_node()->split();
	pTwoEnabledNodeTreeH2->right_node()->left_node()->right_node()->left_node()->set_enabled();
	//Create the GridTreeSet
	GridTreeSet theTwoCellPavingH2( theTrivialGrid, theHeightTwo, pTwoEnabledNodeTreeH2 );

	//Create another binary tree
	BinaryTreeNode * pThreeEnabledNodeTreeH2 = new BinaryTreeNode( *pTwoEnabledNodeTreeH2 );
	pThreeEnabledNodeTreeH2->left_node()->split();
	pThreeEnabledNodeTreeH2->left_node()->left_node()->set_enabled();
	pThreeEnabledNodeTreeH2->right_node()->left_node()->right_node()->right_node()->set_enabled();
	pThreeEnabledNodeTreeH2->right_node()->left_node()->right_node()->left_node()->set_disabled();
	pThreeEnabledNodeTreeH2->right_node()->left_node()->left_node()->split();
	pThreeEnabledNodeTreeH2->right_node()->left_node()->left_node()->right_node()->set_disabled();
	
	//Create another GridTreeSet
	GridTreeSet theThreeCellPavingH2( theTrivialGrid, theHeightTwo, pThreeEnabledNodeTreeH2 );

	//Create another binary tree
	BinaryTreeNode * pThreeEnabledNodeTreeH3 = new BinaryTreeNode(false);
	pThreeEnabledNodeTreeH3->split();
	pThreeEnabledNodeTreeH3->right_node()->split();
	pThreeEnabledNodeTreeH3->right_node()->right_node()->copy_from( pThreeEnabledNodeTreeH2 );
	BinaryWord thePathToSubPavingRoot;
	thePathToSubPavingRoot.push_back(true);
	thePathToSubPavingRoot.push_back(true);
	//Create the GridTreeSubset
	GridTreeSubset theThreeCellSubPavingH3( theTrivialGrid, theHeightThree , thePathToSubPavingRoot, pThreeEnabledNodeTreeH3->right_node()->right_node() );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test restrict operation: GridTreeSet1.restrict( GridTreeSet2 )");
	ARIADNE_PRINT_TEST_COMMENT("The initial GridTreeSet1: ");
	expected_result_arr[0] = new GridCell( theTrivialGrid, 2, make_binary_word("[0,0]") );
	expected_result_arr[1] = new GridCell( theTrivialGrid, 2, make_binary_word("[1,0,0,0]") );
	expected_result_arr[2] = new GridCell( theTrivialGrid, 2, make_binary_word("[1,0,1,1]") );
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theThreeCellPavingH2, 3 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );

	ARIADNE_PRINT_TEST_COMMENT("The initial GridTreeSet2: ");
	expected_result_arr[0] = new GridCell( theTrivialGrid, 2, make_binary_word("[1,0,0]") );
	expected_result_arr[1] = new GridCell( theTrivialGrid, 2, make_binary_word("[1,0,1,0]") );
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theTwoCellPavingH2, 2 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );

	ARIADNE_PRINT_TEST_COMMENT("The result after restrict: ");
	theThreeCellPavingH2.restrict( theTwoCellPavingH2 );
	expected_result_arr[0] = new GridCell( theTrivialGrid, 2, make_binary_word("[1,0,0,0]") );
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theThreeCellPavingH2, 1 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test restrict operation: GridTreeSet.restrict( GridTreeSubset )");
	ARIADNE_PRINT_TEST_COMMENT("The initial GridTreeSet: ");
	expected_result_arr[0] = new GridCell( theTrivialGrid, 2, make_binary_word("[1,0,0]") );
	expected_result_arr[1] = new GridCell( theTrivialGrid, 2, make_binary_word("[1,0,1,0]") );
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theTwoCellPavingH2, 2 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );

	ARIADNE_PRINT_TEST_COMMENT("The initial GridTreeSubset: ");
	expected_result_arr[0] = new GridCell( theTrivialGrid, 3, make_binary_word("[1,1,0,0]") );
	expected_result_arr[1] = new GridCell( theTrivialGrid, 3, make_binary_word("[1,1,1,0,0,0]") );
	expected_result_arr[2] = new GridCell( theTrivialGrid, 3, make_binary_word("[1,1,1,0,1,1]") );
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theThreeCellSubPavingH3, 3 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );

	ARIADNE_PRINT_TEST_COMMENT("The result after restrict: ");
	theTwoCellPavingH2.restrict( theThreeCellSubPavingH3 );
	expected_result_arr[0] = new GridCell( theTrivialGrid, 3, make_binary_word("[1,1,1,0,0,0]") );
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theTwoCellPavingH2, 1 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );
	
	//TODO: Test the case when the GridTreeSet has primary cell of the level 3
	//	The GridTreeSubset is at level 1 and it's primary cell is at level 2
}

void test_remove_one() {
	std::vector< GridCell* > expected_result_arr(3);
	
	//Allocate a trivial Grid
	Grid theTrivialGrid(2, 1.0);
	
	//Define the higths of the primary root cell.
	const uint theHeightTwo = 2;
	const uint theHeightThree = 3;

	//Create the binary tree with two enabled nodes;
	BinaryTreeNode * pTwoEnabledNodeTreeH2 = new BinaryTreeNode(false);
	pTwoEnabledNodeTreeH2->split();
	pTwoEnabledNodeTreeH2->right_node()->split();
	pTwoEnabledNodeTreeH2->right_node()->left_node()->split();
	pTwoEnabledNodeTreeH2->right_node()->left_node()->left_node()->set_enabled();
	pTwoEnabledNodeTreeH2->right_node()->left_node()->right_node()->split();
	pTwoEnabledNodeTreeH2->right_node()->left_node()->right_node()->left_node()->set_enabled();
	//Create the GridTreeSet
	GridTreeSet theTwoCellPavingH2( theTrivialGrid, theHeightTwo, pTwoEnabledNodeTreeH2 );

	//Create another binary tree
	BinaryTreeNode * pThreeEnabledNodeTreeH2 = new BinaryTreeNode( *pTwoEnabledNodeTreeH2 );
	pThreeEnabledNodeTreeH2->left_node()->split();
	pThreeEnabledNodeTreeH2->left_node()->left_node()->set_enabled();
	pThreeEnabledNodeTreeH2->right_node()->left_node()->right_node()->right_node()->set_enabled();
	pThreeEnabledNodeTreeH2->right_node()->left_node()->right_node()->left_node()->set_disabled();
	pThreeEnabledNodeTreeH2->right_node()->left_node()->left_node()->split();
	pThreeEnabledNodeTreeH2->right_node()->left_node()->left_node()->right_node()->set_disabled();
	
	//Create another GridTreeSet
	GridTreeSet theThreeCellPavingH2( theTrivialGrid, theHeightTwo, pThreeEnabledNodeTreeH2 );

	//Create another binary tree
	BinaryTreeNode * pThreeEnabledNodeTreeH3 = new BinaryTreeNode(false);
	pThreeEnabledNodeTreeH3->split();
	pThreeEnabledNodeTreeH3->right_node()->split();
	pThreeEnabledNodeTreeH3->right_node()->right_node()->copy_from( pThreeEnabledNodeTreeH2 );
	BinaryWord thePathToSubPavingRoot;
	thePathToSubPavingRoot.push_back(true);
	thePathToSubPavingRoot.push_back(true);
	//Create the GridTreeSubset
	GridTreeSubset theThreeCellSubPavingH3( theTrivialGrid, theHeightThree , thePathToSubPavingRoot, pThreeEnabledNodeTreeH3->right_node()->right_node() );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test remove operation: GridTreeSet1.remove( GridTreeSet2 )");
	ARIADNE_PRINT_TEST_COMMENT("The initial GridTreeSet1: ");
	expected_result_arr[0] = new GridCell( theTrivialGrid, 2, make_binary_word("[0,0]") );
	expected_result_arr[1] = new GridCell( theTrivialGrid, 2, make_binary_word("[1,0,0,0]") );
	expected_result_arr[2] = new GridCell( theTrivialGrid, 2, make_binary_word("[1,0,1,1]") );
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theThreeCellPavingH2, 3 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );

	ARIADNE_PRINT_TEST_COMMENT("The initial GridTreeSet2: ");
	expected_result_arr[0] = new GridCell( theTrivialGrid, 2, make_binary_word("[1,0,0]") );
	expected_result_arr[1] = new GridCell( theTrivialGrid, 2, make_binary_word("[1,0,1,0]") );
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theTwoCellPavingH2, 2 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );

	ARIADNE_PRINT_TEST_COMMENT("The result after removal: ");
	theThreeCellPavingH2.remove( theTwoCellPavingH2 );
	expected_result_arr[0] = new GridCell( theTrivialGrid, 2, make_binary_word("[0,0]") );
	expected_result_arr[1] = new GridCell( theTrivialGrid, 2, make_binary_word("[1,0,1,1]") );
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theThreeCellPavingH2, 2 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test remove operation: GridTreeSet.remove( GridTreeSubset )");
	ARIADNE_PRINT_TEST_COMMENT("The initial GridTreeSet: ");
	expected_result_arr[0] = new GridCell( theTrivialGrid, 2, make_binary_word("[1,0,0]") );
	expected_result_arr[1] = new GridCell( theTrivialGrid, 2, make_binary_word("[1,0,1,0]") );
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theTwoCellPavingH2, 2 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );

	ARIADNE_PRINT_TEST_COMMENT("The initial GridTreeSubset: ");
	expected_result_arr[0] = new GridCell( theTrivialGrid, 3, make_binary_word("[1,1,0,0]") );
	expected_result_arr[1] = new GridCell( theTrivialGrid, 3, make_binary_word("[1,1,1,0,0,0]") );
	expected_result_arr[2] = new GridCell( theTrivialGrid, 3, make_binary_word("[1,1,1,0,1,1]") );
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theThreeCellSubPavingH3, 3 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );

	ARIADNE_PRINT_TEST_COMMENT("The result after remove: ");
	theTwoCellPavingH2.remove( theThreeCellSubPavingH3 );
	expected_result_arr[0] = new GridCell( theTrivialGrid, 3, make_binary_word("[1,1,1,0,0,1]") );
	expected_result_arr[1] = new GridCell( theTrivialGrid, 3, make_binary_word("[1,1,1,0,1,0]") );
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theTwoCellPavingH2, 2 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );
	
	//TODO: Test the case when the GridTreeSet has primary cell of the level 3
	//	The GridTreeSubset is at level 1 and it's primary cell is at level 2
}

void test_remove_two() {
	std::vector< GridCell* > expected_result_arr(4);
	
	//Allocate a trivial Grid
	Grid theTrivialGrid(2, 1.0);
	
	//Define the higths of the primary root cell.
	const uint theHeightZero = 0;
	const uint theHeightOne = 1;
	const uint theHeightTwo = 2;
	
	//Create a GridTreeSet and its copies
	GridTreeSet theSet01( theTrivialGrid, theHeightOne, make_binary_word("110010100"), make_binary_word("10001") );
	GridTreeSet theSet02( theSet01 );
	GridTreeSet theSet03( theSet01 );
	GridTreeSet theSet04( theSet01 );
	GridTreeSet theSet05( theSet01 );
	GridTreeSet theSet06( theSet01 );
	
	//Create the cells we will be removing and test removals on the copies of the set
	
	ARIADNE_PRINT_TEST_CASE_TITLE("Remove a GridCell (p.c. height=0) form a GridTreeSet(p.c. height=1): The cell is a subset.");
	GridCell lowerPrimaryCellSubset( theTrivialGrid, theHeightZero, make_binary_word("11") );
	theSet01.remove( lowerPrimaryCellSubset );
	expected_result_arr[0] = new GridCell( theTrivialGrid, theHeightOne, make_binary_word("00") );
	expected_result_arr[1] = new GridCell( theTrivialGrid, theHeightOne, make_binary_word("1110") );
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theSet01, 2 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );
	
	ARIADNE_PRINT_TEST_CASE_TITLE("Remove a GridCell (p.c. height=0) form a GridTreeSet(p.c. height=1): The cell does not intersect the set.");
	GridCell lowerPrimaryCellNoIntersection( theTrivialGrid, theHeightZero, make_binary_word("00") );
	theSet02.remove( lowerPrimaryCellNoIntersection );
	expected_result_arr[0] = new GridCell( theTrivialGrid, theHeightOne, make_binary_word("00") );
	expected_result_arr[1] = new GridCell( theTrivialGrid, theHeightOne, make_binary_word("111") );
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theSet02, 2 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );
	
	ARIADNE_PRINT_TEST_CASE_TITLE("Remove a GridCell (p.c. height=0) form a GridTreeSet(p.c. height=1): The cell intersects the set.");
	GridCell lowerPrimaryCellIntersection( theTrivialGrid, theHeightZero, BinaryWord() );
	theSet03.remove( lowerPrimaryCellIntersection );
	expected_result_arr[0] = new GridCell( theTrivialGrid, theHeightOne, make_binary_word("00") );
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theSet03, 1 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );
	
	ARIADNE_PRINT_TEST_CASE_TITLE("Remove a GridCell (p.c. height=2) form a GridTreeSet(p.c. height=1): The cell is a subset.");
	GridCell higherPrimaryCellSubset( theTrivialGrid, theHeightTwo, make_binary_word("000011") );
	theSet04.remove( higherPrimaryCellSubset );
	expected_result_arr[0] = new GridCell( theTrivialGrid, theHeightTwo, make_binary_word("00000") );
	expected_result_arr[1] = new GridCell( theTrivialGrid, theHeightTwo, make_binary_word("000010") );
	expected_result_arr[2] = new GridCell( theTrivialGrid, theHeightTwo, make_binary_word("00111") );
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theSet04, 3 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );
	
	ARIADNE_PRINT_TEST_CASE_TITLE("Remove a GridCell (p.c. height=2) form a GridTreeSet(p.c. height=1): The cell does not intersect the set.");
	GridCell higherPrimaryCellNoIntersection( theTrivialGrid, theHeightTwo, make_binary_word("010") );
	theSet05.remove( higherPrimaryCellNoIntersection );
	expected_result_arr[0] = new GridCell( theTrivialGrid, theHeightTwo, make_binary_word("0000") );
	expected_result_arr[1] = new GridCell( theTrivialGrid, theHeightTwo, make_binary_word("00111") );
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theSet05, 2 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );
	
	ARIADNE_PRINT_TEST_CASE_TITLE("Remove a GridCell (p.c. height=2) form a GridTreeSet(p.c. height=1): The cell intersects the set.");
	GridCell higherPrimaryCellIntersection( theTrivialGrid, theHeightTwo, make_binary_word("0011") );
	theSet06.remove( higherPrimaryCellIntersection );
	expected_result_arr[0] = new GridCell( theTrivialGrid, theHeightTwo, make_binary_word("0000") );
	ARIADNE_TEST_GRID_TREE_SUBSET_ITERATOR( expected_result_arr, theSet06, 1 );
	ARIADNE_CLEAN_TEST__VECTOR( expected_result_arr );
}

void test_cell_subset_subset() {
		
	//Allocate a trivial Grid two dimensional grid
	Grid theTrivialGrid(2, 1.0);
	
	const uint smallHeight = 0;
	const uint mediumHeight = 1;
	const uint bigHeight = 2;

	//Create the cell, will be rooted to the primary cell mediumHeight
	GridCell theCell( theTrivialGrid, mediumHeight, make_binary_word("110") );

	//Create the binary tree, will be rooted to the primary cell of height bigHeight
	BinaryTreeNode * pBinaryTreeRoot = new BinaryTreeNode(true);
	pBinaryTreeRoot->split();
	pBinaryTreeRoot->left_node()->split();
	pBinaryTreeRoot->right_node()->split();
	pBinaryTreeRoot->left_node()->right_node()->set_disabled();
	pBinaryTreeRoot->right_node()->left_node()->set_disabled();
	
	pBinaryTreeRoot->left_node()->left_node()->split();
	pBinaryTreeRoot->left_node()->left_node()->right_node()->split();
	pBinaryTreeRoot->left_node()->left_node()->left_node()->split();
	pBinaryTreeRoot->left_node()->left_node()->right_node()->left_node()->set_disabled();
	pBinaryTreeRoot->left_node()->left_node()->left_node()->right_node()->set_disabled();
	
	//Create the GridTreeSubset, will be rooted to the primary cell of height smallHeight
	GridTreeSubset theSmallSubPaving( theTrivialGrid, smallHeight, make_binary_word("0011"),
						pBinaryTreeRoot->left_node()->left_node()->right_node()->right_node() );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test subset operation GridCell (height=1), GridTreeSubset.mince(2) (height=0)");
	theSmallSubPaving.mince(2);
	ARIADNE_PRINT_TEST_COMMENT("theCell");
	cout << theCell << endl;
	ARIADNE_PRINT_TEST_COMMENT("theSmallSubPaving");
	cout << theSmallSubPaving << endl;
	ARIADNE_TEST_EQUAL( subset( theCell, theSmallSubPaving), false );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test subset operation GridCell (height=1), GridTreeSubset.recombine() (height=0)");
	theSmallSubPaving.recombine();
	ARIADNE_PRINT_TEST_COMMENT("theCell");
	cout << theCell << endl;
	ARIADNE_PRINT_TEST_COMMENT("theSmallSubPaving");
	cout << theSmallSubPaving << endl;
	ARIADNE_TEST_EQUAL( subset( theCell, theSmallSubPaving), false );
	
	//Create the GridTreeSubset, will be rooted to the primary cell of height bigHeight
	GridTreeSubset theBigSubPaving( theTrivialGrid, bigHeight, make_binary_word("0011"),
						pBinaryTreeRoot->left_node()->left_node()->right_node()->right_node() );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test subset operation GridCell (height=1), GridTreeSubset.mince(2) (height=2)");
	theSmallSubPaving.mince(2);
	ARIADNE_PRINT_TEST_COMMENT("theCell");
	cout << theCell << endl;
	ARIADNE_PRINT_TEST_COMMENT("theBigSubPaving");
	cout << theBigSubPaving << endl;
	ARIADNE_TEST_EQUAL( subset( theCell, theBigSubPaving), true );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test subset operation GridCell (height=1), GridTreeSubset.recombine() (height=2)");
	theSmallSubPaving.recombine();
	ARIADNE_PRINT_TEST_COMMENT("theCell");
	cout << theCell << endl;
	ARIADNE_PRINT_TEST_COMMENT("theBigSubPaving");
	cout << theBigSubPaving << endl;
	ARIADNE_TEST_EQUAL( subset( theCell, theBigSubPaving), true );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test subset operation GridCell (height=1), GridTreeSubset.mince(2), left_node()->left_node()->set_disabled() (height=2)");
	theSmallSubPaving.mince(2);
	theSmallSubPaving.binary_tree()->left_node()->left_node()->set_disabled();
	ARIADNE_PRINT_TEST_COMMENT("theCell");
	cout << theCell << endl;
	ARIADNE_PRINT_TEST_COMMENT("theBigSubPaving");
	cout << theBigSubPaving << endl;
	ARIADNE_TEST_EQUAL( subset( theCell, theBigSubPaving), false );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test subset operation GridCell (height=1), GridTreeSubset.mince(2), left_node()->left_node()->set_enabled(), right_node()->make_leaf(false) (height=2)");
	theSmallSubPaving.mince(2);
	theSmallSubPaving.binary_tree()->left_node()->left_node()->set_enabled();
	theSmallSubPaving.binary_tree()->right_node()->make_leaf(false);
	ARIADNE_PRINT_TEST_COMMENT("theCell");
	cout << theCell << endl;
	ARIADNE_PRINT_TEST_COMMENT("theBigSubPaving");
	cout << theBigSubPaving << endl;
	ARIADNE_TEST_EQUAL( subset( theCell, theBigSubPaving), true );


	//Create the GridTreeSub, will be rooted to the primary cell of height smallHeight
	GridTreeSet theSmallPaving( theTrivialGrid, smallHeight, pBinaryTreeRoot->left_node()->left_node()->right_node()->right_node() );
	//Restore the binary tree
	theSmallPaving.binary_tree()->right_node()->set_enabled();
	theSmallPaving.binary_tree()->right_node()->split();
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test subset operation GridCell (height=1), GridTreeSet.mince(2) (after restoring the binary tree) (height=2)");
	theSmallPaving.mince(2);
	ARIADNE_PRINT_TEST_COMMENT("theCell");
	cout << theCell << endl;
	ARIADNE_PRINT_TEST_COMMENT("theSmallPaving");
	cout << theSmallPaving << endl;
	ARIADNE_TEST_EQUAL( subset( theCell, theSmallPaving), true );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test subset operation GridCell (height=1), GridTreeSet.recombine() (height=2)");
	theSmallPaving.recombine();
	ARIADNE_PRINT_TEST_COMMENT("theCell");
	cout << theCell << endl;
	ARIADNE_PRINT_TEST_COMMENT("theSmallPaving");
	cout << theSmallPaving << endl;
	ARIADNE_TEST_EQUAL( subset( theCell, theSmallPaving), true );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test subset operation GridCell (height=1), GridTreeSet.mince(2), left_node()->left_node()->set_disabled() (height=2)");
	theSmallPaving.mince(2);
	theSmallPaving.binary_tree()->left_node()->left_node()->set_disabled();
	ARIADNE_PRINT_TEST_COMMENT("theCell");
	cout << theCell << endl;
	ARIADNE_PRINT_TEST_COMMENT("theSmallPaving");
	cout << theSmallPaving << endl;
	ARIADNE_TEST_EQUAL( subset( theCell, theSmallPaving), false );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test subset operation GridCell (height=1), GridTreeSet.mince(2), left_node()->left_node()->set_enabled(), ...right_node()->make_leaf(false) (height=2)");
	theSmallPaving.mince(2);
	theSmallPaving.binary_tree()->left_node()->left_node()->set_enabled();
	theSmallPaving.binary_tree()->right_node()->make_leaf(false);
	ARIADNE_PRINT_TEST_COMMENT("theCell");
	cout << theCell << endl;
	ARIADNE_PRINT_TEST_COMMENT("theSmallPaving");
	cout << theSmallPaving << endl;
	ARIADNE_TEST_EQUAL( subset( theCell, theSmallPaving), true );
	
	
	//Create the GridTreeSub, will be rooted to the primary cell of height bigHeight
	BinaryTreeNode * pBinaryTreeRootCopy = new BinaryTreeNode( * pBinaryTreeRoot );
	GridTreeSet theBigPaving( theTrivialGrid, bigHeight, pBinaryTreeRootCopy );
	//Restore the binary tree
	pBinaryTreeRootCopy->left_node()->left_node()->right_node()->right_node()->right_node()->set_enabled();
	pBinaryTreeRootCopy->left_node()->left_node()->right_node()->right_node()->right_node()->split();
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test subset operation GridCell (height=1), GridTreeSet.mince(6) (after restoring the binary tree) (height=2)");
	theBigPaving.mince(6);
	ARIADNE_PRINT_TEST_COMMENT("theCell");
	cout << theCell << endl;
	ARIADNE_PRINT_TEST_COMMENT("theBigPaving");
	cout << theBigPaving << endl;
	ARIADNE_TEST_EQUAL( subset( theCell, theBigPaving), true );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test subset operation GridCell (height=1), GridTreeSet.recombine() (height=2)");
	theBigPaving.recombine();
	ARIADNE_PRINT_TEST_COMMENT("theCell");
	cout << theCell << endl;
	ARIADNE_PRINT_TEST_COMMENT("theBigPaving");
	cout << theBigPaving << endl;
	ARIADNE_TEST_EQUAL( subset( theCell, theBigPaving), true );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test subset operation GridCell (height=1), GridTreeSet.mince(6), ...left_node()->left_node()->set_disabled() (height=2)");
	theBigPaving.mince(6);
	theBigPaving.binary_tree()->left_node()->left_node()->right_node()->right_node()->left_node()->left_node()->set_disabled();
	ARIADNE_PRINT_TEST_COMMENT("theCell");
	cout << theCell << endl;
	ARIADNE_PRINT_TEST_COMMENT("theBigPaving");
	cout << theBigPaving << endl;
	ARIADNE_TEST_EQUAL( subset( theCell, theBigPaving), false );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Test subset operation GridCell (height=1), GridTreeSet.mince(6), ...left_node()->left_node()->set_enabled(), ...right_node()->make_leaf(false) (height=2)");
	theBigPaving.mince(6);
	theBigPaving.binary_tree()->left_node()->left_node()->right_node()->right_node()->left_node()->left_node()->set_enabled();
	theBigPaving.binary_tree()->left_node()->left_node()->right_node()->right_node()->right_node()->make_leaf(false);
	ARIADNE_PRINT_TEST_COMMENT("theCell");
	cout << theCell << endl;
	ARIADNE_PRINT_TEST_COMMENT("theBigPaving");
	cout << theBigPaving << endl;
	ARIADNE_TEST_EQUAL( subset( theCell, theBigPaving), true );
}

void test_subsets_join() {

	//Allocate a trivial Grid two dimensional grid
	Grid theTrivialGrid(2, 1.0);
	
	const uint smallHeight = 1;
	const uint bigHeight = 2;
	
	//Make set one
	BinaryWord tree = make_binary_word("1100100");
	BinaryWord leaves = make_binary_word("1010");
	BinaryTreeNode binaryTreeRootOne( tree, leaves );
	GridTreeSubset theSet1( theTrivialGrid, smallHeight, BinaryWord(), &binaryTreeRootOne );
	
	//Make set two
	tree = make_binary_word("1100100");
	leaves = make_binary_word("0101");
	BinaryTreeNode binaryTreeRootTwo( tree, leaves );
	GridTreeSubset theSet2( theTrivialGrid, bigHeight, make_binary_word("00"), &binaryTreeRootTwo );

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Join theSet1 and theSet2, when recombined they should give us the primary cell of height 1 rooted to the primary cell of height 2");
	ARIADNE_PRINT_TEST_COMMENT("theSet1");
	cout << theSet1 << endl;
	ARIADNE_PRINT_TEST_COMMENT("theSet2");
	cout << theSet2 << endl;
	GridTreeSet resultSet = join( theSet1, theSet2 );
	resultSet.recombine();
	GridTreeSet expectedResultSet( theTrivialGrid, bigHeight, make_binary_word("11000"), make_binary_word("1000") );
	ARIADNE_TEST_EQUAL( expectedResultSet, resultSet);
}

void test_subsets_intersection() {

	//Allocate a trivial Grid two dimensional grid
	Grid theTrivialGrid(2, 1.0);
	
	const uint smallHeight = 1;
	const uint bigHeight = 2;
	
	//Make set one
	BinaryWord tree = make_binary_word("1100100");
	BinaryWord leaves = make_binary_word("1010");
	BinaryTreeNode binaryTreeRootOne( tree, leaves );
	GridTreeSubset theSet1( theTrivialGrid, smallHeight, BinaryWord(), &binaryTreeRootOne );
	
	//Make set two
	tree = make_binary_word("1100100");
	leaves = make_binary_word("0111");
	BinaryTreeNode binaryTreeRootTwo( tree, leaves );
	GridTreeSubset theSet2( theTrivialGrid, bigHeight, make_binary_word("00"), &binaryTreeRootTwo );

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Intersect theSet1 and theSet2, should give us a set rooted to the primary cell of height 2");
	ARIADNE_PRINT_TEST_COMMENT("theSet1");
	cout << theSet1 << endl;
	ARIADNE_PRINT_TEST_COMMENT("theSet2");
	cout << theSet2 << endl;
	GridTreeSet resultSet = intersection( theSet1, theSet2 );
	GridTreeSet expectedResultSet( theTrivialGrid, bigHeight, make_binary_word("11110010000"), make_binary_word("001000") );
	ARIADNE_TEST_EQUAL( expectedResultSet, resultSet);
}

void test_subsets_difference() {

	//Allocate a trivial Grid two dimensional grid
	Grid theTrivialGrid(2, 1.0);
	
	const uint smallHeight = 1;
	const uint bigHeight = 2;
	
	//Make set one
	BinaryWord tree = make_binary_word("1100100");
	BinaryWord leaves = make_binary_word("1010");
	BinaryTreeNode binaryTreeRootOne( tree, leaves );
	GridTreeSubset theSet1( theTrivialGrid, smallHeight, BinaryWord(), &binaryTreeRootOne );
	
	//Make set two
	tree = make_binary_word("1100100");
	leaves = make_binary_word("0111");
	BinaryTreeNode binaryTreeRootTwo( tree, leaves );
	GridTreeSubset theSet2( theTrivialGrid, bigHeight, make_binary_word("00"), &binaryTreeRootTwo );

	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Removing theSet1 from theSet2, should give us a set rooted to the primary cell of height 2");
	ARIADNE_PRINT_TEST_COMMENT("theSet1");
	cout << theSet1 << endl;
	ARIADNE_PRINT_TEST_COMMENT("theSet2");
	cout << theSet2 << endl;
	GridTreeSet resultSet = difference( theSet1, theSet2 );
	GridTreeSet expectedResultSet( theTrivialGrid, bigHeight, make_binary_word("11110010000"), make_binary_word("100000") );
	ARIADNE_TEST_EQUAL( expectedResultSet, resultSet);
}

void test_cell_overlap_subset() {
	
	//Allocate a trivial Grid two dimensional grid
	Grid theTrivialGrid(2, 1.0);
	
	const uint heightZero = 0;
	const uint heightOne = 1;
	const uint heightTwo = 2;
	
	BinaryWord tree = make_binary_word("1111001000100");
	BinaryWord leaves = make_binary_word("1001001");
	//Create the set and the subset of this set, they are both rooted to the same primary node of heightTwo 
	GridTreeSet theSet( theTrivialGrid, heightTwo, new BinaryTreeNode( tree, leaves ) );
	//The subset is basically the zero level primary cell
	BinaryWord path = make_binary_word("11");
	GridTreeSubset theSubset( theTrivialGrid, heightOne, path, theSet.binary_tree()->left_node()->left_node()->right_node()->right_node() );

	GridCell theLowCell( theTrivialGrid, heightZero, make_binary_word("111") );		// does intersect with the set and the subset
	GridCell theMediumCellOne( theTrivialGrid, heightOne, make_binary_word("1000") );	// does not intersect with the set and the subset
	GridCell theMediumCellTwo( theTrivialGrid, heightOne, make_binary_word("0000") );	// does intersect with the set but not the subset
	GridCell theHighCellOne( theTrivialGrid, heightTwo, make_binary_word("11") );		// does intersect with the set but not the subset
	GridCell theHighCellTwo( theTrivialGrid, heightTwo, make_binary_word("01") );		// does not intersect with the set and the subset
	GridCell theHighCellThree( theTrivialGrid, heightTwo, make_binary_word("00") );		// does intersect with the set and the subset
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Testing bool overlap( const GridCell& , const GridTreeSubset& )");
	ARIADNE_PRINT_TEST_COMMENT("theSet");
	cout << theSet << endl;
	ARIADNE_PRINT_TEST_COMMENT("theLowCell");
	cout << theLowCell << endl;
	ARIADNE_TEST_EQUAL( overlap( theLowCell, theSet ), true );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Testing bool overlap( const GridCell& , const GridTreeSubset& )");
	ARIADNE_PRINT_TEST_COMMENT("theSet");
	cout << theSet << endl;
	ARIADNE_PRINT_TEST_COMMENT("theMediumCellOne");
	cout << theMediumCellOne << endl;
	ARIADNE_TEST_EQUAL( overlap( theMediumCellOne, theSet ), false );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Testing bool overlap( const GridCell& , const GridTreeSubset& )");
	ARIADNE_PRINT_TEST_COMMENT("theSet");
	cout << theSet << endl;
	ARIADNE_PRINT_TEST_COMMENT("theMediumCellTwo");
	cout << theMediumCellTwo << endl;
	ARIADNE_TEST_EQUAL( overlap( theMediumCellTwo, theSet ), true );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Testing bool overlap( const GridCell& , const GridTreeSubset& )");
	ARIADNE_PRINT_TEST_COMMENT("theSet");
	cout << theSet << endl;
	ARIADNE_PRINT_TEST_COMMENT("theHighCellOne");
	cout << theHighCellOne << endl;
	ARIADNE_TEST_EQUAL( overlap( theHighCellOne, theSet ), true );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Testing bool overlap( const GridCell& , const GridTreeSubset& )");
	ARIADNE_PRINT_TEST_COMMENT("theSet");
	cout << theSet << endl;
	ARIADNE_PRINT_TEST_COMMENT("");
	cout << theHighCellTwo << endl;
	ARIADNE_TEST_EQUAL( overlap( theHighCellTwo, theSet ), false );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Testing bool overlap( const GridCell& , const GridTreeSubset& )");
	ARIADNE_PRINT_TEST_COMMENT("theSet");
	cout << theSet << endl;
	ARIADNE_PRINT_TEST_COMMENT("theHighCellThree");
	cout << theHighCellThree << endl;
	ARIADNE_TEST_EQUAL( overlap( theHighCellThree, theSet ), true );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Testing bool overlap( const GridCell& , const GridTreeSet& )");
	ARIADNE_PRINT_TEST_COMMENT("theSubset");
	cout << theSubset << endl;
	ARIADNE_PRINT_TEST_COMMENT("theLowCell");
	cout << theLowCell << endl;
	ARIADNE_TEST_EQUAL( overlap( theLowCell, theSubset ), true );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Testing bool overlap( const GridCell& , const GridTreeSet& )");
	ARIADNE_PRINT_TEST_COMMENT("theSubset");
	cout << theSubset << endl;
	ARIADNE_PRINT_TEST_COMMENT("theMediumCellOne");
	cout << theMediumCellOne << endl;
	ARIADNE_TEST_EQUAL( overlap( theMediumCellOne, theSubset ), false );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Testing bool overlap( const GridCell& , const GridTreeSet& )");
	ARIADNE_PRINT_TEST_COMMENT("theSubset");
	cout << theSubset << endl;
	ARIADNE_PRINT_TEST_COMMENT("theMediumCellTwo");
	cout << theMediumCellTwo << endl;
	ARIADNE_TEST_EQUAL( overlap( theMediumCellTwo, theSubset ), false );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Testing bool overlap( const GridCell& , const GridTreeSet& )");
	ARIADNE_PRINT_TEST_COMMENT("theSubset");
	cout << theSubset << endl;
	ARIADNE_PRINT_TEST_COMMENT("theHighCellOne");
	cout << theHighCellOne << endl;
	ARIADNE_TEST_EQUAL( overlap( theHighCellOne, theSubset ), false );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Testing bool overlap( const GridCell& , const GridTreeSet& )");
	ARIADNE_PRINT_TEST_COMMENT("theSubset");
	cout << theSubset << endl;
	ARIADNE_PRINT_TEST_COMMENT("theHighCellTwo");
	cout << theHighCellTwo << endl;
	ARIADNE_TEST_EQUAL( overlap( theHighCellTwo, theSubset ), false );
	
	// !!!
	ARIADNE_PRINT_TEST_CASE_TITLE("Testing bool overlap( const GridCell& , const GridTreeSet& )");
	ARIADNE_PRINT_TEST_COMMENT("theSubset");
	cout << theSubset << endl;
	ARIADNE_PRINT_TEST_COMMENT("theHighCellThree");
	cout << theHighCellThree << endl;
	ARIADNE_TEST_EQUAL( overlap( theHighCellThree, theSubset ), true );
	
}

int main() {
	
	test_binary_tree();
	
	test_grid_paving_cursor();

	test_grid_paving_const_iterator();

	test_grid_paving_cell();
	
	test_grid_sub_paving();
	
	test_grid_paving();
	
	test_adjoin_operation_one();
	
	test_adjoin_operation_two();
	
	test_adjoin_operation_three();

	test_adjoin_outer_approximation_operation();

	test_restrict();
	
	test_remove_one();
	
	test_remove_two();
	
	test_cell_subset_subset();
	
	test_subsets_join();
	
	test_subsets_intersection();
	
	test_subsets_difference();
	
	test_cell_overlap_subset();
	
	return ARIADNE_TEST_FAILURES;
}

