/***************************************************************************
 *            test_grid_set.cc
 *
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

#include "set.h"
#include "grid_set.h"
#include "set_interface.h"


#include "test.h"

using namespace Ariadne;
using namespace std;

static int title_counter = 0;

void print_title(const char * pTitle){
	cout << endl << "***" << ++title_counter << ": "<< pTitle << "***" << endl;
	cout.flush();
}

void print_comment(const char * pComment){
	cout << "* COMMENT: " << pComment << "" << endl;
	cout.flush();
}

#define TEST_RESULT_STRING( pMsg, actual_result, expected_result ) \
{ \
	print_comment(pMsg); \
	cout << "* " ; \
	ARIADNE_TEST_EQUAL(actual_result, expected_result); \
	cout.flush(); \
}

#define TEST_BINARY_TREE( pMsg, expected_result, pBinaryTreeRoot) \
{ \
	string actual_result = (pBinaryTreeRoot)->tree_to_string(); \
	TEST_RESULT_STRING( pMsg, actual_result, expected_result ) \
}

#define TEST_BINARY_TREE_NODE( pMsg, expected_result, pBinaryTreeRoot) \
{ \
	string actual_result = (pBinaryTreeRoot)->node_to_string(); \
	TEST_RESULT_STRING( pMsg, actual_result, expected_result ) \
}

#define TEST_TO_STRING_RESULT( pMsg, expected_result, theToStringObject) \
{ \
	string actual_result = theToStringObject.to_string()+"\n"; \
	TEST_RESULT_STRING( pMsg, actual_result, expected_result ) \
}

#define TEST_BINARY_WORD( pMsg, expected_result, theToStringObject) \
{ \
	ostringstream actual_result_os; \
	actual_result_os << theBinaryPath; \
	string actual_result = actual_result_os.str(); \
	TEST_RESULT_STRING( pMsg, actual_result, expected_result ) \
}

void test_binary_tree() {
	string expected_result;
	
	// !!!
	print_title("Allocate an ebabled node and manipulate it");

	BinaryTreeNode * pBinaryTreeRoot = new BinaryTreeNode(true);
	expected_result = "+0";
	TEST_BINARY_TREE("Initial tree node", expected_result, pBinaryTreeRoot );

	pBinaryTreeRoot->set_disabled();
	expected_result = "isLeaf = 1, isEnabled = 0, isDisabled = 1";
	TEST_BINARY_TREE_NODE( "Disable the node", expected_result, pBinaryTreeRoot );

	print_comment("Making the leaf node intermediate should cause the IsALeafNodeException exception");
	ARIADNE_TEST_THROW( pBinaryTreeRoot->set_unknown(), IsALeafNodeException);
	expected_result = "isLeaf = 1, isEnabled = 0, isDisabled = 1";
	TEST_BINARY_TREE_NODE( "The node still have to be disabled", expected_result, pBinaryTreeRoot );
	
	pBinaryTreeRoot->set_enabled();
	expected_result = "isLeaf = 1, isEnabled = 1, isDisabled = 0";
	TEST_BINARY_TREE_NODE( "Enable the node", expected_result, pBinaryTreeRoot );
	
	// !!!
	print_title("Split the enabled node");
	print_comment("Should mark the node as intermediate and create two enabled subnodes");

	pBinaryTreeRoot->split();
	expected_result = "?1+01+0";
	TEST_BINARY_TREE("Splitted node tree", expected_result, pBinaryTreeRoot);
	expected_result = "isLeaf = 0, isEnabled = 0, isDisabled = 0";
	TEST_BINARY_TREE_NODE( "Splitted node", expected_result, pBinaryTreeRoot );

	// !!!
	print_title("Split a splited node");
	print_comment("We expect an unchanged tree, since it is already split");
	pBinaryTreeRoot->mince(1);
	expected_result = "?1+01+0";
	TEST_BINARY_TREE("Unchanged tree", expected_result, pBinaryTreeRoot);
	
	// !!!
	print_title("Copy the node by using the copy constructor");
	print_comment("The entire subtree should be copied");
	BinaryTreeNode theBinaryTreeCopy( *pBinaryTreeRoot );
	expected_result = "?1+01+0";
	TEST_BINARY_TREE("Copied tree", expected_result, &theBinaryTreeCopy);
	
	// !!!
	print_title("Recombine the new tree and see what happens");
	print_comment("The new tree should be reduced to a node");
	print_comment("The old tree should remain unchanged");
	theBinaryTreeCopy.recombine();
	expected_result = "+0";
	TEST_BINARY_TREE("Reduced tree", expected_result, &theBinaryTreeCopy);
	expected_result = "?1+01+0";
	TEST_BINARY_TREE("The original tree", expected_result, pBinaryTreeRoot);
	
	// !!!
	print_title("Disable the left tree leaf, recombine the tree and split the disabled leaf to level 3");
	pBinaryTreeRoot->left_node()->set_disabled();
	expected_result = "?1-01+0";
	TEST_BINARY_TREE("A tree with disabled left node", expected_result, pBinaryTreeRoot);
	pBinaryTreeRoot->recombine();
	expected_result = "?1-01+0";
	TEST_BINARY_TREE("A tree after recombination", expected_result, pBinaryTreeRoot);
	pBinaryTreeRoot->left_node()->mince(3);
	expected_result = "?1-01+0";
	TEST_BINARY_TREE("A tree after splitting the disabled node", expected_result, pBinaryTreeRoot);
	pBinaryTreeRoot->mince(3);
	expected_result = "?1-01?1?1+01+01?1+01+0";
	TEST_BINARY_TREE("A tree after splitting the root node", expected_result, pBinaryTreeRoot);
	print_comment("Testing the depth of the binary tree");
	ARIADNE_TEST_EQUAL(pBinaryTreeRoot->depth(), 3);

	// !!!
	print_title("Disable the left most enabled leaf, recombine the tree");
	pBinaryTreeRoot->right_node()->left_node()->left_node()->set_disabled();
	expected_result = "?1-01?1?1-01+01?1+01+0";
	TEST_BINARY_TREE("A tree with enabled left-most leaf", expected_result, pBinaryTreeRoot);
	pBinaryTreeRoot->recombine();
	expected_result = "?1-01?1?1-01+01+0";
	TEST_BINARY_TREE("A tree after recombination", expected_result, pBinaryTreeRoot);
	
	delete pBinaryTreeRoot;
	
	// !!!
	print_title("Testing the add_enabled( BinaryWord& ) method");
	pBinaryTreeRoot = new BinaryTreeNode(true);
	BinaryWord binaryWord;
	print_comment("Adding enabled node to an enabled node, the former one is defined by an empty path");
	pBinaryTreeRoot->add_enabled( binaryWord );
	expected_result = "+0";
	TEST_BINARY_TREE("The binary tree should stay intact, i.e. consist of one enabled node", expected_result, pBinaryTreeRoot);
	
	print_comment("Adding enabled sub node, path: (true,false) to a disabled node");
	pBinaryTreeRoot->set_disabled();
	binaryWord.push_back(true);
	binaryWord.push_back(false);
	pBinaryTreeRoot->add_enabled( binaryWord );
	expected_result = "?1-01?1+01-0";
	TEST_BINARY_TREE("The binary tree have five nodes, four leafs and one (true, false) is enabled", expected_result, pBinaryTreeRoot);

	print_comment("Adding enabled sub node, path: (true) to a disabled node");
	binaryWord = BinaryWord();
	binaryWord.push_back(true);
	pBinaryTreeRoot->add_enabled( binaryWord );
	expected_result = "?1-01+0";
	TEST_BINARY_TREE("The binary tree have three nodes, two leafs and one (true) is enabled", expected_result, pBinaryTreeRoot);

	print_comment("Adding enabled sub node, path: (true,false) to a disabled node");
	binaryWord = BinaryWord();
	binaryWord.push_back(false);
	binaryWord.push_back(true);
	pBinaryTreeRoot->add_enabled( binaryWord );
	expected_result = "?1?1-01+01+0";
	TEST_BINARY_TREE("The binary tree have five nodes, four leafs and two (false, true) and (true) are enabled", expected_result, pBinaryTreeRoot);
}

/*
* Combines the string output matching the output of the to_string method of GridCell
*/
string createPavingCellOutput(const string& grid, const string& heigth, const string& binary_word, const string& box){
	string result;
	result =  "\n The grid: Grid( "+grid+" )\n";
	result += " Primary root cell's height: "+heigth+"\n";
	result += " The path to the cell from the primary root cell: "+binary_word+"\n";
	result += " The cell's box: "+box+"\n";
	return result;
}

/*
* Combines the string output matching the output of the to_string method of GridTreeSubset
*/
string createSubPavingOutput(const string& grid, const string& heigth, const string& binary_word, const string& box, const string& tree){
	string result = createPavingCellOutput( grid, heigth, binary_word, box );
	result += " The subpaving's tree: "+tree+"\n";
	return result;
}

void test_grid_paving_cursor(){
	string expected_result;
	
	//Allocate the Grid
	Vector<Float> origin(3); origin.set(0,-0.25); origin.set(1,0.25); origin.set(2,1.5);
	Vector<Float> lengths(3); lengths.set(0, 0.25); lengths.set(1, 0.25); lengths.set(2, 0.25);
	Grid theGrid( origin, lengths );

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
	print_title("Check the created Cursor data on a small subpaving");
	expected_result = "\n The underlying subpaving:";
	expected_result += createSubPavingOutput("origin=[-0.25,0.25,1.5], lengths=[0.25,0.25,0.25]","2","01","[[-0.5:0],[0.5:1],[1.25:2.25]]","?1-01+0");
	expected_result += " The current stack index: 0\n";
	expected_result += " The current stack data: \n";
	expected_result += " Element 0: isLeaf = 0, isEnabled = 0, isDisabled = 0\n";
	expected_result += " The current grid cell: ";
	expected_result += createPavingCellOutput("origin=[-0.25,0.25,1.5], lengths=[0.25,0.25,0.25]","2","01","[[-0.5:0],[0.5:1],[1.25:2.25]]");
	TEST_TO_STRING_RESULT("The sub-tree based paving cursor (initial): ", expected_result, theGSPCursorSmall);

	// !!!
	print_title("Check the created Cursor data on a large subpaving");
	expected_result = "\n The underlying subpaving:";
	expected_result += createSubPavingOutput("origin=[-0.25,0.25,1.5], lengths=[0.25,0.25,0.25]", "2","e","[[-0.5:0.5],[0:1],[1.25:2.25]]","?1?1?1+01-01?1-01+01?1?1+01+01?1+01-0");
	expected_result += " The current stack index: 0\n";
	expected_result += " The current stack data: \n";
	expected_result += " Element 0: isLeaf = 0, isEnabled = 0, isDisabled = 0\n";
	expected_result += " The current grid cell: ";
	expected_result += createPavingCellOutput("origin=[-0.25,0.25,1.5], lengths=[0.25,0.25,0.25]","2","e","[[-0.5:0.5],[0:1],[1.25:2.25]]");
	TEST_TO_STRING_RESULT("The root-node based paving cursor (initial): ", expected_result, theGSPCursorLarge);
	
	// !!!
	print_title("Moving up from the Cursor's root node should cause the NotAllowedMoveException exception");
	ARIADNE_TEST_THROW( theGSPCursorSmall.move_up(), NotAllowedMoveException);
	
	// !!!
	print_title("Let's move the cursor of theGSPCursorLarge: left, right, we should be in the root of the smaller subpaving then");
	theGSPCursorLarge.move_left();
	theGSPCursorLarge.move_right();
	const GridTreeSubset tmpSubPaving = *theGSPCursorLarge;
	expected_result = createSubPavingOutput("origin=[-0.25,0.25,1.5], lengths=[0.25,0.25,0.25]","2","01","[[-0.5:0],[0.5:1],[1.25:2.25]]","?1-01+0");
	TEST_TO_STRING_RESULT("The cursor has moved and we have the following subpaving: ", expected_result, tmpSubPaving);
	expected_result = "\n The underlying subpaving:";
	expected_result += createSubPavingOutput("origin=[-0.25,0.25,1.5], lengths=[0.25,0.25,0.25]","2","e","[[-0.5:0.5],[0:1],[1.25:2.25]]","?1?1?1+01-01?1-01+01?1?1+01+01?1+01-0");
	expected_result += " The current stack index: 2\n";
	expected_result += " The current stack data: \n";
	expected_result += " Element 0: isLeaf = 0, isEnabled = 0, isDisabled = 0\n";
	expected_result += " Element 1: isLeaf = 0, isEnabled = 0, isDisabled = 0\n";
	expected_result += " Element 2: isLeaf = 0, isEnabled = 0, isDisabled = 0\n";
	expected_result += " The current grid cell: ";
	expected_result += createPavingCellOutput("origin=[-0.25,0.25,1.5], lengths=[0.25,0.25,0.25]","2","01","[[-0.5:0],[0.5:1],[1.25:2.25]]");
	TEST_TO_STRING_RESULT("The cursor's state is: ", expected_result, theGSPCursorLarge);
	
	// !!!
	print_title("Check the is_enabled/disabled/leaf/root functions, while moving the cursor in the large and small subpavings");
	ARIADNE_TEST_EQUAL(theGSPCursorSmall.is_root(), true);
	ARIADNE_TEST_EQUAL(theGSPCursorSmall.is_leaf(), false);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_root(), false);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_leaf(), false);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_enabled(), false);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_disabled(), false);
	print_comment("Moving to the left, to a disabled node");
	theGSPCursorLarge.move_left();
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_leaf(), true);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_root(), false);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_enabled(), false);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_disabled(), true);
	print_comment("Moving up and to the right, to an enabled node");
	theGSPCursorLarge.move_up();
	theGSPCursorLarge.move_right();
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_leaf(), true);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_root(), false);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_enabled(), true);
	ARIADNE_TEST_EQUAL(theGSPCursorLarge.is_disabled(), false);
	
	// !!!
	print_title("Test that moving to the left and right with the leaf nodes causes the NotAllowedMoveException exception");
	print_comment("Get to the leaf node on the left");
	theGSPCursorSmall.move_left();
	print_comment("Test the present state of theGSPCursorSmall");
	expected_result = "\n The underlying subpaving:";
	expected_result += createSubPavingOutput("origin=[-0.25,0.25,1.5], lengths=[0.25,0.25,0.25]","2","01","[[-0.5:0],[0.5:1],[1.25:2.25]]","?1-01+0");
	expected_result += " The current stack index: 1\n";
	expected_result += " The current stack data: \n";
	expected_result += " Element 0: isLeaf = 0, isEnabled = 0, isDisabled = 0\n";
	expected_result += " Element 1: isLeaf = 1, isEnabled = 0, isDisabled = 1\n";
	expected_result += " The current grid cell: ";
	expected_result += createPavingCellOutput("origin=[-0.25,0.25,1.5], lengths=[0.25,0.25,0.25]","2","010","[[-0.5:0],[0.5:1],[1.25:1.75]]");
	TEST_TO_STRING_RESULT("The cursor's state is: ", expected_result, theGSPCursorSmall);
	print_comment("Try to get to the node on the left");
	ARIADNE_TEST_THROW(theGSPCursorSmall.move_left(), NotAllowedMoveException );
	print_comment("Try to get to the node on the right");
	ARIADNE_TEST_THROW(theGSPCursorSmall.move_right(), NotAllowedMoveException );
	print_comment("The state of theGSPCursorSmall is supposed to remain unchanged");
	TEST_TO_STRING_RESULT("The cursor's state is: ", expected_result, theGSPCursorSmall);

	// !!!
	print_title("Test the set_enabled/set_disabled methods of the cursor");
	delete pRootTreeNode;
	pRootTreeNode = new BinaryTreeNode(true);
	print_comment("Mince the enabled binary tree node to level 2");
	pRootTreeNode->mince(2);
	GridTreeSet theGridTreeSet( theGrid, theHeight, pRootTreeNode );
	//Create the Cursor for the GridTreeSubset
	GridTreeCursor theGPCursor( &theGridTreeSet );
	//Move to the leaf node
	print_comment("Move the cursor to the right most leaf and check it");
	theGPCursor.move_right();
	theGPCursor.move_right();
	ARIADNE_TEST_EQUAL(theGPCursor.is_leaf(), true);
	ARIADNE_TEST_EQUAL(theGPCursor.is_root(), false);
	ARIADNE_TEST_EQUAL(theGPCursor.is_enabled(), true);
	ARIADNE_TEST_EQUAL(theGPCursor.is_disabled(), false);
	print_comment("Set the leaf as disabled");
	theGPCursor.set_disabled();
	ARIADNE_TEST_EQUAL(theGPCursor.is_disabled(), true);
	print_comment("Set the leaf (back) as enabled");
	theGPCursor.set_enabled();
	ARIADNE_TEST_EQUAL(theGPCursor.is_enabled(), true);
	print_comment("Move one node up, to a non-leaf node and try to enable it, the NotALeafNodeException should be thrown");
	theGPCursor.move_up();
	ARIADNE_TEST_THROW( theGPCursor.set_disabled(), NotALeafNodeException);
}

void test_iterator( const string expected_result[], const GridTreeSubset & theGridTreeSubset, const int expected_number_elements ){
	int elements_count = 0;
	for (GridTreeSubset::const_iterator it = theGridTreeSubset.begin(), end = theGridTreeSubset.end(); it != end; it++, elements_count++) {
		TEST_TO_STRING_RESULT("The next iterator node is: ", expected_result[elements_count], (*it) );
		//cout<< (*it).to_string() << endl;
	}
	print_comment("Test that we iterated through the right number of nodes");
	ARIADNE_TEST_EQUAL(elements_count, expected_number_elements);
}

void test_grid_paving_const_iterator(){
	string expected_result[8];
	
	//Allocate the Grid
	Vector<Float> origin(2); origin.set(0, -0.25); origin.set(1, 0.25);
	Vector<Float> lengths(2); lengths.set(0, 0.25); lengths.set(1, 0.25);
	Grid theGrid( origin, lengths );
	
	//Define the higth of the primary root cell.
	const uint theHeight = 2;

	//Create the binary tree;
	BinaryTreeNode * pRootTreeNode = new BinaryTreeNode(true);

	//Create the GridTreeSubset
	GridTreeSubset theGridSubPavingLarge( theGrid, theHeight, BinaryWord(), pRootTreeNode );

	// !!!
	print_title("Test the sequence in which GridPavingIterator goes through the tree leafs ");
	print_comment("The tree depth is 0 and all leaf nodes are enabled");

	expected_result[0] = createPavingCellOutput( "origin=[-0.25,0.25], lengths=[0.25,0.25]","2","e","[[-0.5:0.5],[0:1]]");

	test_iterator( expected_result, theGridSubPavingLarge, 1 );

	// !!!
	print_title("Test the sequence in which GridPavingIterator goes through the tree leafs ");
	print_comment("The tree depth is 1 and all leaf nodes are enabled");

	//Mince the tree to the certain higth
	pRootTreeNode->mince(1);

	expected_result[0] = createPavingCellOutput( "origin=[-0.25,0.25], lengths=[0.25,0.25]","2","0","[[-0.5:0],[0:1]]");

	expected_result[1] = createPavingCellOutput( "origin=[-0.25,0.25], lengths=[0.25,0.25]","2","1","[[0:0.5],[0:1]]");

	test_iterator( expected_result, theGridSubPavingLarge, 2 );

	// !!!
	print_title("Test the sequence in which GridPavingIterator goes through the tree leafs ");
	print_comment("The tree depth is 2 and all leaf nodes are enabled");

	//Mince the tree to the certain higth
	pRootTreeNode->mince(2);

	expected_result[0] = createPavingCellOutput( "origin=[-0.25,0.25], lengths=[0.25,0.25]","2","00","[[-0.5:0],[0:0.5]]");

	expected_result[1] = createPavingCellOutput( "origin=[-0.25,0.25], lengths=[0.25,0.25]","2","01","[[-0.5:0],[0.5:1]]");

	expected_result[2] = createPavingCellOutput( "origin=[-0.25,0.25], lengths=[0.25,0.25]","2","10","[[0:0.5],[0:0.5]]");

	expected_result[3] = createPavingCellOutput( "origin=[-0.25,0.25], lengths=[0.25,0.25]","2","11","[[0:0.5],[0.5:1]]");

	test_iterator( expected_result, theGridSubPavingLarge, 4 );

	// !!!
	print_title("Test the sequence in which GridPavingIterator goes through the tree leafs ");
	print_comment("The tree depth is 3 and all leaf nodes are enabled");

	//Mince the tree to the certain higth
	pRootTreeNode->mince(3);

	expected_result[0] = createPavingCellOutput( "origin=[-0.25,0.25], lengths=[0.25,0.25]","2","000","[[-0.5:-0.25],[0:0.5]]");

	expected_result[1] = createPavingCellOutput( "origin=[-0.25,0.25], lengths=[0.25,0.25]","2","001","[[-0.25:0],[0:0.5]]");

	expected_result[2] = createPavingCellOutput( "origin=[-0.25,0.25], lengths=[0.25,0.25]","2","010","[[-0.5:-0.25],[0.5:1]]");

	expected_result[3] = createPavingCellOutput( "origin=[-0.25,0.25], lengths=[0.25,0.25]","2","011","[[-0.25:0],[0.5:1]]");

	expected_result[4] = createPavingCellOutput( "origin=[-0.25,0.25], lengths=[0.25,0.25]","2","100","[[0:0.25],[0:0.5]]");

	expected_result[5] = createPavingCellOutput( "origin=[-0.25,0.25], lengths=[0.25,0.25]","2","101","[[0.25:0.5],[0:0.5]]");

	expected_result[6] = createPavingCellOutput( "origin=[-0.25,0.25], lengths=[0.25,0.25]","2","110","[[0:0.25],[0.5:1]]");

	expected_result[7] = createPavingCellOutput( "origin=[-0.25,0.25], lengths=[0.25,0.25]","2","111","[[0.25:0.5],[0.5:1]]");

	test_iterator( expected_result, theGridSubPavingLarge, 8 );
	
	// !!!
	print_title("Disable some of the leaf nodes and test GridPavingIterator");
	print_comment("The tree depth is 3 and some tree nodes are disabled");

	//Mark some nodes as disabled
	pRootTreeNode->left_node()->left_node()->left_node()->set_disabled();
	pRootTreeNode->left_node()->right_node()->right_node()->set_disabled();
	pRootTreeNode->right_node()->left_node()->left_node()->set_disabled();
	pRootTreeNode->right_node()->left_node()->right_node()->set_disabled();

	//Reuse the previous result strings
	string expected_result_tmp[4];
	expected_result_tmp[0] = expected_result[1];
	expected_result_tmp[1] = expected_result[2];
	expected_result_tmp[2] = expected_result[6];
	expected_result_tmp[3] = expected_result[7];

	test_iterator( expected_result_tmp, theGridSubPavingLarge, 4 );

	// !!!
	print_title("Recombine the tree and test GridPavingIterator");
	print_comment("The tree depth is 3 but some nodes are at depth 2");
	
	//Recombine the paving
	theGridSubPavingLarge.recombine();
	
	//Reuse some of the previous result strings
	expected_result_tmp[0] = expected_result[1];
	expected_result_tmp[1] = expected_result[2];
	expected_result_tmp[2] = createPavingCellOutput( "origin=[-0.25,0.25], lengths=[0.25,0.25]","2","11","[[0:0.5],[0.5:1]]");
	
	test_iterator( expected_result_tmp, theGridSubPavingLarge, 3 );

	// !!!
	print_title("Mince the tree back to level 3, enable/disable some nodes, recombine and and test GridPavingIterator");
	print_comment("The tree depth is 3 but some nodes are at depth 2");

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
	expected_result_tmp[0] = expected_result[1];
	expected_result_tmp[1] = expected_result[2];
	expected_result_tmp[2] = createPavingCellOutput( "origin=[-0.25,0.25], lengths=[0.25,0.25]","2","10","[[0:0.5],[0:0.5]]");
	expected_result_tmp[3] = expected_result[7];

	test_iterator( expected_result_tmp, theGridSubPavingLarge, 4 );

	// !!!
	print_title("Test how the constant Cursor can be retrieved from the Constant iterator");
	GridTreeSubset::const_iterator it = theGridSubPavingLarge.begin();
	const GridTreeCursor theGPCursor = it.cursor();
	ARIADNE_TEST_EQUAL(theGPCursor.is_leaf(), true);
	ARIADNE_TEST_EQUAL(theGPCursor.is_root(), false);
	ARIADNE_TEST_EQUAL(theGPCursor.is_enabled(), true);
	ARIADNE_TEST_EQUAL(theGPCursor.is_disabled(), false);
	print_comment("Set the leaf as disabled");
	theGPCursor.set_disabled();
	ARIADNE_TEST_EQUAL(theGPCursor.is_disabled(), true);
	print_comment("Set the leaf (back) as enabled");
	theGPCursor.set_enabled();
	ARIADNE_TEST_EQUAL(theGPCursor.is_enabled(), true);
}

void test_grid_sub_paving(){
	string expected_result;
	
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
	print_title("Test Mincing operations of GridTreeSubset on the one dimensional Grid");
	expected_result =  createSubPavingOutput("origin=[0], lengths=[0.5]","2","01","[[0:0.5]]","+0");
	TEST_TO_STRING_RESULT("The initial paving: ", expected_result, theGridSPOneDim);
	
	//Mince to the level two from the root of the paving, this should give
	//us 4 leaf nodes with intervals of width 0.125 (in the original space)
	theGridSPOneDim.mince(2);

	print_comment("Minced the sub-paving to depth 2");
	string expected_result_arr[4];

	expected_result_arr[0] =  createPavingCellOutput( "origin=[0], lengths=[0.5]","2","0100","[[0:0.125]]");

	expected_result_arr[1] =  createPavingCellOutput( "origin=[0], lengths=[0.5]","2","0101","[[0.125:0.25]]");

	expected_result_arr[2] =  createPavingCellOutput( "origin=[0], lengths=[0.5]","2","0110","[[0.25:0.375]]");

	expected_result_arr[3] =  createPavingCellOutput( "origin=[0], lengths=[0.5]","2","0111","[[0.375:0.5]]");

	test_iterator( expected_result_arr, theGridSPOneDim, 4 );
	
	print_comment("Recombine and subdivide the sub-paving to cell width 1.1");
	theGridSPOneDim.recombine();
	theGridSPOneDim.subdivide(1.1);
	TEST_TO_STRING_RESULT("We should have a single node as in the initial sub paving: ", expected_result, theGridSPOneDim);

	print_comment("Subdivide the sub-paving to cell width 0.4, this should give us two sub cells");
	string expected_result_arr_tmp[2];
	expected_result_arr_tmp[0] =  createPavingCellOutput( "origin=[0], lengths=[0.5]","2","010","[[0:0.25]]");

	expected_result_arr_tmp[1] =  createPavingCellOutput( "origin=[0], lengths=[0.5]","2","011","[[0.25:0.5]]");

	theGridSPOneDim.subdivide(0.4);

	test_iterator( expected_result_arr_tmp, theGridSPOneDim, 2 );
	
	print_comment("Subdivide the sub-paving to cell width 0.126, this should give us four sub cells");
	theGridSPOneDim.subdivide(0.126);

	test_iterator( expected_result_arr, theGridSPOneDim, 4 );

	print_comment("Recombine and Subdivide the sub-paving to cell width 0.126, this should give us four sub cells");
	theGridSPOneDim.recombine();
	theGridSPOneDim.subdivide(0.126);

	test_iterator( expected_result_arr, theGridSPOneDim, 4 );

	// !!!
	print_title("Test Mincing operations of GridTreeSubset on the two dimensional Grid");
	//Allocate the Grid, one Dimension
	Vector<Float> originTwo(2); originTwo.set(0, -0.25); originTwo.set(1, 0.5);
	Vector<Float> lengthsTwo(2); lengthsTwo.set(0, 0.25); lengthsTwo.set(1, 0.5);
	const Grid theTwoDimGrid( originTwo, lengthsTwo );

	//Create the GridTreeSubset
	GridTreeSubset theGridSPTwoDim( theTwoDimGrid, theHeight, thePathToSubPavingRoot, pRootTreeNode );

	print_comment("Recombine and Mince the sub-paving to depth 2, this should give us four sub cells");
	theGridSPTwoDim.recombine();
	theGridSPTwoDim.mince(2);

	expected_result_arr[0] =  createPavingCellOutput( "origin=[-0.25,0.5], lengths=[0.25,0.5]","2","0100","[[-0.5:-0.25],[1:1.5]]");

	expected_result_arr[1] =  createPavingCellOutput( "origin=[-0.25,0.5], lengths=[0.25,0.5]","2","0101","[[-0.5:-0.25],[1.5:2]]");

	expected_result_arr[2] =  createPavingCellOutput( "origin=[-0.25,0.5], lengths=[0.25,0.5]","2","0110","[[-0.25:0],[1:1.5]]");

	expected_result_arr[3] =  createPavingCellOutput( "origin=[-0.25,0.5], lengths=[0.25,0.5]","2","0111","[[-0.25:0],[1.5:2]]");

	test_iterator( expected_result_arr, theGridSPTwoDim, 4 );

	print_comment("Recombine and Subdivide the sub-paving to cell width 0.5, this should give us four sub cells");
	//At this moment the coordinate cell widths are: for x -- 0.25 and for y -- 0.5 
	theGridSPTwoDim.recombine();
	theGridSPTwoDim.subdivide(0.51);

	test_iterator( expected_result_arr, theGridSPTwoDim, 4 );

	print_comment("Subdivide the sub-paving to cell width 0.4, this should give us sixteen sub cells");
	//At this moment the coordinate cell widths are: for x -- 0.25 and for y -- 0.5 
	theGridSPTwoDim.subdivide(0.4);

	expected_result = createSubPavingOutput("origin=[-0.25,0.5], lengths=[0.25,0.5]","2","01","[[-0.5:0],[1:2]]",
						"?1?1?1?1+01+01?1+01+01?1?1+01+01?1+01+01?1?1?1+01+01?1+01+01?1?1+01+01?1+01+0");
	TEST_TO_STRING_RESULT("The sub paving with 16 leaf nodes: ", expected_result, theGridSPTwoDim);
	
	print_comment("The depth of the sub-paving's tree should be 4");
	ARIADNE_TEST_EQUAL( theGridSPTwoDim.depth(), 4 );
}

void test_grid_paving(){
	string expected_result;
	
	// !!!
	print_title("Test allocation of a trivial GridTreeSet");
	GridTreeSet * pTrivialPaving = new GridTreeSet(4, true);
	expected_result = createSubPavingOutput("origin=[0,0,0,0], lengths=[1,1,1,1]", "0", "e", "[[0:1],[0:1],[0:1],[0:1]]", "+0");
	TEST_TO_STRING_RESULT("A trivial paving for [0,1]x[0,1]x[0,1]x[0,1], enabled cell: ", expected_result, (*pTrivialPaving) );

	// !!!
	print_title("Test GridTreeSet copy constructor");
	GridTreeSet theTrivialPaving( ( *pTrivialPaving ) );
	pTrivialPaving->mince(1);
	string expected_result_minced = createSubPavingOutput("origin=[0,0,0,0], lengths=[1,1,1,1]", "0", "e", "[[0:1],[0:1],[0:1],[0:1]]", "?1+01+0");
	TEST_TO_STRING_RESULT("Minced trivial paving for [0,1]x[0,1]x[0,1]x[0,1], enabled cell: ", expected_result_minced, (*pTrivialPaving) );
	TEST_TO_STRING_RESULT("A copy of the original paving, should stay unchanged: ", expected_result, theTrivialPaving );

	// !!!
	print_title("Test GridTreeSet (Grid, Box) constructor");
	expected_result = createSubPavingOutput("origin=[-0.25,0.5], lengths=[0.25,0.5]", "4", "e", "[[-1.5:2.5],[-2:6]]", "-0");
	//Allocate the Grid, one Dimension
	const Grid theTwoDimGrid(make_vector<Float>("[-0.25,0.5]"),make_vector<Float>("[0.25,0.5]"));
	//Note: the box is related to the grid, but not to the original space
	GridTreeSet theTwoDimPaving( theTwoDimGrid, make_box("[0,1.5]x[-1.5,3.5]") );
	TEST_TO_STRING_RESULT("The resulting GridTreeSet: ", expected_result, theTwoDimPaving );

	// !!!
	print_title("Test GridTreeSet (Grid, Height, BooleanArray, BooleanArray) constructor");
	expected_result = createSubPavingOutput("origin=[-0.25,0.5], lengths=[0.25,0.5]", "2", "e", "[[-0.5:0.5],[0:2]]", "?1?1+01-01?1?1+01-01+0");
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
	TEST_TO_STRING_RESULT("The resulting GridTreeSet: ", expected_result, theTwoDimTreePaving );
	
}

void test_grid_paving_cell(){
	string expected_result;

	// !!!
	print_title("Test the static methods of the GridCell");
	BinaryWord theBinaryPath;
	
	expected_result = "e";
	theBinaryPath = GridCell::primary_cell_path( 1, 0, 0 );
	TEST_BINARY_WORD( "Dimension: 1, topCellHeight: 0, bottomCellHeight: 0", expected_result , theBinaryPath );
	
	expected_result = "1";
	theBinaryPath = GridCell::primary_cell_path( 1, 1, 0 );
	TEST_BINARY_WORD( "Dimension: 1, topCellHeight: 1, bottomCellHeight: 0", expected_result , theBinaryPath );
	
	expected_result = "01";
	theBinaryPath = GridCell::primary_cell_path( 1, 2, 0 );
	TEST_BINARY_WORD( "Dimension: 1, topCellHeight: 2, bottomCellHeight: 0", expected_result , theBinaryPath );
	
	expected_result = "0";
	theBinaryPath = GridCell::primary_cell_path( 1, 2, 1 );
	TEST_BINARY_WORD( "Dimension: 1, topCellHeight: 2, bottomCellHeight: 1", expected_result , theBinaryPath );
	
	expected_result = "e";
	theBinaryPath = GridCell::primary_cell_path( 1, 2, 2 );
	TEST_BINARY_WORD( "Dimension: 1, topCellHeight: 2, bottomCellHeight: 2", expected_result , theBinaryPath );
}

void test_adjoin_operation_one(){
	string expected_result;
	
	//Allocate a trivial Grid
	Grid theGrid(2, 1.0);
	
	// !!!
	print_title("Test adjoining a GridCell to the GridTreeSet");
	//Define the GridCell that is rooted to a high primary cell
	const int theHigherCellHeight = 2;
	BinaryWord theHigherCellPath;
	theHigherCellPath.push_back(true);
	theHigherCellPath.push_back(false);
	theHigherCellPath.push_back(true);
	theHigherCellPath.push_back(false);
	GridCell theHigherLevelCell( theGrid, theHigherCellHeight, theHigherCellPath );

	print_comment("The GridCell with the primary root cell height = 2");
	expected_result = createPavingCellOutput("origin=[0,0], lengths=[1,1]","2","1010","[[2:3],[-1:0]]");
	TEST_TO_STRING_RESULT("The initial GridCell: ", expected_result, theHigherLevelCell );
	
	//Define the higth of the primary root cell.
	//Create the GridTreeSet with the box is related to the grid, but not to the original space
	GridTreeSet theOneCellPaving( theGrid, true );
	expected_result = createSubPavingOutput("origin=[0,0], lengths=[1,1]", "0", "e", "[[0:1],[0:1]]", "+0");
	print_comment("The GridTreeSet with the primary root cell height = 0");
	TEST_TO_STRING_RESULT("The initial GridTreeSet: ", expected_result, theOneCellPaving );
	
	theOneCellPaving.adjoin( theHigherLevelCell );
	expected_result = createSubPavingOutput("origin=[0,0], lengths=[1,1]", "2", "e", "[[-1:3],[-1:3]]", "?1?1?1-01?1-01+01-01?1?1-01?1+01-01-0");
	TEST_TO_STRING_RESULT("The GridTreeSet after adding the cell: ", expected_result, theOneCellPaving );
}

void test_adjoin_operation_two(){
	string expected_result;
	
	//Allocate a trivial Grid
	Grid theGrid(2, 1.0);
	
	// !!!
	print_title("Test adjoining a GridCell to the GridTreeSet");
	//Define the GridCell that is rooted to the lower primary cell
	const int theLowerCellHeight = 1;
	BinaryWord theLowerCellPath;
	theLowerCellPath.push_back(true);
	theLowerCellPath.push_back(true);
	GridCell theLowerLevelCell( theGrid, theLowerCellHeight, theLowerCellPath );

	print_comment("The GridCell with the primary root cell height = 1");
	expected_result = createPavingCellOutput("origin=[0,0], lengths=[1,1]","1","11","[[0:1],[0:1]]");
	TEST_TO_STRING_RESULT("The initial GridCell: ", expected_result, theLowerLevelCell );
	
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
	expected_result = createSubPavingOutput("origin=[0,0], lengths=[1,1]", "2", "e", "[[-1:3],[-1:3]]", "?1-01?1?1-01?1+01-01-0");
	print_comment("The GridTreeSet with the primary root cell height = 2");
	TEST_TO_STRING_RESULT("The initial GridTreeSet: ", expected_result, theOneCellPaving );
	
	theOneCellPaving.adjoin( theLowerLevelCell );
	expected_result = createSubPavingOutput("origin=[0,0], lengths=[1,1]", "2", "e", "[[-1:3],[-1:3]]", "?1?1?1-01?1-01+01-01?1?1-01?1+01-01-0");
	TEST_TO_STRING_RESULT("The GridTreeSet after adding the cell: ", expected_result, theOneCellPaving );
}

void test_adjoin_operation_three(){
	string expected_result;
	
	//Allocate a trivial Grid
	Grid theGrid(2, 1.0);
	
	// !!!
	print_title("Test adjoining a GridCell to the GridTreeSet");
	//Define the GridCell that is rooted to the same primary cell
	const int theLowerCellHeight = 2;
	BinaryWord theLowerCellPath;
	theLowerCellPath.push_back(false);
	theLowerCellPath.push_back(false);
	theLowerCellPath.push_back(true);
	theLowerCellPath.push_back(true);
	GridCell theLowerLevelCell( theGrid, theLowerCellHeight, theLowerCellPath );

	print_comment("The GridCell with the primary root cell height = 2");
	expected_result = createPavingCellOutput("origin=[0,0], lengths=[1,1]","2","0011","[[0:1],[0:1]]");
	TEST_TO_STRING_RESULT("The initial GridCell: ", expected_result, theLowerLevelCell );
	
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
	expected_result = createSubPavingOutput("origin=[0,0], lengths=[1,1]", "2", "e", "[[-1:3],[-1:3]]", "?1-01?1?1-01?1+01-01-0");
	print_comment("The GridTreeSet with the primary root cell height = 2");
	TEST_TO_STRING_RESULT("The initial GridTreeSet: ", expected_result, theOneCellPaving );
	
	theOneCellPaving.adjoin( theLowerLevelCell );
	expected_result = createSubPavingOutput("origin=[0,0], lengths=[1,1]", "2", "e", "[[-1:3],[-1:3]]", "?1?1?1-01?1-01+01-01?1?1-01?1+01-01-0");
	TEST_TO_STRING_RESULT("The GridTreeSet after adding the cell: ", expected_result, theOneCellPaving );
}

void test_adjoin_outer_approximation_operation(){
	string expected_result;
	
	//Allocate a trivial Grid
	Grid theTrivialGrid(2, 1.0);
	
	// !!!
	print_title("Test adjoining_outer_approximation a SetInterface to the GridTreeSet");
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
	
	//Create the GridTreeSet with the box is related to the grid, but not to the original space
	GridTreeSet theOneCellPaving( theTrivialGrid, theHeight, pRootTreeNode );
	expected_result = createSubPavingOutput("origin=[0,0], lengths=[1,1]", "2", "e", "[[-1:3],[-1:3]]", "?1-01?1?1-01?1+01-01-0");
	print_comment("The GridTreeSet with the primary root cell height = 2");
	TEST_TO_STRING_RESULT("The initial GridTreeSet: ", expected_result, theOneCellPaving );
	
	theOneCellPaving.adjoin_outer_approximation( static_cast<const LocatedSetInterface&>(initialRectangle), 2 );
	
	string tree = "?1?1?1-01?1-01?1?1?1?1?1?1-01+01?1-01+01?1?1+01+01?1+01+01?1?1?1-01+01?1-01+01?1?1+01+01";
	tree += "?1+01+01?1?1?1?1+01-01?1+01-01-01?1?1?1+01-01?1+01-01-01?1?1?1?1?1-01+01?1-01+01?1?1+01+01";
	tree += "?1+01+01?1+01-01?1?1?1?1+01-01?1+01-01-01-01-01-0";
	expected_result = createSubPavingOutput("origin=[0,0], lengths=[1,1]", "4", "e", "[[-5:11],[-5:11]]", tree);
	TEST_TO_STRING_RESULT("The GridTreeSet after adding the cell: ", expected_result, theOneCellPaving );

	print_comment("Recombined GridTreeSet after adding the cell: ");
	string expected_result_arr[16];
	expected_result_arr[0] = createPavingCellOutput( "origin=[0,0], lengths=[1,1]","4","0011000001","[[-1:-0.5],[-0.5:0]]");
	expected_result_arr[1] = createPavingCellOutput( "origin=[0,0], lengths=[1,1]","4","0011000011","[[-0.5:0],[-0.5:0]]");
	expected_result_arr[2] = createPavingCellOutput( "origin=[0,0], lengths=[1,1]","4","00110001","[[-1:0],[0:1]]");
	expected_result_arr[3] = createPavingCellOutput( "origin=[0,0], lengths=[1,1]","4","0011001001","[[0:0.5],[-0.5:0]]");
	expected_result_arr[4] = createPavingCellOutput( "origin=[0,0], lengths=[1,1]","4","0011001011","[[0.5:1],[-0.5:0]]");
	expected_result_arr[5] = createPavingCellOutput( "origin=[0,0], lengths=[1,1]","4","00110011","[[0:1],[0:1]]");
	expected_result_arr[6] = createPavingCellOutput( "origin=[0,0], lengths=[1,1]","4","0011010000","[[-1:-0.5],[1:1.5]]");
	expected_result_arr[7] = createPavingCellOutput( "origin=[0,0], lengths=[1,1]","4","0011010010","[[-0.5:0],[1:1.5]]");
	expected_result_arr[8] = createPavingCellOutput( "origin=[0,0], lengths=[1,1]","4","0011011000","[[0:0.5],[1:1.5]]");
	expected_result_arr[9] = createPavingCellOutput( "origin=[0,0], lengths=[1,1]","4","0011011010","[[0.5:1],[1:1.5]]");
	expected_result_arr[10] = createPavingCellOutput( "origin=[0,0], lengths=[1,1]","4","0011100001","[[1:1.5],[-0.5:0]]");
	expected_result_arr[11] = createPavingCellOutput( "origin=[0,0], lengths=[1,1]","4","0011100011","[[1.5:2],[-0.5:0]]");
	expected_result_arr[12] = createPavingCellOutput( "origin=[0,0], lengths=[1,1]","4","00111001","[[1:2],[0:1]]");
	expected_result_arr[13] = createPavingCellOutput( "origin=[0,0], lengths=[1,1]","4","00111010","[[2:3],[-1:0]]");
	expected_result_arr[14] = createPavingCellOutput( "origin=[0,0], lengths=[1,1]","4","0011110000","[[1:1.5],[1:1.5]]");
	expected_result_arr[15] = createPavingCellOutput( "origin=[0,0], lengths=[1,1]","4","0011110010","[[1.5:2],[1:1.5]]");
	theOneCellPaving.recombine();
	test_iterator( expected_result_arr, theOneCellPaving, 16 );
	
	// !!!
	print_title("Create an outer_approximation of the rectangle on the scaling grid and get the GridTreeSet");
	Grid theScalingGrid(2, 2.0);
	GridTreeSet theOuterApproxGridTreeSet = outer_approximation( static_cast<LocatedSetInterface&>(initialRectangle), theScalingGrid, 2 );
	//IVAN S. ZAPREEV
	//NOTE: The recombination is needed because in the scaling Grid doing
	//	outer_approximation( theScalingGrid, initialRectangle, 2 )
	//will result in subdivisions that are equal to unit cell in the original
	//space in this case we will get much more elements, e.g. the cells
	//  [-1,0]x[0,2], [0,2]x[0,2] will be subdivided as well
	theOuterApproxGridTreeSet.recombine();
	expected_result_arr[0] = createPavingCellOutput( "origin=[0,0], lengths=[2,2]","3","11000011","[[-1:0],[-1:0]]");
	expected_result_arr[1] = createPavingCellOutput( "origin=[0,0], lengths=[2,2]","3","1100011","[[-1:0],[0:2]]");
	expected_result_arr[2] = createPavingCellOutput( "origin=[0,0], lengths=[2,2]","3","11001001","[[0:1],[-1:0]]");
	expected_result_arr[3] = createPavingCellOutput( "origin=[0,0], lengths=[2,2]","3","11001011","[[1:2],[-1:0]]");
	expected_result_arr[4] = createPavingCellOutput( "origin=[0,0], lengths=[2,2]","3","110011","[[0:2],[0:2]]");
	test_iterator( expected_result_arr, theOuterApproxGridTreeSet, 5 );
}

void test_restrict() {
	string expected_result_arr[3];
	
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
	print_title("Test restrict operation: GridTreeSet1.restrict( GridTreeSet2 )");
	print_comment("The initial GridTreeSet1: ");
	expected_result_arr[0] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","2","00","[[-1:1],[-1:1]]");
	expected_result_arr[1] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","2","1000","[[1:2],[-1:0]]");
	expected_result_arr[2] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","2","1011","[[2:3],[0:1]]");
	test_iterator( expected_result_arr, theThreeCellPavingH2, 3 );
	print_comment("The initial GridTreeSet2: ");
	expected_result_arr[0] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","2","100","[[1:2],[-1:1]]");
	expected_result_arr[1] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","2","1010","[[2:3],[-1:0]]");
	test_iterator( expected_result_arr, theTwoCellPavingH2, 2 );
	print_comment("The result after restrict: ");
	theThreeCellPavingH2.restrict( theTwoCellPavingH2 );
	expected_result_arr[0] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","2","1000","[[1:2],[-1:0]]");
	test_iterator( expected_result_arr, theThreeCellPavingH2, 1 );

	// !!!
	print_title("Test restrict operation: GridTreeSet.restrict( GridTreeSubset )");
	print_comment("The initial GridTreeSet: ");
	expected_result_arr[0] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","2","100","[[1:2],[-1:1]]");
	expected_result_arr[1] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","2","1010","[[2:3],[-1:0]]");
	test_iterator( expected_result_arr, theTwoCellPavingH2, 2 );
	print_comment("The initial GridTreeSubset: ");
	expected_result_arr[0] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","3","1100","[[-1:1],[-1:1]]");
	expected_result_arr[1] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","3","111000","[[1:2],[-1:0]]");
	expected_result_arr[2] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","3","111011","[[2:3],[0:1]]");
	test_iterator( expected_result_arr, theThreeCellSubPavingH3, 3 );
	print_comment("The result after restrict: ");
	theTwoCellPavingH2.restrict( theThreeCellSubPavingH3 );
	expected_result_arr[0] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","3","111000","[[1:2],[-1:0]]");
	test_iterator( expected_result_arr, theTwoCellPavingH2, 1 );
	
	//TODO: Test the case when the GridTreeSet has primary cell of the level 3
	//	The GridTreeSubset is at level 1 and it's primary cell is at level 2
}

void test_remove() {
	string expected_result_arr[3];
	
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
	print_title("Test remove operation: GridTreeSet1.remove( GridTreeSet2 )");
	print_comment("The initial GridTreeSet1: ");
	expected_result_arr[0] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","2","00","[[-1:1],[-1:1]]");
	expected_result_arr[1] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","2","1000","[[1:2],[-1:0]]");
	expected_result_arr[2] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","2","1011","[[2:3],[0:1]]");
	test_iterator( expected_result_arr, theThreeCellPavingH2, 3 );
	print_comment("The initial GridTreeSet2: ");
	expected_result_arr[0] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","2","100","[[1:2],[-1:1]]");
	expected_result_arr[1] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","2","1010","[[2:3],[-1:0]]");
	test_iterator( expected_result_arr, theTwoCellPavingH2, 2 );
	print_comment("The result after removal: ");
	theThreeCellPavingH2.remove( theTwoCellPavingH2 );
	expected_result_arr[0] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","2","00","[[-1:1],[-1:1]]");
	expected_result_arr[1] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","2","1011","[[2:3],[0:1]]");
	test_iterator( expected_result_arr, theThreeCellPavingH2, 2 );

	// !!!
	print_title("Test remove operation: GridTreeSet.remove( GridTreeSubset )");
	print_comment("The initial GridTreeSet: ");
	expected_result_arr[0] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","2","100","[[1:2],[-1:1]]");
	expected_result_arr[1] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","2","1010","[[2:3],[-1:0]]");
	test_iterator( expected_result_arr, theTwoCellPavingH2, 2 );
	print_comment("The initial GridTreeSubset: ");
	expected_result_arr[0] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","3","1100","[[-1:1],[-1:1]]");
	expected_result_arr[1] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","3","111000","[[1:2],[-1:0]]");
	expected_result_arr[2] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","3","111011","[[2:3],[0:1]]");
	test_iterator( expected_result_arr, theThreeCellSubPavingH3, 3 );
	print_comment("The result after remove: ");
	theTwoCellPavingH2.remove( theThreeCellSubPavingH3 );
	expected_result_arr[0] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","3","111001","[[1:2],[0:1]]");
	expected_result_arr[1] = createPavingCellOutput("origin=[0,0], lengths=[1,1]","3","111010","[[2:3],[-1:0]]");
	test_iterator( expected_result_arr, theTwoCellPavingH2, 2 );
	
	//TODO: Test the case when the GridTreeSet has primary cell of the level 3
	//	The GridTreeSubset is at level 1 and it's primary cell is at level 2
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
	
	test_remove();
	
	return ARIADNE_TEST_FAILURES;
}

