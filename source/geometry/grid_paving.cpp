/***************************************************************************
 *            grid_paving.cpp
 *
 *  Copyright  2008-12  Ivan S. Zapreev, Pieter Collins
 *
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "../function/functional.hpp"
#include "../config.hpp"

#include <iostream>
#include <iomanip>

#include "../utility/macros.hpp"
#include "../utility/exceptions.hpp"
#include "../utility/stlio.hpp"
#include "../geometry/function_set.hpp"
#include "../geometry/list_set.hpp"
#include "../geometry/grid_paving.hpp"
#include "../geometry/binary_tree.hpp"

#include "../geometry/set_interface.hpp"


namespace Ariadne {

typedef SizeType SizeType;

Bool subset(const GridCell& theCell, const GridTreeSubpaving& theSet);
GridTreePaving intersection( const GridTreeSubpaving& theSet1, const GridTreeSubpaving& theSet2 );
Bool subset( const GridCell& theCellOne, const GridCell& theCellTwo,
             BinaryWord * pPathPrefixOne, BinaryWord * pPathPrefixTwo, Nat * pPrimaryCellHeight);


//****************************************GridAbstractCell*******************************************/

LatticeBoxType GridAbstractCell::primary_cell_lattice_box( const Nat theHeight, const dimension_type dimensions ) {
    Int leftBottomCorner = 0, rightTopCorner = 1;
    //The zero level coordinates are known, so we need to iterate only for higher level primary cells
    for(Nat i = 1; i <= theHeight; i++){
        primary_cell_at_height(i, leftBottomCorner, rightTopCorner);
    }
    // 1.2 Constructing and return the box defining the primary cell (relative to the grid).

    return Vector< ExactIntervalType >( dimensions, ExactIntervalType( leftBottomCorner, rightTopCorner ) );
}

Nat GridAbstractCell::smallest_enclosing_primary_cell_height( const LatticeBoxType& theLatticeBoxType ) {
    const dimension_type dimensions = theLatticeBoxType.size();
    Int leftBottomCorner = 0, rightTopCorner = 1;
    Nat height = 0;
    //The zero level coordinates are known, so we need to iterate only for higher level primary cells
    do{
        //Check if the given box is a subset of a primary cell
        Box< ExactIntervalType > primaryCellBoxType( dimensions, ExactIntervalType( leftBottomCorner, rightTopCorner ) );
        if(  inside(theLatticeBoxType, primaryCellBoxType ) ) {
            //If yes then we are done
            break;
        }
        //Otherwise increase the height and recompute the new borders
        primary_cell_at_height( ++height, leftBottomCorner, rightTopCorner);
    }while( true );

    return height;
}

Nat GridAbstractCell::smallest_enclosing_primary_cell_height( const UpperBoxType& theBoxType, const Grid& theGrid) {
    Vector<ExactIntervalType> theLatticeBoxType( theBoxType.size() );
    //Convert the box to theGrid coordinates
    for( Nat i = 0; i != theBoxType.size(); ++i ) {
        theLatticeBoxType[i] = cast_exact_interval( ( theBoxType[i] - cast_exact(theGrid.origin()[i])) / cast_exact(theGrid.lengths()[i]) );
    }
    //Compute and return the smallest primary cell, enclosing this box on the grid
    return smallest_enclosing_primary_cell_height( theLatticeBoxType );
}

ExactBoxType GridAbstractCell::lattice_box_to_space(const LatticeBoxType & theLatticeBoxType, const Grid& theGrid ){
    const DimensionType dimensions = theGrid.dimension();
    ExactBoxType theTmpBoxType( dimensions );

    Vector<FloatDP> theGridOrigin( theGrid.origin() );
    Vector<FloatDP> theGridLengths( theGrid.lengths() );

    for(Nat current_dimension = 0; current_dimension < dimensions; current_dimension ++){
        const FloatDP theDimLength = theGridLengths[current_dimension];
        const FloatDP theDimOrigin = theGridOrigin[current_dimension];
        //Recompute the new dimension coordinates, detaching them from the grid
        //Compute lower and upper bounds separately, and then set the box lower
        //and upper values simultaneously to prevent lower temporarily higher than upper.
        FloatDP lower = add(approx, theDimOrigin, mul(approx, theDimLength, theLatticeBoxType[current_dimension].lower().raw() ) );
        FloatDP upper = add(approx, theDimOrigin, mul(approx, theDimLength, theLatticeBoxType[current_dimension].upper().raw() ) );
        theTmpBoxType[current_dimension].set(cast_exact(lower),cast_exact(upper));
    }

    return theTmpBoxType;
}

BinaryWord GridAbstractCell::primary_cell_path( const DimensionType dimensions, const Nat topPCellHeight, const Nat bottomPCellHeight) {
    BinaryWord theBinaryPath;

    //The path from one primary cell to another consists of alternating subsequences
    //of length \a dimensions. These subsequences consist either of ones or zeroes.
    //Odd primary cell height means that the first subsequence will consist of all
    //ones. Even primary cell height indicates the the first subsequence will consis
    //of all zeroes. This is due to the way we do the space subdivisions.
    if( topPCellHeight > bottomPCellHeight ){

        for( Nat i = topPCellHeight; i > bottomPCellHeight; i-- ){
            Bool odd_height = (i % 2) != 0;
            for( Nat j = 0; j < dimensions; j++ ){
                theBinaryPath.push_back( odd_height );
            }
        }
    }

    return theBinaryPath;
}

//TODO: Perhaps we can make it more efficient by reducing the height of the words till
//the minimal primary cell height and then comparing them by height and binary words
Bool GridAbstractCell::compare_abstract_grid_cells(const GridAbstractCell * pCellLeft, const GridAbstractCell &cellRight, const Nat comparator ) {
    ARIADNE_ASSERT( pCellLeft->_theGrid == cellRight._theGrid );
    const BinaryWord * pThisWord, * pOtherWord;
    //This is the temporary word to be used if heights are not equal
    BinaryWord rootNodePath;

    if( pCellLeft->_theHeight == cellRight._theHeight ) {
        //if the primary cells are of the same height, then we just compare the original binary words.
        pThisWord = & (pCellLeft->_theWord);
        pOtherWord = & (cellRight._theWord);
    } else {
        //Otherwise we have to re-root the cell with the lowest primary cell
        //to the highest primary cell and then compare the words again,
        if( pCellLeft->_theHeight > cellRight._theHeight ){
            rootNodePath = primary_cell_path( pCellLeft->_theGrid.dimension() , pCellLeft->_theHeight, cellRight._theHeight );
            rootNodePath.append( cellRight._theWord );
            pThisWord = & (pCellLeft->_theWord);
            pOtherWord = & (rootNodePath);
        } else {
            rootNodePath = primary_cell_path( pCellLeft->_theGrid.dimension() , cellRight._theHeight, pCellLeft->_theHeight );
            rootNodePath.append( pCellLeft->_theWord );
            pThisWord = & (rootNodePath);
            pOtherWord = & (cellRight._theWord);
        }
    }
    switch( comparator ){
        case COMPARE_EQUAL : return (*pThisWord) == (*pOtherWord);
        case COMPARE_LESS : return (*pThisWord) < (*pOtherWord);
        default:
            throw InvalidInput("The method's comparator argument should be either GridAbstractCell::COMPARE_EQUAL or GridAbstractCell::COMPARE_LESS.");
    }
}

//********************************************GridCell***********************************************/

GridCell GridCell::split(Bool isRight) const {
    BinaryWord theWord = _theWord;
    theWord.push_back( isRight );
    return GridCell( _theGrid, _theHeight, theWord);
}

Pair<GridCell,GridCell> GridCell::split() const {
    BinaryWord theLeftWord = _theWord;
    theLeftWord.push_back( false );
    BinaryWord theRightWord = _theWord;
    theRightWord.push_back( true );
    return std::make_pair(GridCell( _theGrid, _theHeight, theLeftWord),GridCell( _theGrid, _theHeight, theRightWord));
}

GridCell GridCell::smallest_enclosing_primary_cell( const UpperBoxType& theBoxType, const Grid& theGrid) {
    //Create the GridCell, corresponding to the smallest primary cell, enclosing this box
    return GridCell( theGrid, smallest_enclosing_primary_cell_height(theBoxType, theGrid), BinaryWord() );
}

//Computes the box corresponding the the cell defined by the primary cell and the binary word.
//The resulting box is not relaterd to the original space, but is a lattice box.
// 1. Compute the primary cell located the the height \a theHeight above the zero level,
// 2. Compute the cell defined by the path \a theWord (from the primary cell).
LatticeBoxType GridCell::compute_lattice_box( const DimensionType dimensions, const Nat theHeight, const BinaryWord& theWord ) {
    LatticeBoxType theResultLatticeBoxType( primary_cell_lattice_box( theHeight , dimensions ) );

    //2. Compute the cell on some grid, corresponding to the binary path from the primary cell.
    Nat current_dimension = 0;
    for(Nat i = 0; i < theWord.size(); i++){
        //We move through the dimensions in a linear fasshion
        current_dimension = i % dimensions;
        //Compute the middle point of the box's projection onto
        //the dimension \a current_dimension (relative to the grid)
        FloatDP middlePointInCurrDim = theResultLatticeBoxType[current_dimension].midpoint().raw();
        if( theWord[i] ){
            //Choose the right half
            theResultLatticeBoxType[current_dimension].set_lower( cast_exact(middlePointInCurrDim) );
        } else {
            //Choose the left half
            theResultLatticeBoxType[current_dimension].set_upper( cast_exact(middlePointInCurrDim) );
        }
    }
    return theResultLatticeBoxType;
}

//This method appends \a dimension() zeroes to the binary word defining this cell
//and return a GridOpenCell created with the given grid, primany cell height and
//the newly created word for the low-left cell of the open cell
GridOpenCell GridCell::interior() const {
    BinaryWord theOpenCellWord = _theWord;
    for( Nat i = 0; i < _theGrid.dimension(); i++) {
        theOpenCellWord.push_back(false);
    }
    //The open cell will be defined by the given new word, i.e. the path to the
    //left-bottom sub-quadrant cell but the box and the rest will be the same
    return GridOpenCell( _theGrid, _theHeight, theOpenCellWord, _theBoxType );
}

//NOTE: The cell defined byt the method's arguments is called the base cell.
//NOTE: Here we work with the lattice boxes that are in the grid
GridCell GridCell::neighboringCell( const Grid& theGrid, const Nat theHeight, const BinaryWord& theWord, const DimensionType dim ) {
    const DimensionType dimensions = theGrid.dimension();
    //1. Extend the base cell in the given dimension (dim) by it's half width. This way
    //   we are sure that we get a box that overlaps with the required neighboring cell.
    //NOTE: This box is in the original space, but not on the lattice
    Vector<ExactIntervalType> baseCellBoxInLattice =  GridCell::compute_lattice_box( dimensions, theHeight, theWord );
    const FloatDP upperBorderOverlapping = add(approx, baseCellBoxInLattice[dim].upper().raw(), hlf( baseCellBoxInLattice[dim].width().value_raw() ) );

    //2. Now check if the neighboring cell can be rooted to the given primary cell. For that
    //   we simply use the box computed in 1. and get the primary cell that encloses it.
    //NOTE: In fact, we only need to take about the upper border, because the lower border does not change.
    Int leftBottomCorner = 0, rightTopCorner = 1; Nat height = 0;
    do{
        if( upperBorderOverlapping <= rightTopCorner ) {
            //As soon as we fall into the primary cell we are done
            break;
        }
        //Otherwise increase the height and recompute the new borders
        primary_cell_at_height( ++height, leftBottomCorner, rightTopCorner);
    } while(true);

    //3. If it can not then we take the lowest required primary cell to root this cell to and
    //   to re-route the given base cell of GridOpenCell to that one.
    Nat theBaseCellHeight = theHeight;
    BinaryWord theBaseCellWord = theWord;
    if( height > theBaseCellHeight ) {
        //If we need a higher primary cell then extend the height and the word for the case cell
        theBaseCellWord = primary_cell_path( dimensions, height, theBaseCellHeight );
        theBaseCellWord.append( theWord );
        theBaseCellHeight = height;
    }

    //4. We need to start from the end of the new (extended) word representing the base cell. Then
    //   we go backwards and look for the smallest cell such that the upper border computed in 1.
    //   is less than the cell's box upper border (in the given dimension). This is indicated by incountering
    //   the first zero in the path to the base cell in dimension dim from the end of the path
    Nat position = theBaseCellWord.size() - 1;
    while(true) {
        //Only consider the dimension that we need and look for the first opotrunity to invert the path suffix.
        if( ( position % dimensions == dim ) && !theBaseCellWord[position] ) {
            break;
        }
        if (position == 0)
            break;
        position--;
    }

    //5. When this entry in the word is found from that point on we have to inverse the path in such
    //   a way that every component in the dimension from this point till the end of the word is
    //   inverted. This will provide us with the path to the neighborind cell in the given dimension
    for(Nat index = position; index < theBaseCellWord.size(); index++){
        if( index % dimensions == dim ) {
            //If this element of the path corresponds to the needed dimension then we need to invert it
            theBaseCellWord[index] = !theBaseCellWord[index];
        }
    }

    return GridCell( theGrid, theBaseCellHeight, theBaseCellWord );
}

DimensionType GridCell::dimension() const {
    return this->grid().dimension();
}

//************************************FRIENDS OF GridCell*****************************************/

Bool subset(const GridCell& theCellOne, const GridCell& theCellTwo, BinaryWord * pPathPrefixOne, BinaryWord * pPathPrefixTwo, Nat *pPrimaryCellHeight ) {
    //Test that the Grids are equal
    ARIADNE_ASSERT( theCellOne.grid() == theCellTwo.grid() );

    //Test that the binary words are empty, otherwise the results of computations are undefined
    Bool isDefault = false;
    if( ( pPathPrefixOne == nullptr ) && ( pPathPrefixTwo == nullptr ) ) {
        //If both parameters are null, then the user does not
        //need this data, so we allocate tepmorary objects
        pPathPrefixOne = new BinaryWord();
        pPathPrefixTwo = new BinaryWord();
        isDefault = true;
    } else {
        if( ( pPathPrefixOne == nullptr ) || ( pPathPrefixTwo == nullptr ) ) {
            //TODO: Find some better way to notify the user that
            //only one of the required pointers is not nullptr.
            ARIADNE_ASSERT( false );
        } else {
            //DO NOTHING: Both parameters are not nullptr, so it
            //is user data and isDefault should stay false.
        }
    }

    //Check if theCellOne is a subset of theCellTwo:

    //01 Align the cell's primary cells by finding path prefixes
    Nat primary_cell_height;
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
    if( pPrimaryCellHeight != nullptr ){
        *pPrimaryCellHeight = primary_cell_height;
    }

    //02 Add the rest of the paths to the path prefixes to get the complete paths
    pPathPrefixOne->append( theCellOne.word() );
    pPathPrefixTwo->append( theCellTwo.word() );

    //03 theCellOne is a subset of theCellTwo if pathPrefixTwo is a prefix of pathPrefixOne
    Bool result = pPathPrefixTwo->is_prefix( *pPathPrefixOne );

    //04 If working with default parameters, then release memory
    if( isDefault ) {
        delete pPathPrefixOne;
        delete pPathPrefixTwo;
    }

    return result;
}

OutputStream& operator<<(OutputStream& os, const GridCell& gridPavingCell){
    //Write the grid data to the string stream
    return os << "GridCell( " << gridPavingCell.grid() <<
        ", Primary cell height: " << gridPavingCell.height() <<
        ", Path to the root: " << gridPavingCell.word() <<
        ", ExactBoxType: " << gridPavingCell.box() << " )";
}


//********************************************GridOpenCell******************************************/

//NOTE: In this method first works with the boxes on the lattice, to make
//      computation exact, and then maps them to the original space.
ExactBoxType GridOpenCell::compute_box(const Grid& theGrid, const Nat theHeight, const BinaryWord& theWord) {
    Vector<ExactIntervalType> baseCellBoxInLattice =  GridCell::compute_lattice_box( theGrid.dimension(), theHeight, theWord );
    ExactBoxType openCellBoxInLattice( theGrid.dimension() );

    //Go through all the dimensions, and double the box size in the positive axis direction.
    for( DimensionType dim = 0; dim < theGrid.dimension(); dim++){
        ExactIntervalType openCellBoxInLatticeDimIntervalType;
        ExactIntervalType baseCellBoxInLatticeDimIntervalType = baseCellBoxInLattice[dim];
        FloatDP lower = baseCellBoxInLatticeDimIntervalType.lower().raw();
        FloatDP upper = add( near, baseCellBoxInLatticeDimIntervalType.upper().raw(),
                             sub( near, baseCellBoxInLatticeDimIntervalType.upper().raw(),
                                        baseCellBoxInLatticeDimIntervalType.lower().raw() ) );
        openCellBoxInLatticeDimIntervalType.set(cast_exact(lower),cast_exact(upper));

        openCellBoxInLattice[dim] = openCellBoxInLatticeDimIntervalType;
    }

    return lattice_box_to_space( openCellBoxInLattice, theGrid );
}

GridOpenCell GridOpenCell::split(TernaryChild child) const {
    BinaryWord theNewBaseCellPath = _theWord;
    Nat theNewPrimaryCellHeight;
    if( child == TernaryChild::MIDDLE ) {
        //Return the middle open cell (child == TernaryChild::MIDDLE)
        theNewBaseCellPath.push_back( true );
        theNewPrimaryCellHeight = _theHeight;
    } else {
        if( child == TernaryChild::RIGHT ) {
            //Return the right-most open cell (child == TernaryChild::RIGHT)
            //1. First determine in which dimention we are going to split
            //NOTE: We use theNewBaseCellPath.size() but not theNewBaseCellPath.size()-1 because
            //we want to determine the dimension in which will be the next split, but not the
            //dimension in which we had the last split.
            const DimensionType dim = theNewBaseCellPath.size() % _theGrid.dimension();
            //2. Then get the neighboring cell in this dimension
            GridCell neighboringCell = GridCell::neighboringCell( _theGrid, _theHeight, theNewBaseCellPath, dim );
            //3. Get the neighboring's cell height and word, because heigh might change
            theNewBaseCellPath = neighboringCell.word();
            theNewPrimaryCellHeight = neighboringCell.height();
            //4. Take the left half of the new base cell
            theNewBaseCellPath.push_back( false );
        } else {
            //Return the left-most open cell (child == TernaryChild::LEFT)
            theNewBaseCellPath.push_back( false );
            theNewPrimaryCellHeight = _theHeight;
        }
    }

    //Construct the new open cell and return it
    return GridOpenCell( _theGrid, theNewPrimaryCellHeight, theNewBaseCellPath );
}

GridOpenCell * GridOpenCell::smallest_open_subcell( const GridOpenCell &theOpenCell, const ExactBoxType & theBoxType ) {
    GridOpenCell * pSmallestOpenCell = nullptr;
    //If the box of the givel open cell covers the given box
    if( definitely( theOpenCell.box().covers( theBoxType ) ) ) {
        //First check the left sub cell
        pSmallestOpenCell = smallest_open_subcell( theOpenCell.split( TernaryChild::LEFT ), theBoxType );
        if( pSmallestOpenCell == nullptr ) {
            //If the left sub cell does not cover theBoxType then check the middle sub cell
            pSmallestOpenCell = smallest_open_subcell( theOpenCell.split( TernaryChild::MIDDLE ), theBoxType );
            if( pSmallestOpenCell == nullptr ) {
                //If the middle sub cell does not cover the box then check the right subcell
                pSmallestOpenCell = smallest_open_subcell( theOpenCell.split( TernaryChild::RIGHT ), theBoxType );
                if( pSmallestOpenCell == nullptr ) {
                    //If the right sub cell does not cover the box then the open cell
                    //theOpenCell is the smallest GridOpenCell covering theBoxType.
                    pSmallestOpenCell = new GridOpenCell( theOpenCell.grid(), theOpenCell.height(),
                                                          theOpenCell.word(), theOpenCell.box() );
                }
            }
        }
    }

    //Return the result, is null If the box of theOpenCell does not cover the given box
    return pSmallestOpenCell;
}

GridOpenCell GridOpenCell::outer_approximation( const ExactBoxType & theBoxType, const Grid& theGrid ) {
    //01. First we find the smallest primary GridCell that contain the given ExactBoxType.
    GridCell thePrimaryCell = GridCell::smallest_enclosing_primary_cell( theBoxType, theGrid );

    //02. Second we start subdividing it to find out the root cell for the smallest open cell containing theBoxType
    //NOTE: smallest_open_subcell returns nullptr or a pointer to new open cell, but here we are sure that
    //we can not get nullptr because of the choice of thePrimaryCell, therefore we do the following:
    GridOpenCell * pOpenCell = GridOpenCell::smallest_open_subcell( thePrimaryCell.interior(), theBoxType );
    GridOpenCell theOpenCell = (* pOpenCell); delete pOpenCell; //Deallocate the memory, to avoid the memory leaks
    return theOpenCell;
}

GridTreePaving GridOpenCell::closure() const {
    //01. First we compute the height of the primary cell that encloses the given open cell
    const Nat newHeight = smallest_enclosing_primary_cell_height( _theBoxType, _theGrid );

    //02. Re-route (if needed) the base cell to the new primary cell
    Nat theBaseCellHeight = _theHeight;
    BinaryWord theBaseCellWord = _theWord;
    //If we need a higher primary cell then extend the height and the word for the base cell
    if( newHeight > theBaseCellHeight ) {
        theBaseCellWord = primary_cell_path( _theGrid.dimension(), newHeight, theBaseCellHeight );
        theBaseCellWord.append( _theWord );
        theBaseCellHeight = newHeight;
    }

    //03. Allocate the resulting GridTreeSet with the root at the needed height
    GridTreePaving theResultSet( _theGrid, theBaseCellHeight, new BinaryTreeNode( false ) );

    //04. The preparations are done, now we need to add the base cell to the
    //    resulting GridTreeSet and to compute and add the other neighboring cells.
    BinaryWord tmpWord;
    neighboring_cells( theResultSet.root_cell().height(), theBaseCellWord, tmpWord, theResultSet );

    return theResultSet;
}

Void GridOpenCell::neighboring_cells( const Nat theHeight, const BinaryWord& theBaseCellWord,
                                      BinaryWord& cellPosition, GridTreePaving& theResultSet ) const {
    if( cellPosition.size() < _theGrid.dimension() ) {
        //Choose the left direction in the current dimension
        cellPosition.push_back( false );
        GridOpenCell::neighboring_cells( theHeight, theBaseCellWord, cellPosition, theResultSet );
        cellPosition.pop_back( );
        //Choose the right direction in the current dimension
        cellPosition.push_back( true );
        GridOpenCell::neighboring_cells( theHeight, theBaseCellWord, cellPosition, theResultSet );
        cellPosition.pop_back( );
    } else {
        //We have constructed the cell position relative to the base cell
        //for the case of _theGrid.dimension() dimensional space, now it
        //is time to compute this cell and add it to theResultSet
        theResultSet.adjoin( GridOpenCell::neighboring_cell( _theGrid, theHeight, theBaseCellWord, cellPosition ) );
    }
}

GridCell GridOpenCell::neighboring_cell( const Grid& theGrid, const Nat theHeight,
                                         const BinaryWord& theBaseCellWord, BinaryWord& cellPosition ) {
    const Nat num_dimensions = theGrid.dimension();
    //01. Allocate the Array of size _theGrid.dimensions() in which we will store
    //    the position in the path theBaseCellWord, for each dimension, from which on
    //    we need to inverse the path to get the proper neighboring cell.
    Vector<Nat> invert_position(num_dimensions);
    const Nat NO_INVERSE_POSITION = static_cast<Nat>(theBaseCellWord.size());
    //Initialize the Array with NO_INVERSE_POSITION to make sure that the inversion positions
    //for the dimensions that are not set to one in cellPosition will be undefined. Also,
    //count the required number of iverse dimensions
    Nat inverseDimensionsNumericType = 0;
    for( Nat i = 0; i < num_dimensions; i++ ) {
        invert_position[i] = NO_INVERSE_POSITION;
        inverseDimensionsNumericType = inverseDimensionsNumericType + cellPosition[i];
    }

    //02. Create the path to the neighboring cell and initialize it with the path to the base cell
    BinaryWord theNeighborCellWord = theBaseCellWord;

    //03. We need to start from the end of the new (extended) word representing the base cell. Then
    //    we go backwards and, for each dimension in which we need to move from the base cell, look
    //    for the first zero in the path. This position, for each dimension, will indicate the path
    //    suffix which has to be inverted to get the neighboring cell defined by cellPosition
    Nat firstInversePosition = NO_INVERSE_POSITION;
    if( inverseDimensionsNumericType > 0 ) {
        //If there is a need to do inverses, i.e. we are not adding the base cell itself
        Nat foundNumericTypeOfInverses = 0;
        Nat position = theNeighborCellWord.size() - 1u;
        while(true) {
            //Only consider the dimension that we need and look for the first opotrunity to invert the path suffix.
            Nat dimension = position % num_dimensions;
            //If we need to inverse in this dimension and this is the first found position in
            //this dimension from which on we should inverse then save the position index.
            if( cellPosition[ dimension ] && !theNeighborCellWord[ position ] &&
                ( invert_position[ dimension ] == NO_INVERSE_POSITION ) ) {
                invert_position[ dimension ] = position;
                //Since it will typically be the case the the binary word to the base cell will be
                //longer than the number of dimensions, we also find the first inverse positions.
                if( position < firstInversePosition ) {
                    firstInversePosition = position;
                }
                //Incremenet the number of fount inverses and check if this is all we need, if yes then break
                foundNumericTypeOfInverses++;
                if( foundNumericTypeOfInverses == inverseDimensionsNumericType ) {
                    break;
                }
            }
            if (position == 0) break;
            position--;
        }
    }

    //04. Since now all the inversion positions are found, we need to go through the path again and
    //    inverse it in the needed dimesnions starting from (corresponding) the found positions on.
    //    This will provide us with the path to the neighborind cell in the given dimension
    for( Nat index = firstInversePosition; index < theNeighborCellWord.size(); index++ ) {
        Nat dimension = index % num_dimensions;
        if( cellPosition[ dimension ] && ( index >= invert_position[ dimension ] ) ) {
            theNeighborCellWord[index] = !theNeighborCellWord[index];
        }
    }

    return GridCell( theGrid, theHeight, theNeighborCellWord );
}

Void GridOpenCell::cover_cell_and_borders( const GridCell& theCell, const GridTreePaving& theSet,
                                  BinaryWord& cellPosition, std::vector<GridOpenCell>& result ) {
    const Nat num_dimensions = theCell.grid().dimension();
    if( cellPosition.size() < num_dimensions ) {
        //Choose the left direction in the current dimension
        cellPosition.push_back( false );
        cover_cell_and_borders( theCell, theSet, cellPosition, result );
        cellPosition.pop_back( );
        //Choose the right direction in the current dimension
        cellPosition.push_back( true );
        cover_cell_and_borders( theCell, theSet, cellPosition, result );
        cellPosition.pop_back( );
    } else {
        //We have constructed the cell position relative to the base cell
        //for the case of _theGrid.dimension() dimensional space, now it
        //is time to compute this cell and add it to result
        GridCell neighborCell = GridOpenCell::neighboring_cell( theCell.grid(), theCell.height(), theCell.word(), cellPosition );
        //Check if the found neighboring cell is in theSet
        if( theSet.binary_tree()->is_enabled( neighborCell.word() ) ) {
            //So this cell is the enabled neighbor of theCell therefore we need to cover the boundary
            BinaryWord coverCellBaseWord = theCell.word();
            //Take the given word theCell.word() and then add the directions from the cellPosition
            //The latter start from the first axis till the last one, but the path in coverCellBaseWord
            //currently ends at some other axis so we need to align them when appending cellPosition
            for( Nat i = 0; i < num_dimensions ; i++ ) {
                coverCellBaseWord.push_back( cellPosition[ coverCellBaseWord.size() % num_dimensions ] );
            }
            //Add the resulting cover cell
            result.push_back( GridOpenCell( theCell.grid(), theCell.height(), coverCellBaseWord ) );
        }
    }
}

List<GridOpenCell> GridOpenCell::intersection( const GridOpenCell & theLeftOpenCell, const GridOpenCell & theRightOpenCell ) {
    List<GridOpenCell> result;

    //01 First check if one open cell is a subset of another open cell or if they overlap
    if( definitely( theLeftOpenCell.box().covers( theRightOpenCell.box() ) ) ) {
            //If theRightOpenCell is a subset of theLeftOpenCell
        result.push_back( theRightOpenCell );
    } else {
        if( definitely( theRightOpenCell.box().covers( theLeftOpenCell.box() ) ) ) {
            //If theLeftOpenCell is a subset of theRightOpenCell
            result.push_back( theLeftOpenCell );
        } else {
            if( definitely( theRightOpenCell.box().overlaps( theLeftOpenCell.box() ) ) ) {
                //02 If the open cells overlap then let us get the cells contained by the two open cells
                GridTreePaving theLeftCellSet = theLeftOpenCell.closure();
                GridTreePaving theRightCellSet = theRightOpenCell.closure();

                //03 Then we compute their intersection
                GridTreePaving intersectionSet = Ariadne::intersection( theLeftCellSet, theRightCellSet );
                //NOTE: It seems to me there is no need to recombine the resulting set, it will not reduce anything

                //04 Iterate through all the cells in the intersection and first add their interiors to the resulting set, second
                //check if the two cells have common border and/or vertex and if they do then add extra open cells, lying withing
                //the intersection and covering the common border and/or vertex.
                for ( GridTreeSubpaving::ConstIterator it = intersectionSet.begin(), end = intersectionSet.end(); it != end; it++) {
                    //Cover the interior of the cell and the borders with the cells bordered with
                    //the given one in each  positive direction in each dimension. The borders
                    //are covered only if the neighboring cell is also in the intersection set.
                    BinaryWord tmpWord;
                    cover_cell_and_borders( (*it), intersectionSet, tmpWord, result );
                }
            }
        }
    }

    return result;
}

//***************************************FRIENDS OF GridOpenCell*******************************************/

OutputStream& operator<<(OutputStream& os, const GridOpenCell& theGridOpenCell ) {
    //Write the grid data to the string stream
    return os << "GridOpenCell( " << theGridOpenCell.grid() <<
        ", Primary cell height: " << theGridOpenCell.height() <<
        ", Base cell path: " << theGridOpenCell.word() <<
        ", ExactBoxType (closure): " << theGridOpenCell.box() << " )";
}



//*******************************************GridTreeCursor***************************************/

GridTreeCursor::GridTreeCursor( ) :
    _currentStackIndex(-1), _pSubPaving(0), _theCurrentGridCell( Grid(), 0, BinaryWord() ) {
}

GridTreeCursor::GridTreeCursor(const GridTreeCursor & otherCursor) :
    _currentStackIndex(otherCursor._currentStackIndex), _theStack(otherCursor._theStack),
    _pSubPaving(otherCursor._pSubPaving), _theCurrentGridCell(otherCursor._theCurrentGridCell){
}

GridTreeCursor& GridTreeCursor::operator=(const GridTreeCursor & otherCursor) {
    _currentStackIndex = otherCursor._currentStackIndex;
    _theStack = otherCursor._theStack;
    _pSubPaving = otherCursor._pSubPaving;
    _theCurrentGridCell = otherCursor._theCurrentGridCell;

    return *this;
}

GridTreeCursor::GridTreeCursor(const GridTreeSubpaving * pSubPaving) :
    _currentStackIndex(-1), _pSubPaving(pSubPaving), _theCurrentGridCell( pSubPaving->root_cell() ) {
    //Remember that GridTreeSubset contains GridCell

    //Add the current node to the stack, since we are in it

    //IVAN S ZAPREEV:
    //NOTE: There is no need in allocating _theStack elements,
    //it is all done in the push method
    push(pSubPaving->_pRootTreeNode);
}

GridTreeCursor::~GridTreeCursor() {
    //IVAN S ZAPREEV:
    //WARNING: The subpaving should not be deallocated here!
    //There are no other data feels that we need to deallocate
    _pSubPaving = nullptr;
}

inline Void GridTreeCursor::push( BinaryTreeNode* pLatestNode ){
    //If we are out of free space, increase the Array's capacity
    if( static_cast<Nat>(_currentStackIndex +1) == ( _theStack.size()  ) ){
        _theStack.resize( _theStack.size() + STACK_SIZE_INCREMENT );
    }

    //The index for the tree node to be added
    _currentStackIndex += 1;

    //Put the tree node pointer into the stack
    _theStack[static_cast<Nat>(_currentStackIndex)] = pLatestNode;

}

inline BinaryTreeNode* GridTreeCursor::pop( ){
    BinaryTreeNode* pLastNode = nullptr;

    //If there are non-root nodes in the stack
    if( _currentStackIndex > 0 ){
        //Return the stack element at _currentStackIndex,
        //then decrement the current element's index (here inverted in order)
        _currentStackIndex--;
        return _theStack[ static_cast<Nat>(_currentStackIndex+1) ];
    }

    return pLastNode;
}

Bool GridTreeCursor::is_at_the_root() const {
    //If we are pointing at the zero cell of the Array then it
    //means that we are in the root of the tree
    return (_currentStackIndex == 0);
}

Bool GridTreeCursor::is_enabled() const {
    return _theStack[ static_cast<Nat>(_currentStackIndex) ]->is_enabled();
}

Bool GridTreeCursor::is_disabled() const {
    return _theStack[ static_cast<Nat>(_currentStackIndex) ]->is_disabled();
}

Bool GridTreeCursor::is_leaf() const {
    return _theStack[ static_cast<Nat>(_currentStackIndex) ]->is_leaf();
}

Bool GridTreeCursor::is_root() const {
    return  is_at_the_root();
}

Void GridTreeCursor::set_enabled() const {
    return _theStack[ static_cast<Nat>(_currentStackIndex) ]->set_enabled();
}

Void GridTreeCursor::set_disabled() const {
    return _theStack[ static_cast<Nat>(_currentStackIndex) ]->set_disabled();
}

Bool GridTreeCursor::is_left_child() const{
    //If there is a parent node and the given node is it's left child
    return ( _currentStackIndex > 0 ) && ( _theStack[ static_cast<Nat>(_currentStackIndex - 1) ]->left_node() == _theStack[ static_cast<Nat>(_currentStackIndex) ] );
}

Bool GridTreeCursor::is_right_child() const{
    //If there is a parent node and the given node is it's right child
    return ( _currentStackIndex > 0 ) && ( _theStack[ static_cast<Nat>(_currentStackIndex - 1) ]->right_node() == _theStack[ static_cast<Nat>(_currentStackIndex) ] );
}

GridTreeCursor& GridTreeCursor::move_up() {
    return move(BinaryTreeDirection::UP);
}

GridTreeCursor& GridTreeCursor::move_left() {
    return move(BinaryTreeDirection::LEFT);
}

GridTreeCursor& GridTreeCursor::move_right() {
    return move(BinaryTreeDirection::RIGHT);
}

inline Void GridTreeCursor::updateTheCurrentGridCell( BinaryTreeDirection up_or_left_or_right ) {
    if( up_or_left_or_right == BinaryTreeDirection::UP ){
        //We go "up" in the tree, so we remove the last bit of the path.
        _theCurrentGridCell._theWord.pop_back();
    } else {
        //Here up_or_left_or_right is either LEFT or RIGHT, so true defines moving to the right branch
        _theCurrentGridCell._theWord.push_back( up_or_left_or_right==BinaryTreeDirection::RIGHT );
    }
    _theCurrentGridCell = GridCell( _theCurrentGridCell._theGrid, _theCurrentGridCell._theHeight, _theCurrentGridCell._theWord );
}

inline GridTreeCursor& GridTreeCursor::move(BinaryTreeDirection up_or_left_or_right) {
    if( up_or_left_or_right == BinaryTreeDirection::UP ) {
        if( is_root() ) {
            throw NotAllowedMoveException(ARIADNE_PRETTY_FUNCTION);
        }
        //Remove the current node from the stack
        pop();
    }
    else {
        if( is_leaf() ){
            throw NotAllowedMoveException(ARIADNE_PRETTY_FUNCTION);
        }
        //If we are not in the leaf node then we can go down
        BinaryTreeNode* pNextNode;
        if( up_or_left_or_right == BinaryTreeDirection::RIGHT ) { //move to the right
            pNextNode = _theStack[ static_cast<Nat>(_currentStackIndex) ]->right_node();
        } else { // move to the left
            pNextNode = _theStack[ static_cast<Nat>(_currentStackIndex) ]->left_node();
        }
        //Put the node into the stack
        push(pNextNode);
    }

    //Recompute the PavingGridCell
    updateTheCurrentGridCell( up_or_left_or_right );
    //Return the object back
    return ( * this);
}

Bool GridTreeCursor::operator==(const GridTreeCursor& anotherGridTreeCursor) const {
    Bool areEqual = false;
    if( (this->_currentStackIndex >=0) && (anotherGridTreeCursor._currentStackIndex >=0 ) ){
        areEqual = (this->_theStack[static_cast<Nat>(this->_currentStackIndex)] == anotherGridTreeCursor._theStack[static_cast<Nat>(anotherGridTreeCursor._currentStackIndex)]);
    }
    return areEqual;
}

GridTreeSubpaving GridTreeCursor::operator*() {
    //IVAN S ZAPREEV:
    //NOTE: The first three parameters define the location of the _theStack[ _currentStackIndex ]
    //node with respect to the primary cell of the GridTreeSet.
    return GridTreeSubpaving( _theCurrentGridCell._theGrid,
                           _theCurrentGridCell._theHeight,
                           _theCurrentGridCell._theWord,
                           _theStack[ static_cast<Nat>(_currentStackIndex) ] );
}

const GridTreeSubpaving GridTreeCursor::operator*() const {
    return const_cast<GridTreeCursor&>(*this).operator*();
}

//************************************FRIENDS OF GridTreeCursor*****************************************/

OutputStream& operator<<(OutputStream& os, const GridTreeCursor& theGridTreeCursor) {
    const Int curr_stack_idx = theGridTreeCursor._currentStackIndex;
    os << "GridTreeCursor( " << theGridTreeCursor._pSubPaving <<
        ", Curr. stack index: " << curr_stack_idx  <<
        ", Stack data: [ ";
    for(Nat i = 0; i <= static_cast<Nat>(curr_stack_idx); i++){
        os << theGridTreeCursor._theStack[i]->node_to_string() << ( ( i < static_cast<Nat>(curr_stack_idx) ) ? "" : ", ");
    }
    return os<<" ], " << theGridTreeCursor._theCurrentGridCell << " )";
}

//***************************************GridTreeConstIterator************************************/

GridTreeConstIterator::GridTreeConstIterator(  ) : _pGridTreeCursor()  {
}

GridTreeConstIterator::GridTreeConstIterator( const GridTreeSubpaving * pSubPaving, const ValidatedKleenean firstLastNone )
    :  _pGridTreeCursor(pSubPaving)
{
    if( is_determinate( firstLastNone ) ) {
        //If the first/last enabled node is not found, it means that there are no elements
        //to iterate on, then we switch to the "end Iterator" state
        _is_in_end_state = ! navigate_to( definitely( firstLastNone ) );
    } else {
        //In this case if we do nothing, the cursor in the Iterator will point to the root node
        //Since this can be the only node in the tre we should add a marker that indicates that
        //the Iterator is at the end state.
        _is_in_end_state = true;
    }
}

Void GridTreeConstIterator::find_next_enabled_leaf() {
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

Bool GridTreeConstIterator::navigate_to(Bool firstLast){
    Bool isEnabledLeafFound = false;
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

Void GridTreeConstIterator::increment() {
    //If we are not done iterating
    if( ! _is_in_end_state){
        //We are at some enabled-leaf node and we want to find the next one.
        //The next node is somewhere on the right in the tree.
        find_next_enabled_leaf();
    }
}

Bool GridTreeConstIterator::equal( GridTreeConstIterator const & theOtherIterator) const {
    //Check if both iterators are in the "end Iterator" state
    Bool result = theOtherIterator._is_in_end_state && this->_is_in_end_state;

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

GridTreeConstIterator* GridTreeConstIterator::clone() const {
    return new GridTreeConstIterator(*this);
}

Bool GridTreeConstIterator::equals( ForwardConstantIteratorInterface<GridCell> const & theOtherIterator) const {
    GridTreeConstIterator const* theOtherIteratorPointer = dynamic_cast<GridTreeConstIterator const*>(&theOtherIterator);
    return theOtherIteratorPointer && (this->equal(*theOtherIteratorPointer));
}

Void GridTreeConstIterator::write(OutputStream& os) const {
    os << "GridTreeConstIterator(" << this->cursor() << ")";
}


//*******************************************GridTreeSubset******************************************/

SizeType GridTreeSubpaving::size() const {
    return BinaryTreeNode::count_enabled_leaf_nodes( this->binary_tree() );
}

Void GridTreeSubpaving::set_root_cell(Bool enabled_or_disabled)  {
    this->_pRootTreeNode->set(enabled_or_disabled);
}

Void GridTreeSubpaving::mince_to_tree_depth( const Nat theNewDepth ) {
    _pRootTreeNode->mince( theNewDepth );
}

Void GridTreeSubpaving::recombine() {
    _pRootTreeNode->recombine();
}

Nat GridTreeSubpaving::depth() const {
    return _pRootTreeNode->depth();
}

Bool GridTreeSubpaving::operator==(const GridTreeSubpaving& anotherGridTreeSubset) const {
    return ( this->_theGridCell == anotherGridTreeSubset._theGridCell ) &&
        ( ( * this->_pRootTreeNode ) == ( * anotherGridTreeSubset._pRootTreeNode ) );
}

GridTreeSubpaving::GridTreeSubpaving( const Grid& theGrid, const Nat theHeight,
                                       const BinaryWord& theWord, BinaryTreeNode * pRootTreeNode )
    : _pRootTreeNode(pRootTreeNode), _theGridCell(theGrid, theHeight, theWord) {
}

GridTreeSubpaving::GridTreeSubpaving( const GridTreeSubpaving &otherSubset )
    : _pRootTreeNode(otherSubset._pRootTreeNode),
      _theGridCell(otherSubset._theGridCell) {
}

GridTreeSubpaving::~GridTreeSubpaving() {
    //IVAN S ZAPREEV:
    //WARNING: This method should have no implementation what so ever
    //All the synamically allocatged data should be destroyed from the
    //corresponding Paving object
}

inline Nat GridTreeSubpaving::compute_number_subdiv( FloatDP theWidth, const FloatDP theMaxWidth) const{
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
    //NOTE: FloatDP now uses approximate operators by default
    Nat result = 0;
    if ( theWidth > theMaxWidth ){
        //result = (Nat) ceil( div(approx, log(approx, div(approx, theWidth, theMaxWidth ) ) , log(approx, R(2.0) ) ) );
        result = integer_cast<Nat>(ceil( div( near, log( near, div( near, theWidth, theMaxWidth ) ) , log( near, FloatDP(2.0) ) ) ) );
    }
    return result;
}

Bool GridTreeSubpaving::is_empty() const {
    return this->size() == 0;
}

DimensionType GridTreeSubpaving::dimension( ) const {
    return grid().dimension();
}


const Grid& GridTreeSubpaving::grid() const {
    return this->_theGridCell.grid();
}

GridCell GridTreeSubpaving::root_cell() const {
    return _theGridCell;
}

UpperBoxType GridTreeSubpaving::bounding_box() const {
    if(this->is_empty()) return ExactBoxType(this->dimension());

    GridTreePaving::ConstIterator iter=this->begin();
    UpperBoxType bbox = iter->box();

    for( ; iter!=this->end(); ++iter) {
        UpperBoxType cell = iter->box();
        for(Nat i = 0; i < cell.dimension(); ++i) {
            if(cell[i].lower().raw() < bbox[i].lower().raw()) bbox[i].set_lower(cell[i].lower());
            if(cell[i].upper().raw() > bbox[i].upper().raw()) bbox[i].set_upper(cell[i].upper());
        }
    }

    return bbox;

}

inline Nat GridTreeSubpaving::zero_cell_subdivisions_to_tree_subdivisions( const Nat numSubdivInDim, const Nat primaryCellHeight,
                                                                        const Nat primaryToRootCellPathLength ) const {
    //Here we take the height of the primary cell that the subpaving's root cell is rooted to
    //This height times the number of dimensions is the number of subdivisions to make in
    //the primary cell to reach the level of the zero cell. This plus the numSubdivInDim
    //times the number of dimensions gives us the total number of the primary cell we have
    //to make. But, the given paving is has the root node different from the primary cell,
    //thus we have to subtract the length of the path from the primary cell to the root cell
    //of this subpaving, to get the proper number of subdivisions to make in the binary tree
    Int theNewDepth = ( primaryCellHeight + numSubdivInDim ) * _theGridCell.grid().dimension() - primaryToRootCellPathLength;
    //If the new depth is not positive then we already have the required number
    //of subdivisions so then nothing has to be done, so we return zero!
    return (theNewDepth > 0) ? static_cast<Nat>(theNewDepth) : 0u;
}

Void GridTreeSubpaving::mince( Nat numSubdivInDim ) {
    mince_to_tree_depth( zero_cell_subdivisions_to_tree_subdivisions( numSubdivInDim, _theGridCell.height(), _theGridCell.word().size() ) );
}

GridTreeSubpaving::ConstIterator GridTreeSubpaving::begin() const {
    return GridTreeSubpaving::ConstIterator(this, true);
}

GridTreeSubpaving::ConstIterator GridTreeSubpaving::end() const {
    return GridTreeSubpaving::ConstIterator(this, indeterminate);
}

const BinaryTreeNode * GridTreeSubpaving::binary_tree() const {
    return _pRootTreeNode;
}

GridTreeSubpaving& GridTreeSubpaving::operator=( const GridTreeSubpaving &otherSubset) {
    _pRootTreeNode = otherSubset._pRootTreeNode;
    _theGridCell = otherSubset._theGridCell;

    return *this;
}

inline GridTreeSubpaving* GridTreeSubpaving::_branch(Bool left_or_right) const {
    return new GridTreeSubpaving(this->branch(left_or_right));
}

inline ForwardConstantIteratorInterface<GridCell>* GridTreeSubpaving::_begin() const {
    return new GridTreeConstIterator(this->begin());
}

inline ForwardConstantIteratorInterface<GridCell>* GridTreeSubpaving::_end() const {
    return new GridTreeConstIterator(this->end());
}

Bool GridTreeSubpaving::equals(const SubPavingInterface& paving) const {
    return (*this) == dynamic_cast<const GridTreeSubpaving&>(paving);
}

Bool GridTreeSubpaving::superset(const GridCell& theCell) const {
    return Ariadne::subset(theCell,*this);
}

Bool GridTreeSubpaving::subset(const SubPavingInterface& paving) const {
    return Ariadne::subset(*this, dynamic_cast<const GridTreeSubpaving&>(paving));
}

Bool GridTreeSubpaving::intersects(const SubPavingInterface& paving) const {
    return Ariadne::intersect(*this, dynamic_cast<const GridTreeSubpaving&>(paving));
}

GridTreeSubpaving* GridTreeSubpaving::clone( ) const {
    // Return a GridTreeSet to ensure that memory is copied.
    return new GridTreePaving(this->grid(),this->_pRootTreeNode);
}


GridTreeSubpaving GridTreeSubpaving::branch(Bool left_or_right) const {
    BinaryWord theWord=this->_theGridCell.word();
    theWord.append(left_or_right);
    BinaryTreeNode* pBinaryTreeNode=this->_pRootTreeNode->child_node(left_or_right);
    return GridTreeSubpaving(this->_theGridCell.grid(),this->_theGridCell.height(),theWord,pBinaryTreeNode);
}


Void GridTreeSubpaving::subdivide( FloatDP theMaxCellWidth ) {
    //1. Take the ExactBoxType of this GridTreeSubset's GridCell
    //   I.e. the box that corresponds to the root cell of
    //   the GridTreeSubset in the original space.
    const ExactBoxType& theRootCellBoxType = _theGridCell.box();

    //2. Compute the widths of the box in each dimension and the maximum number
    //   among the number of subdivisions that we need to do in each dimension
    //   in order to make the width in this dimension <= theMaxCellWidth.
    const DimensionType dimensions = _theGridCell.dimension();
    Nat max_num_subdiv_dim = 0, num_subdiv = 0, max_subdiv_dim = 0;

    for(Nat i = 0; i < dimensions; i++){
        //Get the number of required subdivisions in this dimension
        //IVAN S ZAPREEV:
        //NOTE: We compute sub_up because we do not want to have insufficient number of subdivisions
        num_subdiv = compute_number_subdiv( sub(up, theRootCellBoxType[i].upper().raw(), theRootCellBoxType[i].lower().raw() ) , theMaxCellWidth );

        //Compute the max number of subdivisions and the dimension where to do them
        if( num_subdiv >= max_num_subdiv_dim ){
            max_num_subdiv_dim = num_subdiv;
            max_subdiv_dim = i;
        }
    }

    //3. Let the maximum number of subdivisions M has to be done in dimension K  with the total number of
    //   dimensions N: 1 <= K <= N. This means that from this cell down we have to do M splits for dimension K.
    Nat needed_num_tree_subdiv = 0;
    //If we need to subdivide in one of the dimensions then
    if( max_num_subdiv_dim != 0 ){
        //3.1 Compute the dimension C for which we had the last split, we should start with the primary cell which is the root of
        //the GridTreeSet because from this cell we begin subdividing in dimension one by one: 1,2,...,N, then again 1,2,...,N.
        //The path to the root of the sub-paving is given by the binary word, its length gives the number of tree subdivisions:
        const Nat pathLength = _theGridCell.word().size();
        //If pathLength == 0 then there were no subdivisions in the tree, so we assign last_subdiv_dim == -1
        const Int last_subdiv_dim = ( pathLength == 0 ) ? -1 : static_cast<Int>(( pathLength - 1 )  % dimensions);

        //3.2 Compute the needed number of tree subdivisions by first computing how many subdivisions in the tree
        //we need to do to reach and split the dimension K one first time and then we should add the remaining
        //( M - 1 )*N tree subdevisions which will make shure that the K'th dimension is subdivided M times.
        Nat first_subdiv_steps;
        if( last_subdiv_dim == static_cast<Int>(max_subdiv_dim) ) {
            //If last_subdiv_dim == -1 then we will never get here
            first_subdiv_steps = dimensions; // C == K
        } else {
            //If last_subdiv_dim == -1 then we will add a needed extra subdivision
            first_subdiv_steps = static_cast<Nat>(static_cast<Int>(max_subdiv_dim) - last_subdiv_dim); // C < K
            if( last_subdiv_dim > static_cast<Int>(max_subdiv_dim) ) {
                //If last_subdiv_dim == -1 then we will never get here
                first_subdiv_steps = dimensions - first_subdiv_steps; // C > K
            }
        }
        needed_num_tree_subdiv = first_subdiv_steps + ( max_num_subdiv_dim - 1 ) * dimensions;
    }

    //Mince to the computed number of tree levels
    mince_to_tree_depth(needed_num_tree_subdiv);
}

FloatDPApproximation GridTreeSubpaving::measure() const {
    FloatDPApproximation result={0.0, dp};
    for(ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        result+=iter->box().measure().value();
    }
    return result;
}



GridTreeSubpaving::operator ListSet<ExactBoxType>() const {
    ListSet<ExactBoxType> result(this->root_cell().dimension());

    //IVAN S ZAPREEV:
    //NOTE: Push back the boxes, note that BoxListSet uses a vector, that
    //in its turn uses std:vector to store boxes in it (via push_back method),
    //the latter stores the copies of the boxes, not the references to them.
    for (GridTreeSubpaving::ConstIterator it = this->begin(), end = this->end(); it != end; it++ ) {
        result.push_back((*it).box());
    }

    return result;
}

Bool GridTreeSubpaving::superset( const ExactBoxType& theBoxType ) const {
    //Check that the box corresponding to the root node of the set
    //is not disjoint from theBoxType. If it is then the set is not a
    //subset of theBoxType otherwise we need to traverse the tree and check
    //if all it's enabled nodes give boxes that are subsets of theBoxType.

    ARIADNE_ASSERT( theBoxType.dimension() == root_cell().dimension() );

    Bool isASubSet = theBoxType.subset( root_cell().box() );
    if( ! isASubSet ) {
        //If the box is not covered by the cell corresponding to the root node
        //of the set, then clearly theBoxType is not a subset of this set.
        return false;
    } else {
        //Otherwise, is theBoxType is possibly a subset then we try to see furhter
        BinaryWord pathCopy( root_cell().word() );
        return definitely( GridTreeSubpaving::_superset( binary_tree(), grid(), root_cell().height(), pathCopy, theBoxType ) ) ;
    }
}

Bool GridTreeSubpaving::subset( const ExactBoxType& theBoxType ) const {
    //Check that the box corresponding to the root node of the set
    //is not disjoint from theBoxType. If it is then the set is not a
    //subset of theBoxType otherwise we need to traverse the tree and check
    //if all it's enabled nodes give boxes that are subsets of theBoxType.

    ARIADNE_ASSERT( theBoxType.dimension() == root_cell().dimension() );

    BinaryWord pathCopy( root_cell().word() );

    return definitely( GridTreeSubpaving::_subset( binary_tree(), grid(), root_cell().height(), pathCopy, theBoxType ) );
}

Bool GridTreeSubpaving::disjoint( const ExactBoxType& theBoxType ) const {
    //Check that the box corresponding to the root node of the set
    //is not disjoint from theBoxType. If it is then the set is not a
    //subset of theBoxType otherwise we need to traverse the tree and check
    //if all it's enabled nodes give boxes that are subsets of theBoxType.

    ARIADNE_ASSERT( theBoxType.dimension() == root_cell().dimension() );

    BinaryWord pathCopy( root_cell().word() );

    return definitely( GridTreeSubpaving::_disjoint( binary_tree(), grid(), root_cell().height(), pathCopy, theBoxType ) );
}

Bool GridTreeSubpaving::intersects( const ExactBoxType& theBoxType ) const {
    //Check that the box corresponding to the root node of the set
    //is not disjoint from theBoxType. If it is then the set is not a
    //subset of theBoxType otherwise we need to traverse the tree and check
    //if all it's enabled nodes give boxes that are subsets of theBoxType.

    ARIADNE_ASSERT( theBoxType.dimension() == root_cell().dimension() );

    BinaryWord pathCopy( root_cell().word() );

    return definitely( GridTreeSubpaving::_intersects( binary_tree(), grid(), root_cell().height(), pathCopy, theBoxType ) );
}

ValidatedLowerKleenean GridTreeSubpaving::covers( const ExactBoxType& theBoxType ) const {
    //Simply check if theBoxType is covered by the set and then make sure that
    //all tree cells that are not disjoint from theBoxType are enabled

    ARIADNE_ASSERT( theBoxType.dimension() == root_cell().dimension() );

    Bool isASubSet = theBoxType.subset( root_cell().box() );
    if( ! isASubSet ) {
        //If the box is not covered by the cell corresponding to the root node
        //of the set, then clearly theBoxType is not a subset of this set.
        return false;
    } else {
        //Otherwise, is theBoxType is possibly a subset then we try to see furhter
        BinaryWord pathCopy( root_cell().word() );
        return GridTreeSubpaving::_superset( binary_tree(), grid(), root_cell().height(), pathCopy, cast_exact_box(widen(theBoxType)) );
    }
}

ValidatedLowerKleenean GridTreeSubpaving::inside( const ExactBoxType& theBoxType ) const {
    //Check that the box corresponding to the root node of the set
    //is not disjoint from theBoxType. If it is then the set is not a
    //subset of theBoxType otherwise we need to traverse the tree and check
    //if all it's enabled nodes give boxes that are subsets of theBoxType.

    ARIADNE_ASSERT( theBoxType.dimension() == root_cell().dimension() );

    BinaryWord pathCopy( root_cell().word() );

    return GridTreeSubpaving::_subset( binary_tree(), grid(), root_cell().height(), pathCopy, cast_exact_box(narrow(theBoxType)) );
}

ValidatedLowerKleenean GridTreeSubpaving::separated( const ExactBoxType& theBoxType ) const {
    //Simply check if the box does not intersect with the set

    ARIADNE_ASSERT( theBoxType.dimension() == root_cell().dimension() );

    BinaryWord pathCopy( root_cell().word() );

    return GridTreeSubpaving::_disjoint( binary_tree(), grid(), root_cell().height(), pathCopy, cast_exact_box(widen(theBoxType)) );
}

ValidatedLowerKleenean GridTreeSubpaving::overlaps( const ExactBoxType& theBoxType ) const {
    //Check if the box of the root cell overlaps with theBoxType,
    //if not then theBoxType does not intersect with the cell,
    //otherwise we need to find at least one enabled node
    //in the binary tree, such that it's box overlaps theBoxType.

    ARIADNE_ASSERT( theBoxType.dimension() == root_cell().dimension() );

    BinaryWord pathCopy( root_cell().word() );

    return GridTreeSubpaving::_intersects( binary_tree(), grid(), root_cell().height(), pathCopy, cast_exact_box(narrow(theBoxType)) );
}

ValidatedKleenean GridTreeSubpaving::_superset( const BinaryTreeNode* pCurrentNode, const Grid& theGrid,
                                             const Nat theHeight, BinaryWord &theWord, const ExactBoxType& theBoxType ) {
    ValidatedKleenean result;

    //Check if the current node's cell intersects with theBoxType
    ExactBoxType theCellsBoxType = GridCell::compute_box( theGrid, theHeight, theWord );
    ValidatedKleenean doIntersect = theCellsBoxType.overlaps( theBoxType );

    if( ! definitely(doIntersect) ) {
        //If theBoxType does not intersect with the cell then for the covering relation
        //is is not important if we add or remove this cell, so we return true
        result = true;
    } else {
        //If the cell possibly intersects with theBoxType then
        if( pCurrentNode->is_leaf() ) {
            if( pCurrentNode->is_enabled() ){
                //If this is an enabled node that possibly or definitely intersects
                //with theBoxType then the covering property is fine, so we return true
                result = true;
            } else {
                //If the node is disabled then if it definitely intersects with theBoxType
                //we have to report false, as there is not covering, in case it is only
                //possible intersecting with theBoxType then we return possibly. The latter
                //is because we are not completely sure, if the intersection does not
                //have place then the covering property is not broken.
                result = ! doIntersect;
            }
        } else {
            //The node is not a leaf so we need to go down and see if the cell
            //falls into sub cells for which we can sort things out
            theWord.push_back(false);
            const ValidatedKleenean result_left = GridTreeSubpaving::_superset( pCurrentNode->left_node(), theGrid, theHeight, theWord, theBoxType );
            theWord.pop_back();

            if( definitely(not result_left) ) {
                //If there is definitely no covering property then this is all
                //we need to know so there is not need to check the other branch.
                result = false;
            } else {
                //If the covering property holds or is possible, then we still
                //need to check the second branch because it can change the outcome.
                theWord.push_back(true);
                const ValidatedKleenean result_right = GridTreeSubpaving::_superset( pCurrentNode->right_node(), theGrid, theHeight, theWord, theBoxType );
                theWord.pop_back();

                if( definitely(not result_right) ) {
                    //IF: The right sub-node reports false, then the result is false
                    //NOTE: We already sorted out the case of ( (!result_left) == true) before
                    result = false;
                } else {
                    if( definitely(result_left && result_right) ) {
                        //ELSE: If both sub-nodes report true then it is true
                        result = true;
                    } else {
                        if( is_indeterminate(result_left) || is_indeterminate(result_right) ) {
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

ValidatedKleenean GridTreeSubpaving::_subset( const BinaryTreeNode* pCurrentNode, const Grid& theGrid,
                                           const Nat theHeight, BinaryWord &theWord, const ExactBoxType& theBoxType ) {
    ValidatedKleenean result;

    //Check if the current node overlaps with theBoxType
    ExactBoxType theCellsBoxType = GridCell::compute_box( theGrid, theHeight, theWord );
    ValidatedKleenean isASubset = theCellsBoxType.subset( theBoxType );

    if( definitely(isASubset) ){
        //It does not matter if pCurrentNode has enableds leaves or not we already know that the cell
        //corresponding to this node is geometrically a subset of theBoxType.
        result = true;
    } else {
        //If the cell corresponding to pCurrentNode is not a subset, then
        if( pCurrentNode->is_leaf() && definitely(not isASubset) ) {
            //If pCurrentNode is a leaf node and geometrically the cell (corresponding to the node theCellsBoxType) is not a
            //subset of theBoxType, then: if it is enabled then pCurrentNode is not a subset of theBoxType but otherwise it is.
            result = ! pCurrentNode->is_enabled();
        } else {
            if( pCurrentNode->is_leaf() && is_indeterminate( isASubset ) ) {
                //If we are in a leaf node but we do not know for sure if the given cell
                //is a subset of theBoxType then we can only check if it is enabled or not
                if( pCurrentNode->is_enabled() ){
                    //For an enabled-leaf node (a filled cell) we do not know if it is a subset of theBoxType
                    result = indeterminate;
                } else {
                    //The node is disabled, so it represents an empty set, which is a subset of any set
                    result = true;
                }
            } else {
                //The node is not a leaf, and we either know that the cell of pCurrentNode is not a geometrical subset
                //of theBoxType or we are not sure that it is, This means that we can do recursion to sort things out.
                theWord.push_back(false);
                const ValidatedKleenean result_left = GridTreeSubpaving::_subset( pCurrentNode->left_node(), theGrid, theHeight, theWord, theBoxType );
                theWord.pop_back();

                if( definitely(not result_left) ) {
                    //If the left branch is not a subset, then there is no need to check the right one
                    result = false;
                } else {
                    //if we still do not know the answer, then we check the right branch
                    theWord.push_back(true);
                    const ValidatedKleenean result_right = GridTreeSubpaving::_subset( pCurrentNode->right_node(), theGrid, theHeight, theWord, theBoxType );
                    theWord.pop_back();

                    if( definitely(not result_right) ) {
                        //IF: The right sub-node reports false, then the result is false
                        //NOTE: We already sorted out the case of ( (!result_left) == true) before
                        result = false;
                    } else {
                        if( definitely(result_left && result_right) ) {
                            //ELSE: If both sub-nodes report true then it is true
                            result = true;
                        } else {
                            if( is_indeterminate(result_left) || is_indeterminate(result_right) ) {
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

ValidatedKleenean GridTreeSubpaving::_disjoint( const BinaryTreeNode* pCurrentNode, const Grid& theGrid,
                                             const Nat theHeight, BinaryWord &theWord, const ExactBoxType& theBoxType ) {
    ValidatedKleenean intersect;

    //Check if the current node overlaps with theBoxType
    ExactBoxType theCellsBoxType = GridCell::compute_box( theGrid, theHeight, theWord );
    ValidatedKleenean doPossiblyIntersect = !theCellsBoxType.disjoint( theBoxType );

    if( possibly(doPossiblyIntersect) ) {
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
            const ValidatedKleenean intersect_left = GridTreeSubpaving::_intersects( pCurrentNode->left_node(), theGrid, theHeight, theWord, theBoxType );
            theWord.pop_back();

            //
            //WARNING: I know how to write a shorter code, like:
            //          if( ! ( intersect = definitely( intersect_left ) ) ) {
            //              ...
            //          }
            // I DO NOT DO THIS ON PERPOSE, KEEP THE CODE EASILY UNDERSTANDABLE!
            //
            if( definitely(intersect_left) ) {
                //If we definitely have intersection for the left branch then answer is true
                intersect = true;
            } else {
                //If we still not sure/ or do not know then try to search further, i.e. check the right node
                theWord.push_back(true);
                const ValidatedKleenean intersect_right = GridTreeSubpaving::_intersects( pCurrentNode->right_node(), theGrid, theHeight, theWord, theBoxType );
                theWord.pop_back();
                if( definitely(intersect_right) ) {
                    //If we definitely have intersection for the right branch then answer is true
                    intersect = true;
                } else {
                    //Now either we have indeterminate answers or we do not know, if one
                    //if the answers for one of the branches was indeterminate then we
                    //report indeterminate, otherwise it is definitely false.
                    if( is_indeterminate( intersect_left ) || is_indeterminate( intersect_right )  ) {
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

ValidatedKleenean GridTreeSubpaving::_intersects( const BinaryTreeNode* pCurrentNode, const Grid& theGrid,
                                               const Nat theHeight, BinaryWord &theWord, const ExactBoxType& theBoxType ) {
    ValidatedKleenean result;

    //Check if the current node overlaps with theBoxType
    ExactBoxType theCellsBoxType = GridCell::compute_box( theGrid, theHeight, theWord );
    ValidatedKleenean doPossiblyIntersect = theCellsBoxType.overlaps( theBoxType );

    if( possibly(doPossiblyIntersect) ) {
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
            const ValidatedKleenean result_left = GridTreeSubpaving::_intersects( pCurrentNode->left_node(), theGrid, theHeight, theWord, theBoxType );
            theWord.pop_back();

            //
            //WARNING: I know how to write a shorter code, like:
            //          if( ! ( result = definitely( result_left ) ) ) {
            //              ...
            //          }
            // I DO NOT DO THIS ON PERPOSE, KEEP THE CODE EASILY UNDERSTANDABLE!
            //
            if( definitely(result_left) ) {
                //If we definitely have intersection for the left branch then answer is true
                result = true;
            } else {
                //If we still not sure/ or do not know then try to search further, i.e. check the right node
                theWord.push_back(true);
                const ValidatedKleenean result_right = GridTreeSubpaving::_intersects( pCurrentNode->right_node(), theGrid, theHeight, theWord, theBoxType );
                theWord.pop_back();
                if( definitely(result_right) ) {
                    //If we definitely have intersection for the right branch then answer is true
                    result = true;
                } else {
                    //Now either we have indeterminate answers or we do not know, if one
                    //if the answers for one of the branches was indeterminate then we
                    //report indeterminate, otherwise it is definitely false.
                    if( is_indeterminate( result_left ) || is_indeterminate( result_right )  ) {
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

Void GridTreeSubpaving::draw(CanvasInterface& theGraphic, const Projection2d& theProjection) const {
    for(GridTreeSubpaving::ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        iter->box().draw(theGraphic,theProjection);
    }
}

OutputStream&
GridTreeSubpaving::write(OutputStream& os) const
{
    return os << (*this);
}

//************************************FRIENDS OF GridTreeSubset*****************************************/

Bool subset( const GridCell& theCell, const GridTreeSubpaving& theSet ) {
    Bool result = false;

    //Test that the Grids are equal
    ARIADNE_ASSERT( theCell.grid() == theSet.grid() );

    //Test if theCell is a subset of theSet, first check
    //if the cell can be a subset of the given tree.
    BinaryWord pathPrefixCell, pathPrefixSet;
    if( subset( theCell, theSet.root_cell(), &pathPrefixCell, &pathPrefixSet, nullptr) ) {
        //It can and thus pathPrefixSet is a prefix of pathPrefixCell. Both
        //paths start in the same primary cell. Also note that pathPrefixSet
        //is a path from primary_cell_height to the root node of the binary
        //tree in theSet. Therefore, removing 0..pathPrefixSet.size() elements
        //from pathPrefixCell will give us a path to theCell in the tree of theSet.
        pathPrefixCell.erase( pathPrefixCell.begin(), pathPrefixCell.begin() + static_cast<long>(pathPrefixSet.size()) );

        //Check that the cell given by pathPrefixCell is enabled in the tree of theSet
        result = theSet.binary_tree()->is_enabled( pathPrefixCell );
    } else {
        //DO NOTHING: the cell is a strict superset of the tree
    }
    return result;
}

Bool intersect( const GridCell& theCell, const GridTreeSubpaving& theSet ) {
    Bool result = false;

    //Test that the Grids are equal
    ARIADNE_ASSERT( theCell.grid() == theSet.grid() );

    //If the primary cell of the theCell is lower that that of theSet
    //Then we re-root theCell to the primary cell theSet.root_cell().height()
    const GridCell * pWorkGridCell;
    const Nat theSetsPCellHeight = theSet.root_cell().height();
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
    pathFromPCellCellToSetsRootNode.append( theSet.root_cell().word() );

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
            for( Nat i = pathFromPCellCellToSetsRootNode.size(); i < workCellWord.size(); i++ ) {
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

// This is a helper method, it receives two GridTreeSubset elements
// then it computes the primary cell that is common to them in a sence that
// these sets can be rooted to it. After that the method updates the paths
// with the information about the paths from the found primary cell to the
// root binary tree nodes of both sets.
static Void common_primary_cell_path(const GridTreeSubpaving& theSet1, const GridTreeSubpaving& theSet2, BinaryWord &pathCommonPCtoRC1, BinaryWord &pathCommonPCtoRC2 ) {
    //Get the root cells for the subsets
    GridCell rootCell1 = theSet1.root_cell();
    GridCell rootCell2 = theSet2.root_cell();
    //Get the heights of the primary cells for both subsets
    const Nat heightPC1 = rootCell1.height();
    const Nat heightPC2 = rootCell2.height();
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

// This method locates the node in the tree rooted to pSuperTreeRootNode that corresponds to
//  the path pathFromSuperToSub. In case we encounter a leaf node then we stop and return the node
static const BinaryTreeNode * locate_node( const BinaryTreeNode * pSuperTreeRootNode, const BinaryWord& pathFromSuperToSub ){
    //Locate the node in the pSuperTreeRootNode such that it corresponds to pSubTreeRootNode
    const BinaryTreeNode * pCurrentSuperTreeNode = pSuperTreeRootNode;
    for( Nat i = 0; i < pathFromSuperToSub.size(); i++ ) {
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

// \brief This is a helper functiuon for Bool subset( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 )
// pSuperTreeRootNode is the root tree node, pathFromSuperToSub is the path from this root node to the root node
// pSubTreeRootNode. In general we see if the set represented by pSubTreeRootNode is a subset of pSuperTreeRootNode.
// If pSubTreeRootNode has no enabled leaf nodes then the result is always true, else if pSuperTreeRootNode has no
// enabled leaf nodes then the result is always false.
static Bool subset(const BinaryTreeNode * pSubTreeRootNode, const BinaryTreeNode * pSuperTreeRootNode, const BinaryWord & pathFromSuperToSub) {
    Bool result = false;

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

// This is a helper functiuon for Bool subset( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 )
// pSuperTreeRootNode is the root tree node, pathFromSuperToSub is the path from this root node to the root node
// pSubTreeRootNode. In general we see if the set represented by pSuperTreeRootNode is a subset of pSubTreeRootNode.
// If pSuperTreeRootNode has no enabled leaf nodes then the result is always true, else if pSubTreeRootNode has no
// enabled leaf nodes then the result is always false.
static Bool subset( const BinaryTreeNode * pSuperTreeRootNode, const BinaryWord & pathFromSuperToSub, const BinaryTreeNode * pSubTreeRootNode ) {
    Bool result = false;

    //First we iterate throught the path pathFromSuperToSub trying to reach the common node with pSubTreeRootNode
    //Since we want to know if pSuperTreeRootNode is a subset of pSubTreeRootNode, the branches of pSuperTreeRootNode
    //that we ommit traveling the path pathFromSuperToSub, should contain no enabled leaf nodes. Because otherwise
    //pSuperTreeRootNode is not a subset of pSubTreeRootNode. Apart from this, once we encounter a leaf node on the
    //path we stop because the we can already decide if pSuperTreeRootNode is a subset of pSubTreeRootNode.
    Nat path_element = 0;
    Bool areExtraLeavesDisabled = true;
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
                    //     Bool subset(const BinaryTreeNode *, const BinaryTreeNode * )
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

Bool subset( const GridTreeSubpaving& theSet1, const GridTreeSubpaving& theSet2 ) {
    Bool result = false;

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

// This is a helper functiuon for Bool overlap( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 )
// pSuperTreeRootNode is the root tree node, pathFromSuperToSub is the path from this root node to the root node
// pSubTreeRootNode. In general we see if the sets represented by pSuperTreeRootNode and pSubTreeRootNode overlap.
// If at least one of pSuperTreeRootNode and pSubTreeRootNode have no enabled leaf nodes then the result is always false.
static Bool intersect(const BinaryTreeNode * pSuperTreeRootNode, const BinaryWord & pathFromSuperToSub, const BinaryTreeNode * pSubTreeRootNode) {
    Bool result = false;

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
            result = BinaryTreeNode::intersect( pCurrentSuperTreeNode, pSubTreeRootNode );
        }
    }

    return result;
}

Bool intersect( const GridTreeSubpaving& theSet1, const GridTreeSubpaving& theSet2 ) {
    Bool result = false;

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
        result = intersect( theSet1.binary_tree(), pathCommonPCtoRC2, theSet2.binary_tree() );
    } else {
        if( pathCommonPCtoRC2.is_prefix( pathCommonPCtoRC1 ) ){
            //theSet1 is located somewhere within the bounding box of theSet2
            pathCommonPCtoRC1.erase_prefix( pathCommonPCtoRC2.size() );
            result = intersect( theSet2.binary_tree(), pathCommonPCtoRC1, theSet1.binary_tree() );
        } else {
            //The sets do not overlap
            result = false;
        }
    }

    return result;
}

OutputStream& operator<<(OutputStream& os, const GridTreeSubpaving& theGridTreeSubset) {
    return os << "GridTreeSubset( Primary cell: " << theGridTreeSubset.root_cell() << ", " << (*theGridTreeSubset.binary_tree()) <<" )";
}


//********************************************GridTreeSet*********************************************/

GridCell GridTreePaving::smallest_enclosing_primary_cell( const UpperBoxType& theBoxType ) const {
    ARIADNE_ASSERT_MSG( this->dimension() == theBoxType.dimension(), "Cannot find enclosing cell for ExactBoxType  " << theBoxType << " for GridTreeSet with grid " << this->grid() );

    return GridCell::smallest_enclosing_primary_cell( theBoxType, this->grid() );
}

Void GridTreePaving::adjoin( const GridCell& theCell ) {
    ARIADNE_ASSERT_MSG( this->grid() == theCell.grid(), "Cannot adjoin GridCell with grid "<<theCell.grid()<<" to GridTreeSet with grid "<<this->grid() );
    Bool has_stopped = false;
    //Align the paving and the cell
    BinaryTreeNode* pBinaryTreeNode = align_with_cell( theCell.height(), true, false, has_stopped );

    //If we are not trying to adjoin something into an enabled sub cell of the paving
    if( ! has_stopped ){
        //Add the enabled cell to the binary tree.
        pBinaryTreeNode->add_enabled( theCell.word() );
    }
}

Void GridTreePaving::adjoin( const GridTreeSubpaving& theOtherSubPaving ) {
    ARIADNE_ASSERT_MSG( this->grid() == theOtherSubPaving.root_cell().grid(), "Cannot adjoin GridTreeSubset with grid "<<theOtherSubPaving.root_cell().grid()<<" to GridTreeSet with grid "<<this->grid() );

    Bool has_stopped = false;
    //Align the paving and the cell
    BinaryTreeNode* pBinaryTreeNode = align_with_cell( theOtherSubPaving.root_cell().height(), true, false, has_stopped );

    //If we are not trying to adjoin something into an enabled sub cell of the paving
    if( ! has_stopped ){
        //Now, the pBinaryTreeNode of this paving corresponds to the primary cell common with theOtherSubPaving.
        //The theOtherSubPaving's root node is defined the path theOtherSubPaving.word() which starts in the
        //node corresponding to the common primary cell.
        pBinaryTreeNode->add_enabled( theOtherSubPaving.binary_tree(), theOtherSubPaving.root_cell().word() );
    }
}

Void GridTreePaving::adjoin( const SubPavingInterface& theOtherSubPaving ) {
    this->adjoin(dynamic_cast<const GridTreeSubpaving&>(theOtherSubPaving));
}

Void GridTreePaving::restrict( const SubPavingInterface& theOtherSubPaving ) {
    this->restrict(dynamic_cast<const GridTreeSubpaving&>(theOtherSubPaving));
}

Void GridTreePaving::remove( const SubPavingInterface& theOtherSubPaving ) {
    this->remove(dynamic_cast<const GridTreeSubpaving&>(theOtherSubPaving));
}


GridTreePaving::GridTreePaving( ) : GridTreeSubpaving( Grid(), 0, BinaryWord(), new BinaryTreeNode( false ) ){
}

GridTreePaving::GridTreePaving( const Grid& theGrid, const Bool enable  ) :
    GridTreeSubpaving( theGrid, 0, BinaryWord(), new BinaryTreeNode( enable ) ){
}

GridTreePaving::GridTreePaving( const Grid& theGrid, const Nat theHeight, BinaryTreeNode * pRootTreeNode ) :
    GridTreeSubpaving( theGrid, theHeight, BinaryWord(), pRootTreeNode ){
}

GridTreePaving::GridTreePaving( const GridCell& theGridCell  ) :
    GridTreeSubpaving( theGridCell.grid(), theGridCell.height(), BinaryWord(), new BinaryTreeNode( false ) ){
    this->adjoin(theGridCell);
}

GridTreePaving::GridTreePaving( const Nat theDimension, const Bool enable ) :
    GridTreeSubpaving( Grid( theDimension, FloatDP(1.0) ), 0, BinaryWord(), new BinaryTreeNode( enable )) {
    //We want a [0,1]x...[0,1] cell in N dimensional space with no sxaling or shift of coordinates:
    //1. Create a new non scaling grid with no shift of the coordinates
    //2. The height of the primary cell is zero, since is is [0,1]x...[0,1] itself
    //3. The binary word that describes the path from the primary cell to the root
    //   of the tree is empty, because any paving always has a primary cell as a root
    //4. A new disabled binary tree node, gives us the root for the paving tree
}

GridTreePaving::GridTreePaving(const Grid& theGrid, const ExactBoxType & theLatticeBoxType ) :
    GridTreeSubpaving( theGrid, GridCell::smallest_enclosing_primary_cell_height( theLatticeBoxType ),
                    BinaryWord(), new BinaryTreeNode( false ) ) {
    //1. The main point here is that we have to compute the smallest primary cell that contains theBoundingBoxType
    //2. This cell is defined by it's height and becomes the root of the GridTreeSet
    //3. ExactPoint 2. implies that the word to the root of GridTreeSubset should be set to
    //   empty and we have only one disabled node in the binary tree
}

GridTreePaving::GridTreePaving( const Grid& theGrid, Nat theHeight, const BooleanArray& theTree, const BooleanArray& theEnabledCells ) :
    GridTreeSubpaving( theGrid, theHeight, BinaryWord(), new BinaryTreeNode( theTree, theEnabledCells ) ) {
    //Use the super class constructor and the binary tree constructed from the arrays: theTree and theEnabledCells
}

GridTreePaving::GridTreePaving( const GridTreePaving & theGridTreeSet ) :
    GridTreeSubpaving( theGridTreeSet._theGridCell.grid(), theGridTreeSet._theGridCell.height(),
                    theGridTreeSet._theGridCell.word(), new BinaryTreeNode( *theGridTreeSet._pRootTreeNode )) {
    //Call the super constructor: Create an exact copy of the tree, copy the bounding box
}

GridTreePaving& GridTreePaving::operator=( const GridTreePaving & theGridTreeSet ) {
    //Delete the old tree and make a new one
    if(this!=&theGridTreeSet) {
        if( GridTreeSubpaving::_pRootTreeNode != nullptr){
            delete GridTreeSubpaving::_pRootTreeNode;
            GridTreeSubpaving::_pRootTreeNode = nullptr;
        }
        static_cast<GridTreeSubpaving&>(*this) =
            GridTreeSubpaving( theGridTreeSet._theGridCell.grid(), theGridTreeSet._theGridCell.height(),
                            theGridTreeSet._theGridCell.word(), new BinaryTreeNode( *theGridTreeSet._pRootTreeNode ));
    }
    return *this;
}

GridTreePaving* GridTreePaving::clone() const {
    return new GridTreePaving( *this );
}

GridTreePaving::~GridTreePaving() {
    if( GridTreeSubpaving::_pRootTreeNode != nullptr){
        delete GridTreeSubpaving::_pRootTreeNode;
        GridTreeSubpaving::_pRootTreeNode = nullptr;
    }
}

Void GridTreePaving::up_to_primary_cell( const Nat toPCellHeight ){
    const Nat fromPCellHeight = this->root_cell().height();

    //The primary cell of this paving is lower then the one in the other paving so this
    //paving's has to be rerooted to another primary cell and then we merge the pavings.
    //1. Compute the path
    BinaryWord primaryCellPath = GridCell::primary_cell_path( this->root_cell().grid().dimension(), toPCellHeight, fromPCellHeight );
    //2. Substitute the root node of the paiving with the extended tree
    this->_pRootTreeNode = BinaryTreeNode::prepend_tree( primaryCellPath, this->_pRootTreeNode );
    //3. Update the GridCell that corresponds to the root of this GridTreeSubset
    this->_theGridCell = GridCell( this->_theGridCell.grid(), toPCellHeight, BinaryWord() );
}

BinaryTreeNode* GridTreePaving::align_with_cell( const Nat otherPavingPCellHeight, const Bool stop_on_enabled,
                                              const Bool stop_on_disabled, Bool & has_stopped ) {
    const Nat thisPavingPCellHeight = this->root_cell().height();

    //The current root node of the GridTreeSet
    BinaryTreeNode * pBinaryTreeNode = this->_pRootTreeNode;

    if( thisPavingPCellHeight > otherPavingPCellHeight ){

        //The primary cell of this paving is higher then the one of the other paving.
        //1. We locate the path to the primary cell node common with the other paving
        BinaryWord primaryCellPath = GridCell::primary_cell_path( this->root_cell().grid().dimension(), thisPavingPCellHeight, otherPavingPCellHeight );

        //2. Locate the binary tree node corresponding to this primnary cell
        Nat position = 0;
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

Void GridTreePaving::_adjoin_outer_approximation( const Grid & theGrid, BinaryTreeNode * pBinaryTreeNode, const Nat primary_cell_height,
                                               const Nat max_mince_depth,  const CompactSetInterface& theSet, BinaryWord * pPath ){
    //Compute the cell correspomding to the current node
    GridCell theCurrentCell( theGrid, primary_cell_height, *pPath );

    const OpenSetInterface* pOpenSet=dynamic_cast<const OpenSetInterface*>(static_cast<const SetInterfaceBase*>(&theSet));

    if( definitely( theSet.separated( theCurrentCell.box() ) ) ) {
        //DO NOTHING: We are in the node whoes representation in the original space is
        //disjoint from pSet and thus there will be nothing added to this cell.
    } else if( pOpenSet && definitely( pOpenSet->covers( theCurrentCell.box() ) ) ) {
        pBinaryTreeNode->make_leaf(true);
    } else {
        //This node's cell is not disjoint from theSet, thus it is possible to adjoin elements
        //of its outer approximation, unless this node is already and enabled leaf node.
        if( pBinaryTreeNode->is_enabled() ){ //NOTE: A non-leaf node can not be enabled so this check suffices
            //DO NOTHING: If it is enabled, then we can not add anything new to it.
        } else {
            //If the node is not enabled, so may be we can add something from the outer approximation of theSet.
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
    //Return to the previous level, since the initial evaluate is made
    //with the empty word, we check that it is not yet empty.
    if( pPath->size() > 0 ) {
        pPath->pop_back();
    }
}

Void GridTreePaving::_adjoin_outer_approximation( const Grid & theGrid, BinaryTreeNode * pBinaryTreeNode, const Nat primary_cell_height,
                                               const Nat max_mince_depth,  const ValidatedCompactSetInterface& theSet, BinaryWord * pPath ){
    //Compute the cell correspomding to the current node
    GridCell theCurrentCell( theGrid, primary_cell_height, *pPath );

    const OpenSetInterface* pOpenSet=dynamic_cast<const OpenSetInterface*>(static_cast<const SetInterfaceBase*>(&theSet));

    if( definitely( theSet.separated( theCurrentCell.box() ) ) ) {
        //DO NOTHING: We are in the node whoes representation in the original space is
        //disjoint from pSet and thus there will be nothing added to this cell.
    } else if( pOpenSet && definitely( pOpenSet->covers( theCurrentCell.box() ) ) ) {
        pBinaryTreeNode->make_leaf(true);
    } else {
        //This node's cell is not disjoint from theSet, thus it is possible to adjoin elements
        //of its outer approximation, unless this node is already and enabled leaf node.
        if( pBinaryTreeNode->is_enabled() ){ //NOTE: A non-leaf node can not be enabled so this check suffices
            //DO NOTHING: If it is enabled, then we can not add anything new to it.
        } else {
            //If the node is not enabled, so may be we can add something from the outer approximation of theSet.
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
    //Return to the previous level, since the initial evaluate is made
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
Void GridTreePaving::_adjoin_lower_approximation( const Grid & theGrid, BinaryTreeNode * pBinaryTreeNode, const Nat primary_cell_height,
                                               const Nat max_mince_depth,  const OvertSetInterface& theSet, BinaryWord * pPath ){
    //Compute the cell correspomding to the current node
    GridCell theCurrentCell( theGrid, primary_cell_height, *pPath );

    if( definitely( theSet.overlaps( theCurrentCell.box() ) ) ) {
        if( pPath->size() >= max_mince_depth ) {
            //We should not mince any further.
            //If the cell is not a leaf, then some subset is enabled,
            //so the lower approximation does not add any information.
            //If the cell is a leaf, we mark it as enabled.
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
    //Return to the previous level, since the initial evaluate is made
    //with the empty word, we check that it is not yet empty.
    if( pPath->size() > 0 ) {
        pPath->pop_back();
    }
}

Void GridTreePaving::_adjoin_lower_approximation( const Grid & theGrid, BinaryTreeNode * pBinaryTreeNode, const Nat primary_cell_height,
                                               const Nat max_mince_depth,  const OpenSetInterface& theSet, BinaryWord * pPath ){
    //Compute the cell corresponding to the current node
    GridCell theCurrentCell( theGrid, primary_cell_height, *pPath );

    if( definitely( theSet.covers( theCurrentCell.box() ) ) ) {
        pBinaryTreeNode->make_leaf(true);
        pBinaryTreeNode->mince( max_mince_depth - pPath->size() );
    } else if ( definitely( theSet.overlaps( theCurrentCell.box() ) ) ) {
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
    //Return to the previous level, since the initial evaluate is made
    //with the empty word, we check that it is not yet empty.
    if( pPath->size() > 0 ) {
        pPath->pop_back();
    }
}

Void GridTreePaving::adjoin_over_approximation( const ExactBoxType& theBoxType, const Nat numSubdivInDim ) {
    // FIXME: This adjoins an outer approximation; change to ensure only overlapping cells are adjoined
    for(SizeType i=0; i!=theBoxType.dimension(); ++i) {
        if(theBoxType[i].lower()>=theBoxType[i].upper()) {
            ARIADNE_THROW(std::runtime_error,"GridTreeSet::adjoin_over_approximation(ExactBoxType,Nat)","ExactBoxType "<<theBoxType<<" has empty interior.");
        }
    }
    this->adjoin_outer_approximation( theBoxType, numSubdivInDim );
}

Void GridTreePaving::adjoin_outer_approximation( const UpperBoxType& theBoxType, const Nat numSubdivInDim ) {
    ExactBoxSet theBoxSet=cast_exact_box(theBoxType);
    CompactSetInterface const& theSet=theBoxSet;
    this->adjoin_outer_approximation(theSet,numSubdivInDim);
}

Void GridTreePaving::adjoin_outer_approximation( const CompactSetInterface& theSet, const Nat numSubdivInDim ) {
    Grid theGrid( this->root_cell().grid() );
    ARIADNE_ASSERT( theSet.dimension() == this->root_cell().dimension() );

    //1. Compute the smallest GridCell (corresponding to the primary cell)
    //   that encloses the theSet (after it is mapped onto theGrid).
    const Nat height = GridCell::smallest_enclosing_primary_cell_height( theSet.bounding_box(), theGrid );
    //Compute the height of the primary cell for the outer approximation stepping up by the number of dimensions
    const Nat outer_approx_primary_cell_height = height + theGrid.dimension();

    //2. Align this paving and paving enclosing the provided set
    Bool has_stopped = false;
    BinaryTreeNode* pBinaryTreeNode = align_with_cell( outer_approx_primary_cell_height, true, false, has_stopped );

    //If the outer aproximation of the bounding box of the provided set is enclosed
    //in an enabled cell of this paving, then there is nothing to be done. The latter
    //is because adjoining the outer approx of the set will not change this paving.
    if( ! has_stopped ){
        //Compute the depth to which we must mince the outer approximation of the adjoining set.
        //This depth is relative to the root of the constructed paving, which has been alligned
        //with the binary tree node pBinaryTreeNode.
        const Nat max_mince_depth = zero_cell_subdivisions_to_tree_subdivisions( numSubdivInDim, outer_approx_primary_cell_height, 0 );

        //Adjoin the outer approximation, computing it on the fly.
        BinaryWord * pEmptyPath = new BinaryWord();
        _adjoin_outer_approximation( GridTreeSubpaving::_theGridCell.grid(), pBinaryTreeNode, outer_approx_primary_cell_height,
                                     max_mince_depth, theSet, pEmptyPath );

        delete pEmptyPath;
    }
}

Void GridTreePaving::adjoin_outer_approximation( const ValidatedCompactSetInterface& theSet, const Nat numSubdivInDim ) {
    Grid theGrid( this->root_cell().grid() );
    ARIADNE_ASSERT( theSet.dimension() == this->root_cell().dimension() );

    //1. Compute the smallest GridCell (corresponding to the primary cell)
    //   that encloses the theSet (after it is mapped onto theGrid).
    const Nat height = GridCell::smallest_enclosing_primary_cell_height( theSet.bounding_box(), theGrid );
    //Compute the height of the primary cell for the outer approximation stepping up by the number of dimensions
    const Nat outer_approx_primary_cell_height = height + theGrid.dimension();

    //2. Align this paving and paving enclosing the provided set
    Bool has_stopped = false;
    BinaryTreeNode* pBinaryTreeNode = align_with_cell( outer_approx_primary_cell_height, true, false, has_stopped );

    //If the outer aproximation of the bounding box of the provided set is enclosed
    //in an enabled cell of this paving, then there is nothing to be done. The latter
    //is because adjoining the outer approx of the set will not change this paving.
    if( ! has_stopped ){
        //Compute the depth to which we must mince the outer approximation of the adjoining set.
        //This depth is relative to the root of the constructed paving, which has been alligned
        //with the binary tree node pBinaryTreeNode.
        const Nat max_mince_depth = zero_cell_subdivisions_to_tree_subdivisions( numSubdivInDim, outer_approx_primary_cell_height, 0 );

        //Adjoin the outer approximation, computing it on the fly.
        BinaryWord * pEmptyPath = new BinaryWord();
        _adjoin_outer_approximation( GridTreeSubpaving::_theGridCell.grid(), pBinaryTreeNode, outer_approx_primary_cell_height,
                                     max_mince_depth, theSet, pEmptyPath );

        delete pEmptyPath;
    }
}

// TODO:Think of another representation in terms of covers but not pavings, then the implementation
// will be different, this is why, for now we do not fix these things.
Void GridTreePaving::adjoin_lower_approximation( const LocatedSetInterface& theSet, const Nat numSubdivInDim ) {
    this->adjoin_lower_approximation( theSet, cast_exact_box(theSet.bounding_box()), numSubdivInDim );
}

Void GridTreePaving::adjoin_lower_approximation( const OvertSetInterface& theSet, const Nat height, const Nat numSubdivInDim ) {
    Grid theGrid( this->root_cell().grid() );
    ARIADNE_ASSERT( theSet.dimension() == this->root_cell().dimension() );

    //Align this paving and paving at the given height
    Bool has_stopped = false;
    BinaryTreeNode* pBinaryTreeNode = align_with_cell( height, true, false, has_stopped );

    //If the lower aproximation of the bounding box of the provided set is enclosed
    //in an enabled cell of this paving, then there is nothing to be done. The latter
    //is because adjoining the outer approx of the set will not change this paving.
    if( ! has_stopped ){
        //Compute the depth to which we must mince the outer approximation of the adjoining set.
        //This depth is relative to the root of the constructed paving, which has been alligned
        //with the binary tree node pBinaryTreeNode.
        const Nat max_mince_depth = zero_cell_subdivisions_to_tree_subdivisions( numSubdivInDim, height, 0 );

        //Adjoin the outer approximation, computing it on the fly.
        BinaryWord * pEmptyPath = new BinaryWord();
        //const RegularSetInterface* theRegularVersionOfSet = dynamic_cast<const RegularSetInterface*>(&theSet);
        const OpenSetInterface* theOpenVersionOfSet = dynamic_cast<const OpenSetInterface*>(&theSet);
        //const LocatedSetInterface* theLocatedVersionOfSet = dynamic_cast<const LocatedSetInterface*>(&theSet);
        const OvertSetInterface* theOvertVersionOfSet = dynamic_cast<const OvertSetInterface*>(&theSet);
        if( theOpenVersionOfSet ) {
            _adjoin_lower_approximation( GridTreeSubpaving::_theGridCell.grid(), pBinaryTreeNode, height, max_mince_depth, *theOpenVersionOfSet, pEmptyPath );
        } else {
            _adjoin_lower_approximation( GridTreeSubpaving::_theGridCell.grid(), pBinaryTreeNode, height, max_mince_depth, *theOvertVersionOfSet, pEmptyPath );
        }
        delete pEmptyPath;
    }
}

Void GridTreePaving::adjoin_lower_approximation( const OvertSetInterface& theSet, const ExactBoxType& theBoundingBoxType, const Nat numSubdivInDim ) {
    Grid theGrid( this->root_cell().grid() );
    ARIADNE_ASSERT( theSet.dimension() == this->root_cell().dimension() );
    ARIADNE_ASSERT( theBoundingBoxType.dimension() == this->root_cell().dimension() );

    //Compute the smallest the primary cell that encloses the theSet (after it is mapped onto theGrid).
    const Nat height = GridCell::smallest_enclosing_primary_cell_height( theBoundingBoxType, theGrid );

    //Adjoin the lower approximation with the bounding cell being the primary cell at the given height.
    adjoin_lower_approximation( theSet, height, numSubdivInDim );
}

Void GridTreePaving::_adjoin_inner_approximation( const Grid & theGrid, BinaryTreeNode * pBinaryTreeNode, const Nat primary_cell_height,
                                               const Nat max_mince_depth, const OpenSetInterface& theSet, BinaryWord * pPath ) {
    //Compute the cell corresponding to the current node
    GridCell theCurrentCell( theGrid, primary_cell_height, *pPath );

    // If the set if closed, then we can use the separated(Box) method to speed-up computation
    const ClosedSetInterface* theClosedSet = dynamic_cast<const ClosedSetInterface*>(&theSet);

    if( ! pBinaryTreeNode->is_enabled() ) {
        //If this it is not an enabled leaf node then we can add something to it.
        if( definitely( theSet.covers( theCurrentCell.box() ) ) ) {
            //If this node's box is a subset of theSet then it belongs to the inner approximation
            //Thus we need to make it an enabled leaf and then do the mincing to the maximum depth
            pBinaryTreeNode->make_leaf( true );
        } else if ( theClosedSet && !definitely( theClosedSet->separated( theCurrentCell.box() ) ) ) {
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

    //Return to the previous level, since the initial evaluate is made
    //with the empty word, we check that it is not yet empty.
    if( pPath->size() > 0 ) {
        pPath->pop_back();
    }
}

Void GridTreePaving::adjoin_inner_approximation( const OpenSetInterface& theSet, const Nat height, const Nat numSubdivInDim ) {
    Grid theGrid( this->root_cell().grid() );
    ARIADNE_ASSERT( theSet.dimension() == this->root_cell().dimension() );

    //Align this paving and paving at the given height
    Bool has_stopped = false;
    BinaryTreeNode* pBinaryTreeNode = align_with_cell( height, true, false, has_stopped );

    //If the inner aproximation of the bounding box of the provided set is enclosed
    //in an enabled cell of this paving, then there is nothing to be done. The latter
    //is because adjoining the inner approx of the set will not change this paving.
    if( ! has_stopped ){
        //Compute the depth to which we must mince the inner approximation of the adjoining set.
        //This depth is relative to the root of the constructed paving, which has been alligned
        //with the binary tree node pBinaryTreeNode.
        const Nat max_mince_depth = zero_cell_subdivisions_to_tree_subdivisions( numSubdivInDim, height, 0 );

        //Adjoin the inner approximation, computing it on the fly.
        BinaryWord * pEmptyPath = new BinaryWord();
        _adjoin_inner_approximation( GridTreeSubpaving::_theGridCell.grid(), pBinaryTreeNode, height, max_mince_depth, theSet, pEmptyPath );
        delete pEmptyPath;
    }
}

Void GridTreePaving::adjoin_inner_approximation( const OpenSetInterface& theSet, const ExactBoxType& theBoundingBoxType, const Nat numSubdivInDim ) {
    Grid theGrid( this->root_cell().grid() );
    ARIADNE_ASSERT( theSet.dimension() == this->root_cell().dimension() );
    ARIADNE_ASSERT( theBoundingBoxType.dimension() == this->root_cell().dimension() );

    //Compute the smallest the primary cell that encloses the theSet (after it is mapped onto theGrid).
    const Nat height = GridCell::smallest_enclosing_primary_cell_height( theBoundingBoxType, theGrid );

    //Compute the height of the smallest primary cell that encloses the box
    //Note that, since we will need to adjoin inner approximation singleton by
    //the theBoundingBoxType, it is enough to take this cell's height and not to
    //go higher. Remember that the inner approximations consists of the cells
    //that are subsets of the set theSet Adjoin the inner approximation with
    //the bounding cell being the primary cell at the given height.
    adjoin_inner_approximation( theSet, height, numSubdivInDim );
}

Void GridTreePaving::adjoin_inner_approximation( const SetInterface& theSet, const Nat numSubdivInDim ) {
    OpenSetInterface const& theOpenSet = theSet;
    ExactBoxType theBoundingBox = cast_exact_box(theSet.bounding_box());
    this -> adjoin_inner_approximation(theOpenSet, theBoundingBox, numSubdivInDim);
}

Void GridTreePaving::adjoin_inner_approximation( const LowerBoxType& theBoxType, const Nat numSubdivInDim ) {
    ExactBoxSet theBoxSet=cast_exact_box(theBoxType);
    this->adjoin_inner_approximation(theBoxSet,theBoxSet,numSubdivInDim);
}

Void GridTreePaving::restrict_to_lower( const GridTreeSubpaving& theOtherSubPaving ){
    //The root of the binary tree of the current Paving
    BinaryTreeNode * pBinaryTreeNode = this->_pRootTreeNode;

    //The primary cell of this paving is higher then the one of the other paving.
    //1. We locate the path to the primary cell node common with the other paving
    BinaryWord rootNodePath = GridCell::primary_cell_path( this->root_cell().grid().dimension(), this->root_cell().height(), theOtherSubPaving.root_cell().height() );

    //2. Add the suffix path from the primary cell to the root node of
    //theOtherSubPaving. This is needed in order to be able to reach this root.
    rootNodePath.append( theOtherSubPaving.root_cell().word() );

    //3. Restrict this binary tree to the other one assuming the path prefix rootNodePath
    Nat position = 0;
    //This will point to the children nodes that do not have chance for being in the restriction
    BinaryTreeNode * pBranchToDisable = nullptr;
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

Void GridTreePaving::remove_from_lower( const GridTreeSubpaving& theOtherSubPaving ){
    //The root of the binary tree of the current Paving
    BinaryTreeNode * pBinaryTreeNode = this->_pRootTreeNode;

    //The primary cell of this paving is higher then the one of the other paving.
    //1. We locate the path to the primary cell node common with the other paving
    BinaryWord rootNodePath = GridCell::primary_cell_path( this->root_cell().grid().dimension(), this->root_cell().height(), theOtherSubPaving.root_cell().height() );

    //2. Add the suffix path from the primary cell to the root node of
    //theOtherSubPaving. This is needed in order to be able to reach this root.
    rootNodePath.append( theOtherSubPaving.root_cell().word() );

    //3. Remove theOtherSubPaving from this binary tree assuming the path prefix rootNodePath
    Nat position = 0;
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

Void GridTreePaving::clear(  ) {
    // TODO: Pieter: This implementation may not be optimal. Please could you
    // check this, Ivan.
    *this=GridTreePaving(this->grid());
}


Void GridTreePaving::restrict( const GridTreeSubpaving& theOtherSubPaving ) {
    const Nat thisPavingPCellHeight = this->root_cell().height();
    const Nat otherPavingPCellHeight = theOtherSubPaving.root_cell().height();

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

Void GridTreePaving::remove( const GridCell& theCell ) {
    ARIADNE_ASSERT( this->grid() == theCell.grid() );

    //If needed, extend the tree of this paving and then find it's the
    //primary cell common with the primary cell of the provided GridCell
    //If we encounter a disabled node then we do not move on, since then
    //there is nothing to remove, otherwise we split the node and go down
    Bool has_stopped = false;
    BinaryTreeNode* pCurrentPrimaryCell = align_with_cell( theCell.height(), false, true, has_stopped );

    if( ! has_stopped ) {
        //Follow theCell.word() path in the tree rooted to pCommonPrimaryCell,
        //do that until we encounter a leaf node, then stop
        BinaryWord path = theCell.word();
        Nat position = 0;
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

Void GridTreePaving::remove( const GridTreeSubpaving& theOtherSubPaving ) {
    const Nat thisPavingPCellHeight = this->root_cell().height();
    const Nat otherPavingPCellHeight = theOtherSubPaving.root_cell().height();

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

Void GridTreePaving::restrict_to_height( const Nat theHeight ) {
    const Nat thisPavingPCellHeight = this->root_cell().height();

    if( thisPavingPCellHeight > theHeight){
        // ARIADNE_WARN("restricting GridTreeSet of height " << this->root_cell().height() << " to height " << theHeight << ".\n");

        BinaryWord pathToPCell = GridCell::primary_cell_path( this->dimension(), thisPavingPCellHeight, theHeight );

        //Go throught the tree and disable all the leaves that
        //are not rooted to the primary cell defined by this path
        BinaryTreeNode * pCurrentNode = _pRootTreeNode;
        for( Nat i = 0; i < pathToPCell.size(); i++ ) {
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

//************************************FRIENDS OF GridTreeSet*****************************************/

GridTreePaving outer_approximation(const ExactBoxType& theBoxType, const Grid& theGrid, const Nat depth) {
    ExactBoxSet theBoxSet=theBoxType;
    CompactSetInterface const& theSet=theBoxSet;
    return outer_approximation(theSet,theGrid,depth);
}

GridTreePaving outer_approximation(const ExactBoxType& theBoxType, const Nat depth) {
    return outer_approximation(theBoxType, Grid(theBoxType.dimension()), depth);
}

GridTreePaving outer_approximation( const CompactSetInterface& theSet, const Grid& theGrid, const Nat numSubdivInDim ) {
    GridTreePaving result( theGrid );
    result.adjoin_outer_approximation( theSet, numSubdivInDim );
    return result;
}

GridTreePaving outer_approximation( const CompactSetInterface& theSet, const Nat numSubdivInDim ) {
    Grid theGrid( theSet.dimension() );
    return outer_approximation( theSet, theGrid, numSubdivInDim );
}

GridTreePaving outer_approximation( const ValidatedCompactSetInterface& theSet, const Grid& theGrid, const Nat numSubdivInDim ) {
    GridTreePaving result( theGrid );
    result.adjoin_outer_approximation( theSet, numSubdivInDim );
    return result;
}

GridTreePaving outer_approximation( const ValidatedCompactSetInterface& theSet, const Nat numSubdivInDim ) {
    Grid theGrid( theSet.dimension() );
    return outer_approximation( theSet, theGrid, numSubdivInDim );
}

GridTreePaving inner_approximation( const OpenSetInterface& theSet, const Grid& theGrid, const Nat height, const Nat numSubdivInDim ) {
    GridTreePaving result( theGrid );
    result.adjoin_inner_approximation( theSet, height, numSubdivInDim );
    return result;
}

GridTreePaving inner_approximation( const SetInterface& theSet, const Grid& theGrid, const Nat numSubdivInDim ) {
    GridTreePaving result( theGrid );
    result.adjoin_inner_approximation( theSet, numSubdivInDim );
    return result;
}

GridTreePaving inner_approximation( const OpenSetInterface& theSet, const Grid& theGrid, const ExactBoxType& bounding_box, const Nat numSubdivInDim ) {
    GridTreePaving result( theGrid );
    result.adjoin_inner_approximation( theSet, bounding_box, numSubdivInDim );
    return result;
}

GridTreePaving join( const GridTreeSubpaving& theSet1, const GridTreeSubpaving& theSet2 ) {
    //Test that the Grids are equal
    ARIADNE_ASSERT( theSet1.grid() == theSet2.grid() );

    //Compute the highest primary cell
    const Nat heightSet1 = theSet1.root_cell().height();
    const Nat heightSet2 = theSet2.root_cell().height();
    const Nat maxPCHeight = ( heightSet1 <  heightSet2 ) ? heightSet2 : heightSet1;

    //Create the resulting GridTreeSet
    GridTreePaving resultSet( theSet1.grid(), maxPCHeight, new BinaryTreeNode() );

    //Adjoin the sets
    resultSet.adjoin( theSet1 );
    resultSet.adjoin( theSet2 );

    //Return the result
    return resultSet;
}

GridTreePaving intersection( const GridTreeSubpaving& theSet1, const GridTreeSubpaving& theSet2 ) {
    //Test that the Grids are equal
    ARIADNE_ASSERT( theSet1.grid() == theSet2.grid() );

    //Compute the highest primary cell
    const Nat heightSet1 = theSet1.root_cell().height();
    const Nat heightSet2 = theSet2.root_cell().height();
    const Nat maxPCHeight = ( heightSet1 <  heightSet2 ) ? heightSet2 : heightSet1;

    //Create the resulting GridTreeSet
    GridTreePaving resultSet( theSet1.grid(), maxPCHeight, new BinaryTreeNode() );

    //Adjoin the first set
    resultSet.adjoin( theSet1 );
    //Intersect the result with the second set
    resultSet.restrict( theSet2 );

    //Return the result
    return resultSet;
}

GridTreePaving difference( const GridTreeSubpaving& theSet1, const GridTreeSubpaving& theSet2 ) {
    //Test that the Grids are equal
    ARIADNE_ASSERT( theSet1.grid() == theSet2.grid() );

    //Compute the highest primary cell
    const Nat heightSet1 = theSet1.root_cell().height();
    const Nat heightSet2 = theSet2.root_cell().height();
    const Nat maxPCHeight = ( heightSet1 <  heightSet2 ) ? heightSet2 : heightSet1;

    //Create the resulting GridTreeSet
    GridTreePaving resultSet( theSet1.grid(), maxPCHeight, new BinaryTreeNode() );

    //Adjoin the first set
    resultSet.adjoin( theSet1 );
    //Remove the second set from the result set
    resultSet.remove( theSet2 );

    //Return the result
    return resultSet;
}

Void draw(CanvasInterface& theGraphic, const Projection2d& theProjection, const GridCell& theGridCell) {
    theGridCell.box().draw(theGraphic,theProjection);
}

Void draw(CanvasInterface& theGraphic, const Projection2d& theProjection, const GridTreePaving& theGridTreeSet) {
    for(GridTreePaving::ConstIterator iter=theGridTreeSet.begin(); iter!=theGridTreeSet.end(); ++iter) {
        iter->box().draw(theGraphic,theProjection);
    }
}

Void draw(CanvasInterface& theGraphic, const Projection2d& theProjection, const CompactSetInterface& theSet) {
    static const Int DRAWING_DEPTH=16;
    draw(theGraphic,theProjection,outer_approximation(theSet,Grid(theSet.dimension()),DRAWING_DEPTH));
}

OutputStream& operator<<(OutputStream& os, const GridTreePaving& theGridTreeSet) {
    const GridTreeSubpaving& theGridTreeSubset = theGridTreeSet;
    return os << "GridTreeSet( " << theGridTreeSubset << " )";
}

} // namespace Ariadne

