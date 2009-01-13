/***************************************************************************
 *            binary_word.h
 *
 *  Copyright  2004-8  Alberto Casagrande, Pieter Collins, Ivan S. Zapreev
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
 
/*! \file binary_word.h
 *  \brief Binary words and sets of binary words.
 *
 * Binary words are generally useful objects. 
 * This file provides a class of binary words of 
 * arbitrary size, complementing the STL bitset and
 * Vector<bool> classes. The main purpose of giving 
 * this class is to have a known implementation of 
 * binary words which may help in giving efficient
 * implementations of the list classes described below.
 * The most efficient implementation may, in fact, be an
 * STL Vector<bool>. Hence another purpose of this class
 * is to provide a known interface for binary words.
 *
 * We give a number of classes representing sets of binary words. 
 * These classes have been designed to be optimised for storage, and for 
 * sequential iteration through the elements. Hence the preferred method of
 * access is through iterators, which need only conform to the requirements
 * of an InputIterator. These  classes may be immutable. 
 * The fundamental operations are to join two sets,
 * and to convert between different representations.
 * 
 * A special class is given representing prefix-free sets of words
 * using a binary tree representation. 
 *
 * The main intended use of these classes in Ariadne is to represent 
 * sets as unions of rectangles on a grid. See the file grid_box.h
 * for the implementation of this representation. 
 * The basic idea is that a inside of Euclidean space contained in some 
 * cuboid R can be represented by a union of sets obtained by repeatedly
 * subdividing R in two along the coordinate axes. 
 * The nth subdivision occurs along a coordinate specified by the Grid 
 * class, and whether the lower or upper set is taken depends the value 
 * of the nth element of a BinaryWord. 
 * In this way, sets of points can be represented efficiently by just listing
 * in which half of the remaining rectangle the subdivision occurs in.
 *
 * NOTE TO ALBERTO : This describes the basic ideas. Hope it's comprehensible. 
 * I haven't finished the implementation of all features. There are loose ends
 * in how we choose subdivision directions. I hope all this effort will result
 * in a better implementation than "pointer-based" quad trees; it should be
 * much more memory-efficient, but may be unacceptably slow. On the other hand,
 * for big grids or fine subdivisions, memory may well be a limiting factor.
 */

#ifndef ARIADNE_BINARY_WORD_H
#define ARIADNE_BINARY_WORD_H

#include <limits>
#include <vector>
#include <iosfwd>
#include <stdexcept>

#include "macros.h"
#include "stlio.h"

namespace Ariadne {

/*! \brief A statically-allocated binary word of fixed maximum length.
 */
class BinaryWord : public std::vector<bool> {
  private:
    typedef unsigned char byte_type;
    // The number of bits per byte
    static const size_type _bits_per_byte=std::numeric_limits<byte_type>::digits;
  public:
    typedef std::vector<bool>::const_iterator const_iterator;

    //! \brief Default constructor makes an empty word. 
    BinaryWord() : std::vector<bool>() { }

    //! Comparison operator. 
    bool operator<(const BinaryWord& w) const;

    //! \brief The number of machine bytes the word takes up. 
    size_type bytes() const { return (size()+_bits_per_byte-1)/_bits_per_byte; }

    //! \brief Sets the last bit to \a x.  
    void set_back (bool x) { (*this)[this->size()-1]=x; }
            
    //! \brief Erases the prefix of the word of size \a prefix_size
    void erase_prefix(const int prefix_size);
            
    //! \brief true if this word is a prefix of the other word \a b.
    bool is_prefix(const BinaryWord& b) const;

    //! \brief true if the word is a subword of the other word. 
    bool is_subword(const BinaryWord& b) const;

    //! \brief appends the binary word \a binaryWord
    void append( const BinaryWord& binaryWord );
};

inline void BinaryWord::erase_prefix(const int prefix_size) {
    const int total_size = size();
    //Erase prefix elements until the word is empry or we erased prefix_size elements
    for( int i = 0; ( i < total_size ) && ( i < prefix_size ) ; i++ ) {
        erase( begin() );
    }
}
    
inline bool  BinaryWord::operator<(const BinaryWord& w) const { 
    for( size_type i = 0; i != std::min( this->size(),w.size() ); ++i ) {
        if( (*this)[i] != w[i] ) {
            return (*this)[i] < w[i];
        }
    }
    return this->size() < w.size(); 
}

inline bool BinaryWord::is_prefix( const BinaryWord& b ) const {
    if( this->size() > b.size() ) {
        return false;
    }
    for( size_type i=0; i != this->size(); ++i ) {
        if( (*this)[i] != b[i] ) {
            return false;
        }
    }
    return true;
}

inline bool BinaryWord::is_subword(const BinaryWord& b) const {
    if( this->size() > b.size() ) {
        return false;
    }
    for( size_type i = 0; i != ( b.size() - this->size() + 1 ); ++i ) { 
        size_type j=0;
        while( j != this->size() && (*this)[j] == b[i+j] ) {
            ++j;
        }
        if( j == this->size() ) {
            return true;
        }
    }
    return false;
}

inline void BinaryWord::append( const BinaryWord& binaryWord ) {
    for( size_t i = 0; i < binaryWord.size() ; i++ ){
        this->push_back( binaryWord[i] );
    }
}

/*! \brief Provides a method to fill a BinaryWord with a new data form the stream
 *  The data format is either "[0,1,0,1,1]" or equivalently "01011".
 */
inline std::istream& operator>>(std::istream& input_stream, BinaryWord& binary_word) {
    //Clear the binary word, as we assume that we input
    //new data and the old one is not needed.
    binary_word.clear();
    //Read the first stream symbol to check the input format
    char symbol;
    input_stream >> symbol;
    if( symbol == '[' ) {
        //If it is a standard std:vector input then put
        //the symbol back and use it's input operator
        std::vector<bool> & bool_vector = binary_word;
        input_stream.putback( symbol );
        input_stream >> bool_vector;
    } else {
        //Otherwise read a sequence of 0 and 1.
        while( input_stream && ( symbol =='0' || symbol =='1' ) ) {
            binary_word.push_back( symbol == '0' ? 0 : 1);
            input_stream.get( symbol );
        }
        input_stream.putback( symbol );
    }
    return input_stream;
}

/*! \brief allocates a new binary word created from a string.
 *  The data format is either "[0,1,0,1,1]" or equivalently "01011".
 */
inline BinaryWord make_binary_word(const std::string& string_data) {
    BinaryWord binary_word;
    std::stringstream string_input_stream( string_data );
    string_input_stream >> binary_word;
    return binary_word;
}
    
/*! \brief Serializes data into a stream.
 * The data format is "01011" or "e" for an empty word.
 */
inline std::ostream& operator<<(std::ostream& os, const BinaryWord& bw) {
    if( bw.empty() ) {
        return os << 'e';
    }
    for(size_t i=0; i!=bw.size(); ++i) {
        os << bw[i];
    }
    return os;
}

} // namespace Ariadne
  
#endif /* ARIADNE_BINARY_WORD_H */
