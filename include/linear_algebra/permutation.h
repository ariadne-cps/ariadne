/***************************************************************************
 *            permutation.h
 *
 *  April 6, 2007, 3:04 PM
 *  Copyright  2004-7  Alberto Casagrande, Pieter Collins, Fons Kuijk
 *  casagrande@dimi.uniud.it Pieter.Collins@cwi.nl Fons.Kuijk@cwi.nl
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

/*! \file permutation.h
 *  \brief Permutation and Permutation operations.
 */

#ifndef ARIADNE_PERMUTATION_H
#define	ARIADNE_PERMUTATION_H

#include "base/array.h"

namespace Ariadne {
  
    
    /*! \brief A permutation. */
    class Permutation {
     private:
      array<uint> vec;
     public:
      /*! \brief Default constructor constructs permutation of zero values. */
      Permutation() : vec() { }
          
      /*! \brief Construct the identity permutation of \a n objects. */
      Permutation(const size_type& n) : vec(n) {
        for (size_type i=0; i<n; i++) { vec[i]=i; }
      }
          
      /*! \brief Construct from a beginning and end pointer. */
      Permutation(const uint* begin, const uint* end) : vec(end-begin) {
        for (size_type i=0; i<vec.size(); i++) { vec[i]=begin[i]; }
      }
      
      /*! \brief Destructor. */ 
      ~Permutation() { }
          
      /*!\brief A pointer to beginning of the permutation. */
      inline uint* begin() { return vec.begin(); }
          
      /*!\brief A pointer to end of the permutation. */
      inline uint* end() { return vec.end(); }
          
      /*!\brief The size of the permutation. */
          inline size_type size() { return vec.size(); }
          
      /*! \brief Swap the \a i<sup>th</sup> and \a j<sup>th</sup> elements. */
      inline void swap(unsigned int i, unsigned int j) {
        size_type sz = vec.size();
        if (i != j && i < sz && j < sz) { std::swap(vec[i], vec[j]); }
      }
          
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const {
        size_type sz = vec.size();
        if (sz > 0) {
          os << "[" << vec[0];
          for (uint i = 1; i < sz; i++)
            os << "," << vec[i];
          os << "]";
        }
        return os;
      }
      
      /*! \brief Return a reference to \f$\pi(i)\f$. */ 
      inline uint& operator[](const size_type& i) { return vec[i]; }
      
      /*! \brief Return \f$\pi(i)\f$.  */ 
      inline const uint& operator[](const size_type& i) const { return vec[i]; }
      
      /*! \brief Get the preimage \f$\pi^{-1}(j)\f$ of \a j, or \f$-1\f$ if \a j is not found. */
      inline int getindex(uint j) {
        size_type sz = vec.size();
        for(uint i = 0; i < sz; ++i) {
          if (vec[i] == j) { return i; }
        }
        return -1;
      }
      
#ifdef Doxygen
      /*!\brief Stream output operator. */
      friend std::ostream& operator<<(std::ostream& os, const Permutation& v);
#endif
    };
    
    inline std::ostream& operator<<(std::ostream& os, const Permutation& v) {
      return v.write(os);
    }
    
    
    
} // namespace Ariadne

#endif	/* ARIADNE_PERMUTATION_H */

