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

namespace Ariadne {
  namespace LinearAlgebra {
    
    class Permutation {
      
      private:
        std::vector<int> vec;
        
        public:
          /*! \brief Construct  permutation. */
          Permutation() : vec() { }
          
          /*! \brief Construct the default permutation of size \a n. */
          Permutation(const size_type& n) : vec(n, 0) {
            for (size_type i=0; i<n; i++)
              vec[i]=i;
          }
          
          ~Permutation() { }
          
          /*!\brief Return size of Permutation. */
          inline size_type size() { return vec.size(); }
          
          /*! \brief Swap two elements. */
          inline void swap(unsigned int a, unsigned int b) {
            size_type sz = vec.size();
            if (a != b && a < sz && b < sz)
              std::swap(vec[a], vec[b]);
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
          
          inline  int& operator[](const size_type& i) { return vec[i]; }
          
          inline int getindex(int index) {
            size_type sz = vec.size();
            if (uint(index) > sz) return -1;
            
            if (sz > 0)
              for (uint i = 0; i < sz; i++)
                if (vec[i] == index) return i;
            return -1;
          }
          
    };
    
    inline std::ostream& operator<<(std::ostream& os, const Permutation& v) {
      return v.write(os);
    }
    
    
    
  } // namespace LinearAlgebra
} // namespace Ariadne

#endif	/* ARIADNE_PERMUTATION_H */

