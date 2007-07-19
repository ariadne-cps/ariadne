/***************************************************************************
 *            affine_model.h
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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
 
/*! \file system/affine_model.h
 *  \brief An affine approximation to a function.
 */
 
#ifndef ARIADNE_AFFINE_MODEL_H
#define ARIADNE_AFFINE_MODEL_H

#include <vector>
#include <string>
#include <exception>
#include <stdexcept>

#include "../base/types.h"
#include "../base/array.h"
#include "../numeric/numerical_traits.h"
#include "../linear_algebra/declarations.h"

namespace Ariadne {
  namespace System {
      
    template<class R> class AffineModel;
    template<class R> AffineModel<R> operator+(const AffineModel<R>&, const AffineModel<R>&);
    template<class R> AffineModel<R> operator*(const AffineModel<R>&, const AffineModel<R>&);
    template<class R> AffineModel<R> inverse(const AffineModel<R>&);
    template<class R> AffineModel<R> implicit(const AffineModel<R>&);




    /*!\ingroup System
     * \ingroup DiscreteTime
     * \brief Concrete class for functions.
     */
    template<class R>
    class AffineModel
    {
      typedef typename Numeric::traits<R>::arithmetic_type A; 
      typedef typename Numeric::traits<R>::interval_type I; 
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      
      /*! \brief Destructor. */
      ~AffineModel();
     
      /*! \brief Constructor. */
      AffineModel(dimension_type rd, dimension_type ad, size_type s);
     
      /*! \brief Set an element in degree 0 to x. */
      template<class X> set(dimension_type i, const X& x);

      /*! \brief Set an element in degree 1 to x. */
      set(dimension_type i, dimension_type j, const R& x);
      /*! \brief Set an element in degree 1 to x. */
      set(dimension_type i, dimension_type j, const I& x);

     private:
      dimension_type _result_dimension;
      dimension_type _argument_dimension;
      static const size_type _degree=1;
      size_type _smoothness;
      I* _intervals_begin;
      R* _reals_begin;
    };
   

    template<class R> inline
    AffineModel<R>::instantiate()
    {
      AffineModel<R>* am=0;
      *am + *am;
      *am * *am;
      compose(*am,*am);
      inverse(*am);
      implicit(*am);
    }


    template<class R> inline
    AffineModel<R>::~AffineModel() 
    {
      delete[] _intervals_begin;
    }


    template<class R> inline
    AffineModel<R>::AffineModel(dimension_type rd, dimension_type ad, size_type s) 
      : _result_dimension(rd), _argument_dimension(ad), _smoothness(s),
    {
      if(s=0) {
        R* ptr=new R(rd*(ad+2));
        _intervals_begin=reinterpret_cast<I*>(ptr);
        _reals_begin=ptr+rd;
      } else if(s=1) {
        R* ptr=new R(2*rd*(ad+1));
        _intervals_begin=reinterpret_cast<I*>(ptr);
        _reals_begin=ptr+(rd*(ad+1));
      } else {
        throw std::runtime_error("An affine model must have smoothness 0 or 1");
      }
    }


    template<class R> template<class X> inline
    void AffineModel<R>::set(dimension_type i, const X& x) 
    {
      assert(i<this->_rd);
      this->_rptr[i]=x;
    }

    template<class R> inline
    void AffineModel<R>::set(dimension_type i, dimension_type j, const R& x) 
    {
      assert(i<this->_rd);
      assert(j<this->_ad);
      if(this->_s==0) {
        this->_rptr[i+this->_rd*(j+1)]=x;
      } else {
        this->_iptr[i+this->_rd*(j+1)]=x;
      }
    }

    template<class R> inline
    void AffineModel<R>::set(dimension_type i, dimension_type j, const I& x) 
    {
      assert(i<this->_rd);
      assert(j<this->_ad);
      if(this->_s==0) {
        throw std::runtime_error("Error: Attempting to assign an interval value to a real.");
      } else {
        this->_iptr[i+this->_rd*(j+1)]=x;
      }
    }


   
    template<class R> inline
    AffineModel<R>::write(std::ostream& os) 
    {
      os << "AffineModel( ";
      for(uint i=0; i!=this->_rd; ++i) {
        os << (i==0 ? '[' : ',') << this->_iptr[i]; 
      }
      os << "], ";
      if(s=0) {
        for(uint i=0; i!=this->_rd; ++i) {
          for(uint j=0; j!=this->_ad; ++j) {
            os << ( i==0 ? '[' : (j==0 ? ';' : ',') )
               << this->_iptr[i+(j+1)*+this->_rd]; 
          }
        }
        os << ']';
      } else {
        for(uint i=0; i!=this->_rd; ++i) {
          for(uint j=0; j!=this->_ad; ++j) {
            os << ( i==0 ? '[' : (j==0 ? ';' : ',') )
               << this->_rptr[i+(j+1)*+this->_rd]; 
          }
        }
        os << ']';
      }
      os << " ]";
    }

    template<class R> inline
    AffineModel<R>
    operator+(const AffineModel<R>& am1, const AffineModel<R>& am2) 
    {
      assert(am1._rd==am2._rd && am1._ad==am2._ad);
      dimension_type rd=am1._rd;
      dimension_type ad=am1._rd;
      size_type s=std::min(am1._s,am2._s);
      AffineModel<R> res(rd,ad,s);
      if(res.smoothness==1) {
        for(size_type i=0; i!=rd*(sd+1); ++i) {
          res._iptr[i]=am1._iptr[i]+am2._iptr[i];
        }
      } else {
        if(std::max(am1._s,am2._s)==0) {
          for(size_type i=0; i!=rd; ++i) {
            for(size_type j=0; i!=ad; ++i) {
              res._rptr[i]=am1._iptr[i]+am2._iptr[i];
            }
          }
        } else {
        }
      }
      return result;
    }

    template<class R> inline
    std::istream& operator>>(std::istream& is, AffineModel<R>& am) {
      return f.read(is);
    }

    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const AffineModel<R>& am) {
      return f.write(os);
    };



  }
}

#endif /* ARIADNE_AFFINE_MODEL_H */
