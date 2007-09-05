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
#include "../output/logging.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

namespace Ariadne {
  namespace Function {
      
    template<class R> class AffineModel;
    template<class R> AffineModel<R> operator+(const AffineModel<R>&, const AffineModel<R>&);
    template<class R> AffineModel<R> operator*(const AffineModel<R>&, const AffineModel<R>&);
    template<class R> AffineModel<R> reduce(const AffineModel<R>&, size_type);
    template<class R> AffineModel<R> compose(const AffineModel<R>&, const AffineModel<R>&);
    template<class R> AffineModel<R> inverse(const AffineModel<R>&);
    template<class R> AffineModel<R> implicit(const AffineModel<R>&);


    /*!\ingroup FunctionModel
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
      AffineModel();
     
      /*! \brief Constructor. */
      AffineModel(dimension_type rd, dimension_type ad, size_type s);
     
      /*! \brief Constructor. */
      AffineModel(const LinearAlgebra::Vector<I>& a0, const LinearAlgebra::Matrix<R>& a1);
     
      /*! \brief Constructor. */
      AffineModel(const LinearAlgebra::Vector<I>& a0, const LinearAlgebra::Matrix<I>& a1);
     
      /*! \brief Copy constructor. */
      AffineModel(const AffineModel<R>& am);
     
      /*! \brief Assignement operator. */
      AffineModel<R>& operator=(const AffineModel<R>& am);

      /*! \brief Set an element in degree 0 to x. */
      template<class X> void set(dimension_type i, const X& x);

      /*! \brief Set an element in degree 1 to x. */
      void set(dimension_type i, dimension_type j, const R& x);
      /*! \brief Set an element in degree 1 to x. */
      void set(dimension_type i, dimension_type j, const I& x);

      /*! \brief The smoothness of the model. */
      size_type smoothness() const;

      /*! \brief The centre of the model. */
      LinearAlgebra::VectorSlice<const I> centre() const;

      /*! \brief The jacobian derivative of the model. (Smoothness at least 1). */
      LinearAlgebra::MatrixSlice<const I> jacobian() const;

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;

      template<class X> friend AffineModel<X> reduce(const AffineModel<X>&, size_type); 
      template<class X> friend AffineModel<X> inverse(const AffineModel<X>&, size_type); 
      template<class X> friend AffineModel<X> implicit(const AffineModel<X>&, size_type); 
      template<class X> friend AffineModel<X> operator+(const AffineModel<X>&, const AffineModel<X>&); 
      template<class X> friend AffineModel<X> operator*(const AffineModel<X>&, const AffineModel<X>&); 
      template<class X> friend AffineModel<X> compose(const AffineModel<X>&, const AffineModel<X>&);

     private:
      static void instantiate();
      static size_type real_array_size(dimension_type rd, dimension_type ad, size_type s);
      size_type real_array_size() const;
      void allocate(dimension_type rd, dimension_type ad, size_type s);
      void reallocate(dimension_type rd, dimension_type ad, size_type s);
      void deallocate();
      void assign(const AffineModel<R>& am);
      LinearAlgebra::MatrixSlice<const R> approximate_jacobian() const;

     private:
      dimension_type _rd;
      dimension_type _ad;
      static const size_type _d=1;
      size_type _s;
      I* _iptr;
      R* _rptr;
    };
   

    template<class R> inline
    void AffineModel<R>::instantiate()
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
      ARIADNE_LOG(3,"AffineModel<R>::~AffineModel()\n");
      this->deallocate();
    }


    template<class R> inline
    AffineModel<R>::AffineModel() 
      : _rd(0), _ad(0), _s(0)
    {
      R* ptr=new R(0);
      this->_iptr=reinterpret_cast<I*>(ptr);
      this->_rptr=ptr;
    }

    template<class R> inline
    AffineModel<R>::AffineModel(dimension_type rd, dimension_type ad, size_type s) 
    {
      ARIADNE_LOG(3,"AffineModel<R>::AffineModel(dimension_type rd, dimension_type ad, size_type s)\n");
      this->allocate(rd,ad,s);
    }


    template<class R> inline
    AffineModel<R>::AffineModel(const LinearAlgebra::Vector<I>& a0, const LinearAlgebra::Matrix<R>& a1) 
      : _rd(a0.size()), _ad(a1.number_of_columns()), _s(0)
    {
      ARIADNE_LOG(3,"AffineModel<R>::AffineModel(Vector<I> a0, Matrix<R> a1)\n");
      R* ptr = new R[_rd*(2+_ad)];
      this->_iptr = reinterpret_cast<I*>(ptr);
      this->_rptr = ptr+_rd;

      LinearAlgebra::VectorSlice<I>(this->_rd,this->_iptr,1u) = a0;
      LinearAlgebra::MatrixSlice<R>(this->_rd,this->_ad,this->_rptr+this->_rd,1u,this->_rd) = a1;
    }


    template<class R> inline
    AffineModel<R>::AffineModel(const LinearAlgebra::Vector<I>& a0, const LinearAlgebra::Matrix<I>& a1) 
      : _rd(a0.size()), _ad(a1.number_of_columns()), _s(1)
    {
      ARIADNE_LOG(3,"AffineModel<R>::AffineModel(Vector<I> a0, Matrix<I> a1)\n");
      R* ptr = new R[2*_rd*(1+_ad)];
      this->_iptr = reinterpret_cast<I*>(ptr);
      this->_rptr = ptr+_rd*(1+_ad);

      LinearAlgebra::VectorSlice<I>(this->_rd,this->_iptr,1) = a0;
      LinearAlgebra::MatrixSlice<I>(this->_rd,this->_ad,this->_iptr+this->_rd,1,this->_rd) = a1;
    }


    template<class R> inline
    AffineModel<R>::AffineModel(const AffineModel<R>& am) 
    {
      ARIADNE_LOG(3,"AffineModel<R>::AffineModel(AffineModel<R> am)\n");
      this->allocate(am._rd,am._ad,am._s);
      this->assign(am);
    }


    template<class R> inline
    AffineModel<R>& AffineModel<R>::operator=(const AffineModel<R>& am) 
    {
      ARIADNE_LOG(3,"AffineModel<R>::operator=(AffineModel<R> am)\n");
      if(this!=&am) {
        this->reallocate(am._rd,am._ad,am._s);
        this->assign(am);
      }
      return *this;
    }


    template<class R> inline
    void AffineModel<R>::allocate(dimension_type rd, dimension_type ad, size_type s)
    {
      assert(s==0 or s==1);
      this->_rd=rd;
      this->_ad=ad;
      this->_s=s;

      size_type n=2*_rd+(1+_s)*_rd*_ad;
      R* ptr=new R[n];
      if(_s==0) {
        this->_iptr=reinterpret_cast<I*>(ptr);
        this->_rptr=ptr+_rd;
      } else {
        this->_iptr=reinterpret_cast<I*>(ptr);
        this->_rptr=ptr+(_rd*(_ad+1));
      } 
    }


    template<class R> inline
    void AffineModel<R>::reallocate(dimension_type rd, dimension_type ad, size_type s)
    {
      if(real_array_size(this->_rd,this->_ad,this->_s) != real_array_size(rd,ad,s)) {
        this->deallocate();
        this->allocate(rd,ad,s);
      }
    }


    template<class R> inline
    void AffineModel<R>::deallocate()
    {
      R* ptr=reinterpret_cast<R*>(this->_iptr);
      delete[] ptr;
    }


    template<class R> inline
    size_type AffineModel<R>::real_array_size(dimension_type rd, dimension_type ad, size_type s) 
    {
      return 2*rd+(1+s)*rd*ad;
    }


    template<class R> inline
    size_type AffineModel<R>::real_array_size() const
    {
      return 2*this->_rd+(1+this->_s)*this->_rd*this->_ad;
    }


    template<class R> inline
    void AffineModel<R>::assign(const AffineModel<R>& am)
    {
      R* ptr=reinterpret_cast<R*>(this->_iptr);
      const R* amptr=reinterpret_cast<R*>(am._iptr);
      for(size_type i=0; i!=2*_rd+(1+_s)*_rd*_ad; ++i) {
        ptr[i]=amptr[i];
      } 
    }


    template<class R> inline
    size_type AffineModel<R>::smoothness()  const
    {
      return this->_s;
    }


    template<class R> inline
    LinearAlgebra::VectorSlice<const typename AffineModel<R>::I> 
    AffineModel<R>::centre() const
    {
      return LinearAlgebra::VectorSlice<const I>(this->_rd,this->_iptr,1u);
    }

    template<class R> inline
    LinearAlgebra::MatrixSlice<const typename AffineModel<R>::I> 
    AffineModel<R>::jacobian() const
    {
      assert(this->_s>=1);
      return LinearAlgebra::MatrixSlice<const I>(this->_rd,this->_ad,this->_iptr+this->_rd,1u,this->_rd);
    }

    template<class R> inline
    LinearAlgebra::MatrixSlice<const R> 
    AffineModel<R>::approximate_jacobian() const
    {
      assert(this->_s<1);
      return LinearAlgebra::MatrixSlice<const R>(this->_rd,this->_ad,this->_rptr+this->_rd,1u,this->_rd);
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
    AffineModel<R>
    operator+(const AffineModel<R>& am1, const AffineModel<R>& am2) 
    {
      ARIADNE_LOG(3,"operator+(AffineModel am1, AffineModel am2)\n");
      assert(am1._rd==am2._rd && am1._ad==am2._ad);
      Numeric::Interval<R> tmp;
      dimension_type rd=am1._rd;
      dimension_type ad=am1._ad;
      size_type s=std::min(am1._s,am2._s);
      AffineModel<R> res(rd,ad,s);
      for(size_type i=0; i!=rd; ++i) {
        res._iptr[i]=am1._iptr[i]+am2._iptr[i];
      }
      if(res.smoothness()==1) {
        for(size_type i=rd; i!=rd*(ad+1u); ++i) {
          res._iptr[i]=am1._iptr[i]+am2._iptr[i];
        }
      } else {
        if(std::max(am1._s,am2._s)==0) {
          for(size_type i=0; i!=rd; ++i) {
            for(size_type j=i; j!=i+rd*ad; j+=rd) {
              tmp=am1._rptr[j]+am2._rptr[j];
              res._rptr[j]=tmp.midpoint();
              res._iptr[i]+=(tmp-tmp.midpoint());
            }
          }
        } else {
          assert(false);
        }
      }
      return res;
    }



    template<class R> inline
    AffineModel<R>
    compose(const AffineModel<R>& am1, const AffineModel<R>& am2) 
    {
      ARIADNE_LOG(3,"AffineModel compose(AffineModel am1, AffineModel am2)\n");
      //FIXME: Use slices
      typedef typename AffineModel<R>::I I;
      using namespace LinearAlgebra;
      assert(am1._ad==am2._rd);
      Vector<I> b0;
      Matrix<I> A0;
      Vector< I> b1=am1.centre();
      Vector< I> b2=am2.centre();
      if(am1._s==0) {
        Matrix< R> A1=am1.approximate_jacobian();
        b0=b1+A1*b2;
        if(am2._s==0) {
          Matrix< R> A2=am2.approximate_jacobian();
          A0=A1*A2;
        } else {
          Matrix< I> A2=am2.jacobian();
          A0=A1*A2;
        }
      } else {
        Matrix< I> A1=am1.jacobian();
        b0=b1+A1*b2;
        if(am2._s==0) {
          Matrix< R> A2=am2.approximate_jacobian();
          A0=A1*A2;
        } else {
          Matrix< I> A2=am2.jacobian();
          A0=A1*A2;
        }
      }
      AffineModel<R> res(b0,A0);
      return reduce(res,std::min(am1._s,am2._s));
    }

        
        

    template<class R> inline
    AffineModel<R>
    reduce(const AffineModel<R>& am, size_type s) 
    {
      assert(s<=am._s);
      if(s==am._s) { return am;  }
      
      const dimension_type rd=am._rd;
      const dimension_type ad=am._ad;
      
      AffineModel<R> res(rd,ad,s);
      typename AffineModel<R>::I tmp;
      for(size_type i=0; i!=rd; ++i) {
        res._iptr[i]=am._iptr[i];
        for(size_type j=i+rd; j!=i+rd*(ad+1); j+=rd) {
          res._rptr[j]=am._iptr[j].midpoint();
          res._iptr[i]+=(am._iptr[j]-res._rptr[j]);
        }
      }
      return res;
    }
          


    template<class R> inline
    AffineModel<R>
    operator*(const AffineModel<R>& am1, const AffineModel<R>& am2) 
    {
      ARIADNE_LOG(3,"operator+(AffineModel am1, AffineModel am2)\n");
      assert(am1._rd==am2._rd && am1._ad==am2._ad);
      Numeric::Interval<R> tmp;
      dimension_type rd=am1._rd;
      dimension_type ad=am1._ad;
      size_type s=std::min(am1._s,am2._s);
      AffineModel<R> res(rd,ad,s);
      for(size_type i=0; i!=rd; ++i) {
        res._iptr[i]=am1._iptr[i]+am2._iptr[i];
      }
      if(res.smoothness()==1) {
        for(size_type i=rd; i!=rd*(ad+1u); ++i) {
          res._iptr[i]=am1._iptr[i]+am2._iptr[i];
        }
      } else {
        if(std::max(am1._s,am2._s)==0) {
          for(size_type i=0; i!=rd; ++i) {
            for(size_type j=i; j!=i+rd*ad; j+=rd) {
              tmp=am1._rptr[j]+am2._rptr[j];
              res._rptr[j]=tmp.midpoint();
              res._iptr[i]+=(tmp-tmp.midpoint());
            }
          }
        } else {
          assert(false);
        }
      }
      return res;
    }

    template<class R> inline
    AffineModel<R> 
    inverse(const AffineModel<R>&) 
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R> inline
    AffineModel<R> 
    implicit(const AffineModel<R>&)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R> inline
    std::ostream& AffineModel<R>::write(std::ostream& os) const
    {
      os << "AffineModel( a0=" << LinearAlgebra::VectorSlice<I>(this->_rd,this->_iptr,1) << ", a1=";
      if(this->_s==0) {
        os  << LinearAlgebra::MatrixSlice<R>(this->_rd,this->_ad,this->_rptr+this->_rd,1,this->_rd);
      } else {
        os << LinearAlgebra::MatrixSlice<I>(this->_rd,this->_ad,this->_iptr+this->_rd,1,this->_rd);
      }
      return os << " )";
    }


  /*
  template<class R> inline
    std::istream& operator>>(std::istream& is, AffineModel<R>& am) {
      return am.read(is);
    }
  */

    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const AffineModel<R>& am) {
      return am.write(os);
    };



  }
}

#endif /* ARIADNE_AFFINE_MODEL_H */
